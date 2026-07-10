#!/usr/bin/env python
"""scVI / sysVI count-based integration for the seurat-pipeline step 04.

Corrects over an arbitrary batch column (--batch-key), e.g. chemistry, donor, or
sequencing run. Especially suited to large cross-technology gaps (e.g. a probe panel
vs whole-transcriptome) that linear methods can't bridge.

Plain-file contract (no h5ad round-trip): reads a Matrix Market counts file plus
gene / barcode / metadata sidecars written by the R side, trains a count-based
deep integration model with the batch key held out, and writes the batch-corrected
latent representation back as a CSV keyed by barcode. The R helper aligns that CSV
to the Seurat object by barcode and installs it as a DimReduc.

Input counts orientation: genes (rows) x cells (cols) — i.e. Matrix::writeMM of a
Seurat counts layer. Transposed to cells x genes for AnnData here.
"""
import argparse
import json
import sys

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.sparse as sp


def log(msg):
    print(f"[scvi_integrate] {msg}", file=sys.stderr, flush=True)


def parse_args():
    p = argparse.ArgumentParser(description="scVI/sysVI integration bridge")
    p.add_argument("--counts", required=True, help="Matrix Market file, genes x cells")
    p.add_argument("--genes", required=True, help="gene names, one per line (rows of counts)")
    p.add_argument("--barcodes", required=True, help="cell barcodes, one per line (cols of counts)")
    p.add_argument("--metadata", required=True, help="CSV indexed by barcode; must hold batch-key column")
    p.add_argument("--out", required=True, help="output latent CSV (barcode + Z1..Zn)")
    p.add_argument("--summary", default=None, help="optional JSON run summary out")
    p.add_argument("--model", default="scvi", choices=["scvi", "sysvi"])
    p.add_argument("--batch-key", default="batch",
                   help="metadata column to correct over (e.g. chemistry, donor, run)")
    p.add_argument("--n-latent", type=int, default=30)
    p.add_argument("--n-layers", type=int, default=2)
    p.add_argument("--gene-likelihood", default="nb", choices=["nb", "zinb", "poisson"])
    p.add_argument("--dispersion", default="gene", choices=["gene", "gene-batch", "gene-cell"])
    p.add_argument("--n-hvg", type=int, default=2000)
    p.add_argument("--min-cells-per-batch", type=int, default=3,
                   help="shared gene space: keep genes detected in >= this many cells in EVERY batch")
    p.add_argument("--include-genes", default="", help="comma-separated genes to force into the HVG set")
    p.add_argument("--max-epochs", type=int, default=None, help="null/omit -> scvi cell-count heuristic")
    p.add_argument("--accelerator", default="cpu", help="cpu (default, deterministic) | mps | gpu | auto")
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


def load_anndata(args):
    log(f"reading counts {args.counts}")
    mat = scipy.io.mmread(args.counts)          # genes x cells (sparse COO)
    mat = sp.csr_matrix(mat).T.tocsr()          # -> cells x genes
    genes = [l.strip() for l in open(args.genes)]
    barcodes = [l.strip() for l in open(args.barcodes)]
    if mat.shape != (len(barcodes), len(genes)):
        raise SystemExit(
            f"shape mismatch: counts is {mat.shape} (cells x genes) but got "
            f"{len(barcodes)} barcodes and {len(genes)} genes"
        )
    meta = pd.read_csv(args.metadata, index_col=0)
    meta = meta.reindex(barcodes)               # align metadata rows to matrix cols
    if meta.isnull().any().any():
        missing = int(meta.isnull().any(axis=1).sum())
        raise SystemExit(f"{missing} barcodes missing from metadata (index must match barcodes.txt)")

    adata = sc.AnnData(X=mat, obs=meta.copy())
    adata.obs_names = barcodes
    adata.var_names = genes
    adata.layers["counts"] = adata.X.copy()     # scVI trains on raw counts from here
    if args.batch_key not in adata.obs.columns:
        raise SystemExit(f"batch key '{args.batch_key}' not in metadata columns {list(adata.obs.columns)}")
    adata.obs[args.batch_key] = adata.obs[args.batch_key].astype("category")
    log(f"loaded {adata.n_obs} cells x {adata.n_vars} genes; "
        f"batches: {dict(adata.obs[args.batch_key].value_counts())}")
    return adata


def filter_shared_genes(adata, args):
    """Restrict to the gene space shared across batches: genes detected in >= min_cells
    in EVERY batch category. Matters when batches don't measure the same genes — e.g. a
    probe panel vs whole-transcriptome (Flex vs 3'), where panel-absent genes are
    structurally all-zero in one batch and become perfect batch labels. Such genes must
    be dropped before HVG selection so the model doesn't 'correct' a signal we handed it
    for free. Harmless (a no-op) when all batches share the same gene space."""
    bk = args.batch_key
    detected = (adata.layers["counts"] > 0)
    keep = np.ones(adata.n_vars, dtype=bool)
    per_batch = {}
    for b in adata.obs[bk].cat.categories:
        rows = (adata.obs[bk] == b).to_numpy()
        n_cells = np.asarray(detected[rows].sum(axis=0)).ravel()
        per_batch[str(b)] = int((n_cells >= args.min_cells_per_batch).sum())
        keep &= (n_cells >= args.min_cells_per_batch)
    log(f"shared gene space (>= {args.min_cells_per_batch} cells/batch): "
        f"{int(keep.sum())} of {adata.n_vars} genes "
        f"(per-batch detected: {per_batch})")
    if keep.sum() == 0:
        raise SystemExit("no genes shared across batches at this threshold")
    return adata[:, keep].copy(), {"n_genes_union": int(adata.n_vars),
                                   "n_genes_shared": int(keep.sum())}


def select_genes(adata, args):
    n_hvg = min(args.n_hvg, adata.n_vars)
    sc.pp.highly_variable_genes(
        adata, flavor="seurat_v3", n_top_genes=n_hvg,
        batch_key=args.batch_key, layer="counts", subset=False,
    )
    hvg = set(adata.var_names[adata.var["highly_variable"]])
    forced = [g.strip() for g in args.include_genes.split(",") if g.strip()]
    present = [g for g in forced if g in adata.var_names]
    dropped = [g for g in forced if g not in adata.var_names]
    if dropped:
        log(f"WARNING: {len(dropped)} include-genes absent from shared gene space: {dropped}")
    keep = sorted(hvg.union(present))
    log(f"HVG {len(hvg)} + forced markers {len(present)} -> {len(keep)} genes for the model")
    sub = adata[:, keep].copy()
    return sub, {"n_hvg_selected": len(hvg), "markers_forced": present, "markers_dropped": dropped,
                 "n_genes_model": len(keep)}


def run_scvi(adata, args):
    import scvi
    scvi.settings.seed = args.seed
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=args.batch_key)
    model = scvi.model.SCVI(
        adata, n_latent=args.n_latent, n_layers=args.n_layers,
        gene_likelihood=args.gene_likelihood, dispersion=args.dispersion,
    )
    model.train(max_epochs=args.max_epochs, accelerator=args.accelerator)
    latent = model.get_latent_representation()
    loss = None
    try:
        loss = float(model.history["elbo_train"].iloc[-1].item())
    except Exception:
        pass
    return latent, {"elbo_train_final": loss, "n_epochs": int(model.history["elbo_train"].shape[0])
                    if "elbo_train" in model.history else None}


def run_sysvi(adata, args):
    import scvi
    scvi.settings.seed = args.seed
    if not hasattr(scvi.external, "SysVI"):
        raise SystemExit("SysVI not available in this scvi-tools build; upgrade scvi-tools or use --model scvi")
    scvi.external.SysVI.setup_anndata(adata, layer="counts", batch_keys=[args.batch_key])
    model = scvi.external.SysVI(adata, n_latent=args.n_latent)
    model.train(max_epochs=args.max_epochs or 400, accelerator=args.accelerator)
    return model.get_latent_representation(), {"model": "sysvi"}


def main():
    args = parse_args()
    adata = load_anndata(args)
    adata, shared_info = filter_shared_genes(adata, args)
    adata, gene_info = select_genes(adata, args)
    gene_info = {**shared_info, **gene_info}
    log(f"training {args.model} (n_latent={args.n_latent}, n_layers={args.n_layers}, "
        f"likelihood={args.gene_likelihood}, dispersion={args.dispersion})")
    latent, train_info = (run_sysvi if args.model == "sysvi" else run_scvi)(adata, args)

    cols = [f"Z{i+1}" for i in range(latent.shape[1])]
    out = pd.DataFrame(latent, index=adata.obs_names, columns=cols)
    out.index.name = "barcode"
    out.to_csv(args.out)
    log(f"wrote latent {out.shape} -> {args.out}")

    if args.summary:
        summary = {"model": args.model, "n_cells": int(adata.n_obs), "batch_key": args.batch_key,
                   "n_latent": args.n_latent, "n_layers": args.n_layers,
                   "gene_likelihood": args.gene_likelihood, "dispersion": args.dispersion,
                   "seed": args.seed, **gene_info, **train_info}
        json.dump(summary, open(args.summary, "w"), indent=2)
        log(f"wrote summary -> {args.summary}")


if __name__ == "__main__":
    main()
