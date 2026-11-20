#!/usr/bin/env python3
"""Export a reference atlas .h5ad (e.g. a CELLxGENE download) for build_seurat_reference.R.

Takes a label-stratified subsample (at most --max-per-label cells per label, then at
most --max-cells overall), uses raw integer counts (adata.raw.X if present, else X),
switches gene IDs to symbols (--symbol-column, default feature_name, as on CELLxGENE)
summing duplicated symbols, and writes:
  counts.mtx (genes x cells), genes.txt, barcodes.txt, obs.csv, export_summary.json

Run with the bundled env: scripts/scvi/.venv/bin/python scripts/refs/export_h5ad.py ...
"""
import argparse
import json
import os

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--label-columns", required=True, help="comma-separated obs columns to keep as labels")
    ap.add_argument("--stratify-by", default=None, help="obs column to stratify on (default: first label column)")
    ap.add_argument("--keep-obs", default="", help="extra comma-separated obs columns to export (e.g. donor_id,assay)")
    ap.add_argument("--symbol-column", default="feature_name")
    ap.add_argument("--max-per-label", type=int, default=2000)
    ap.add_argument("--max-cells", type=int, default=100000)
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()

    labels = [c.strip() for c in a.label_columns.split(",") if c.strip()]
    keep = labels + [c.strip() for c in a.keep_obs.split(",") if c.strip() and c.strip() not in labels]
    strat = a.stratify_by or labels[0]
    os.makedirs(a.out_dir, exist_ok=True)

    adata = ad.read_h5ad(a.h5ad)
    missing = [c for c in set(keep + [strat]) if c not in adata.obs.columns]
    if missing:
        raise SystemExit(f"obs columns not found: {missing}; available: {list(adata.obs.columns)}")

    rng = np.random.default_rng(a.seed)
    groups = adata.obs.groupby(strat, observed=True).indices
    idx = np.concatenate([rng.choice(ix, size=min(len(ix), a.max_per_label), replace=False)
                          for ix in groups.values()])
    if len(idx) > a.max_cells:
        idx = rng.choice(idx, size=a.max_cells, replace=False)
    idx = np.sort(idx)

    use_raw = adata.raw is not None
    X = (adata.raw.X if use_raw else adata.X)[idx]
    var = adata.raw.var if use_raw else adata.var
    X = sp.csr_matrix(X)
    sample = X.data[: min(len(X.data), 100000)]
    if len(sample) and not np.allclose(sample, np.round(sample)):
        raise SystemExit("counts are not integers (normalized data?); a Seurat SCT reference needs raw counts")

    symbols = var[a.symbol_column].astype(str).values if a.symbol_column in var.columns else np.asarray(var.index, dtype=str)
    uniq, inv = np.unique(symbols, return_inverse=True)
    agg = sp.csr_matrix((np.ones(len(inv)), (np.arange(len(inv)), inv)), shape=(len(inv), len(uniq)))
    Xg = (X @ agg).tocsc()                                   # cells x unique symbols, duplicates summed

    obs = adata.obs.iloc[idx][keep].copy()
    barcodes = [f"ref{i}" for i in range(len(idx))]          # source obs names can repeat across datasets
    obs.insert(0, "barcode", barcodes)
    scipy.io.mmwrite(os.path.join(a.out_dir, "counts.mtx"), Xg.T.tocoo().astype(np.int32))
    with open(os.path.join(a.out_dir, "genes.txt"), "w") as fh:
        fh.write("\n".join(uniq) + "\n")
    with open(os.path.join(a.out_dir, "barcodes.txt"), "w") as fh:
        fh.write("\n".join(barcodes) + "\n")
    obs.to_csv(os.path.join(a.out_dir, "obs.csv"), index=False)
    summary = {
        "source": os.path.abspath(a.h5ad), "source_bytes": os.path.getsize(a.h5ad),
        "n_cells_source": int(adata.n_obs), "n_cells_exported": int(len(idx)), "n_genes": int(len(uniq)),
        "counts_from": "raw.X" if use_raw else "X", "stratify_by": strat,
        "max_per_label": a.max_per_label, "max_cells": a.max_cells, "seed": a.seed,
        "label_counts": {str(k): int(v) for k, v in obs[strat].value_counts().items()},
    }
    with open(os.path.join(a.out_dir, "export_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"exported {len(idx)} of {adata.n_obs} cells, {len(uniq)} genes -> {a.out_dir}")


if __name__ == "__main__":
    main()
