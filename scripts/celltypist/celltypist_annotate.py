#!/usr/bin/env python3
"""CellTypist annotation for seurat-pipeline step 05 (annotation.celltypist).

Called by run_celltypist() in utils.R with raw counts exported as MatrixMarket
(genes x cells). Normalizes to 10,000 counts per cell + log1p (what CellTypist
expects), runs each model, and writes one CSV per model
(<out-dir>/celltypist_<model>.csv: barcode, label, conf[, majority]) plus a
summary with each model's gene coverage on this data.

Uses the bundled env in scripts/scvi/ (celltypist is one of its dependencies).
Models are downloaded on first use into CellTypist's model folder
(~/.celltypist/data/models, or $CELLTYPIST_FOLDER).
"""
import argparse
import json

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import anndata as ad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--barcodes", required=True)
    ap.add_argument("--models", required=True, help="comma-separated CellTypist model names")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--majority-voting", action="store_true")
    a = ap.parse_args()

    import celltypist
    from celltypist import models as ct_models

    X = scipy.io.mmread(a.counts).T.tocsr()          # cells x genes
    genes = [l.rstrip("\n") for l in open(a.genes)]
    barcodes = [l.rstrip("\n") for l in open(a.barcodes)]
    adata = ad.AnnData(X=X.astype(np.float32), obs=pd.DataFrame(index=barcodes), var=pd.DataFrame(index=genes))
    adata.var_names_make_unique()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    summary = {"celltypist_version": celltypist.__version__, "majority_voting": a.majority_voting,
               "n_cells": adata.n_obs, "n_genes": adata.n_vars, "models": {}}
    present = set(adata.var_names)
    for name in [m.strip() for m in a.models.split(",") if m.strip()]:
        fname = name if name.endswith(".pkl") else name + ".pkl"
        key = fname[:-4]
        ct_models.download_models(model=fname, force_update=False)
        model = ct_models.Model.load(model=fname)
        feats = list(model.features)
        n_present = sum(f in present for f in feats)
        pred = celltypist.annotate(adata, model=model, majority_voting=a.majority_voting)
        labels = pred.predicted_labels
        out = pd.DataFrame({"barcode": barcodes,
                            "label": labels["predicted_labels"].astype(str).values,
                            "conf": np.round(pred.probability_matrix.max(axis=1).values, 4)})
        if a.majority_voting and "majority_voting" in labels:
            out["majority"] = labels["majority_voting"].astype(str).values
        out.to_csv(f"{a.out_dir}/celltypist_{key}.csv", index=False)
        desc = getattr(model, "description", None)
        summary["models"][key] = {
            "description": desc.get("details") if isinstance(desc, dict) else None,
            "n_model_genes": len(feats),
            "n_model_genes_present": n_present,
            "gene_coverage_pct": round(100 * n_present / max(len(feats), 1), 1),
            "n_labels": len(model.cell_types),
            "n_labels_predicted": int(out["label"].nunique()),
        }
        print(f"{key}: {n_present}/{len(feats)} model genes present; "
              f"{summary['models'][key]['n_labels_predicted']} labels predicted", flush=True)

    with open(a.summary, "w") as fh:
        json.dump(summary, fh, indent=2)


if __name__ == "__main__":
    main()
