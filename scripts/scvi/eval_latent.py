#!/usr/bin/env python
"""Score an integration latent on two general criteria: (1) cross-batch mixing and
(2) survival of a marker-defined (often rare/confounded) population. Reads a latent CSV
plus the exported counts/metadata, clusters on the latent, and prints a compact JSON
summary line. Batch key and marker panel are arguments, so it works for any project.

Example:
  eval_latent.py --latent out/scvi/nlayers3/scvi_latent.csv --io out/scvi/nlayers3 \
                 --batch-key chemistry --markers FCGR3B,CSF3R,CXCR2,S100A8,S100A9
"""
import argparse, json, os
import numpy as np, pandas as pd, scanpy as sc, scipy.io, scipy.sparse as sp

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--latent", required=True)
    ap.add_argument("--io", default=".", help="dir with counts.mtx/genes.txt/barcodes.txt/metadata.csv")
    ap.add_argument("--batch-key", default="batch", help="metadata column integration corrected over")
    ap.add_argument("--markers", default="", help="comma-separated marker genes for the population to track")
    ap.add_argument("--label", default=None)
    ap.add_argument("--res", type=float, default=0.5)
    ap.add_argument("--summary-json", default=None, help="run summary JSON for the config line")
    a = ap.parse_args()
    label = a.label or os.path.basename(os.path.dirname(a.latent)) or "run"
    bk = a.batch_key

    lat = pd.read_csv(a.latent, index_col=0)
    meta = pd.read_csv(f"{a.io}/metadata.csv", index_col=0).reindex(lat.index)
    if bk not in meta.columns:
        raise SystemExit(f"batch key '{bk}' not in metadata columns {list(meta.columns)}")
    ad = sc.AnnData(X=lat.values, obs=meta.copy())
    ad.obs_names = lat.index; ad.obsm["X_lat"] = lat.values
    ad.obs[bk] = ad.obs[bk].astype("category")

    sc.pp.neighbors(ad, use_rep="X_lat", n_neighbors=30)
    sc.tl.leiden(ad, resolution=a.res, flavor="igraph", n_iterations=2, directed=False, random_state=0)

    # (1) mixing: per-cluster batch purity = dominant-category share (general to N batches)
    comp = pd.crosstab(ad.obs["leiden"], ad.obs[bk])
    total = comp.sum(1)
    purity = comp.max(1) / total
    mixed = int(total[purity < 0.8].sum()); tot = int(total.sum())
    pct_mixed = round(100 * mixed / tot, 1)

    # (2) local kNN mixing: mean fraction of a cell's neighbours from a different batch
    conn = ad.obsp["distances"].tocsr(); code = ad.obs[bk].cat.codes.values
    of = [np.mean(code[conn.indices[conn.indptr[i]:conn.indptr[i+1]]] != code[i])
          for i in range(ad.n_obs) if conn.indptr[i+1] > conn.indptr[i]]
    knn_other = round(float(np.mean(of)), 3)

    out = {"label": label, "n_clusters": int(ad.obs["leiden"].nunique()),
           "pct_cells_well_mixed": pct_mixed, "knn_other_batch_frac": knn_other}

    # (3) marker population survival (optional)
    markers = [g.strip() for g in a.markers.split(",") if g.strip()]
    if markers:
        mat = sp.csr_matrix(scipy.io.mmread(f"{a.io}/counts.mtx")).T.tocsr()
        genes = [l.strip() for l in open(f"{a.io}/genes.txt")]
        barcodes = [l.strip() for l in open(f"{a.io}/barcodes.txt")]
        cd = sc.AnnData(X=mat); cd.var_names = genes; cd.obs_names = barcodes
        cd = cd[lat.index].copy(); sc.pp.normalize_total(cd, target_sum=1e4); sc.pp.log1p(cd)
        present = [g for g in markers if g in cd.var_names]
        sc.tl.score_genes(cd, present, score_name="pop")
        ad.obs["pop"] = cd.obs["pop"].values
        by = ad.obs.groupby("leiden", observed=True)["pop"].mean().sort_values(ascending=False)
        top = by.index[0]; sub = ad.obs[ad.obs["leiden"] == top]
        dom = sub[bk].value_counts()
        out.update({"markers_present": present, "pop_cluster_n": int(len(sub)),
                    "pop_cluster_dominant_batch": f"{dom.index[0]} {round(100*dom.iloc[0]/len(sub))}%",
                    "pop_score_margin": round(float(by.iloc[0] - by.iloc[1]), 2)})

    if a.summary_json and os.path.exists(a.summary_json):
        s = json.load(open(a.summary_json))
        out["config"] = f"n_latent={s.get('n_latent')} n_layers={s.get('n_layers')} disp={s.get('dispersion')}"
    print(json.dumps(out))

if __name__ == "__main__":
    main()
