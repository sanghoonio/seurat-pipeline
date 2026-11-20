# Bundled Python env: scVI / sysVI, CellTypist, reference export

One uv-managed env serves every Python part of the pipeline:

- **scVI / sysVI integration** (step 06, `integration.method: scvi`): count-based deep
  integration for large cross-technology batch effects (e.g. 10x Flex probe panel vs 3'
  whole transcriptome) that Harmony / RPCA can't bridge.
- **CellTypist annotation** (step 07, `annotation.celltypist.models`), via
  `scripts/celltypist/celltypist_annotate.py`.
- **Reference export** for a custom Seurat reference, `scripts/refs/export_h5ad.py`.

The core pipeline is R and loads none of this unless a config uses one of them. The R
side shells out with `system2()`; Python never links into the R process.

## One-time setup

```bash
cd scripts/scvi
uv sync          # creates ./.venv with scvi-tools, torch, scanpy, celltypist (large; CPU is fine)
```

The `.venv/` is gitignored and provisioned per clone. Configs need no paths: the R
helpers default `python`/`script` to this env and the bundled scripts.

## Files

| file | purpose |
|---|---|
| `scvi_integrate.py` | mtx + sidecars in -> shared-gene filter -> HVG (+ forced markers) -> SCVI/SysVI training -> latent CSV out. `--model scvi\|sysvi`. |
| `pyproject.toml` | env spec (scvi-tools, scanpy, celltypist, ...). |
| `eval_latent.py` | score a latent on cross-batch mixing and survival of a marker-defined rare population (a cheap comparison before a full step 06 run). |

UMAPs are drawn by the R steps; this subsystem only trains latents and scores them.

## scVI config

`06_cluster` block, `integration.scvi`: `n_latent`, `n_layers` (capacity),
`gene_likelihood`, `dispersion`, `n_hvg`, `min_cells_per_batch` (shared gene space),
`include_genes` (force markers into the model), `max_epochs`, `seed`.

## How the R-Python contract works

`run_scvi_integration` (utils.R) writes raw counts (`counts.mtx`), `genes.txt`,
`barcodes.txt` and `metadata.csv`, runs `scvi_integrate.py`, and reads back
`scvi_latent.csv`, aligned to cells by barcode.

The cache lives in the folder of the 05_normalize run the 06 run reads:
`outputs/05_normalize/<label>/scvi/<key>/`. The key hashes the scVI parameters, the
input cell and gene names, and the 05 run's run_id. So 06 runs that differ only in UMAP
or clustering reuse one trained model, and rerunning the 05 run (which clears its
folder) or anything upstream of it forces retraining.
