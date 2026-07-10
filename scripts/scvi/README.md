# scVI / sysVI integration (opt-in)

Count-based deep integration for step 04, for **large cross-technology batch effects**
(e.g. 10x Flex probe-panel vs 3′ whole-transcriptome) that Harmony/rpca can't bridge.

This is an **optional Python subsystem**. The core pipeline is pure R and loads none of
this unless a config sets `integration.method: scvi`. The R side
(`utils.R::run_scvi_integration`) shells out to `scvi_integrate.py` here via `system2()`;
it never links Python into the R process.

## One-time setup

```bash
cd scripts/scvi
uv sync          # creates ./.venv with scvi-tools + torch (large; CPU/MPS is fine)
```

The `.venv/` is gitignored and provisioned per clone. Nothing else is needed — the R
helper defaults `integration.scvi.python`/`.script` to this bundled env, so a project's
config only needs `method: scvi` plus any tunables.

## Files

| file | purpose |
|---|---|
| `scvi_integrate.py` | mtx + sidecars in → shared-gene filter → HVG(∪markers) → SCVI/SysVI train → latent CSV out. `--model scvi\|sysvi`. |
| `pyproject.toml` | uv env spec (scvi-tools, scanpy, torch, …). |
| `eval_latent.py` | score a latent on cross-batch mixing + a marker-defined rare population's survival (cheap config comparison before a full step 04 run). |

UMAP *visualization* is left to the R pipeline (`04`/`05` emit `*_umap_*.png` from the
clustered object via `DimPlot`); this subsystem only trains latents and scores them.

## Config knobs

See `config.template.yaml` → `integration.scvi`. Key ones: `n_latent`, `n_layers`
(capacity), `gene_likelihood`, `dispersion`, `n_hvg`, `min_cells_per_batch` (shared gene
space), `include_genes` (force markers into the model), `run_label` (tags coexisting
runs), `max_epochs`.

## How the R↔Python contract works

`run_scvi_integration` writes raw counts (`counts.mtx`) + `genes.txt` / `barcodes.txt` /
`metadata.csv` to `<output_dir>/scvi/<run_label>/`, runs this script, and reads back
`scvi_latent.csv` aligned to cells by barcode. It caches on the call parameters, so a
re-run with unchanged params skips the (slow) training.
