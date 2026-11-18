# seurat-pipeline

A reproducible, config-driven Seurat v5 pipeline for single-cell RNA-seq.
Drop in any 10x Cell Ranger output, edit one YAML file, run one script —
get a clustered, annotated dataset plus structured JSON summaries that
downstream tools (or LLMs) can read directly.

Repo: https://github.com/sanghoonio/seurat-pipeline

> **Using [Claude Code](https://claude.com/claude-code)?** This repo ships
> with a `CLAUDE.md`. Skip the README and just say:
> *"Set this pipeline up for the Cell Ranger output at `/path/to/run` and
> iterate the QC for me in autonomous mode."* See
> [Letting Claude Code drive](#letting-claude-code-drive) below.

## Pipeline

Five sequential, modular steps. Each step is a standalone R script that
reads a single YAML config, writes its own RDS checkpoint and a JSON
summary, and can be re-run independently.

| Step | Script | What it does | Key outputs |
|------|--------|--------------|-------------|
| 1 | `01_load_qc_metrics.R` | Load Cell Ranger H5 files, merge samples, compute per-cell QC metrics (`nFeature`, `nCount`, `percent.mt`, `percent.ribo`, complexity) | `seurat_raw.rds`, `qc_metrics.csv`, `qc_summary.json` |
| 2 | `02_ambient.R` | Estimate and remove ambient RNA contamination per sample (SoupX or DecontX, `ambient.method`) | `seurat_decontaminated.rds`, `ambient_summary.json` |
| 3 | `03_filter.R` | Apply QC thresholds, optional doublet detection (scDblFinder) | `seurat_filtered.rds`, `filter_summary.json` |
| 4 | `04_normalize_cluster.R` | SCTransform / LogNormalize, PCA, optional batch integration (Harmony / RPCA / CCA, or opt-in scVI for large cross-technology gaps — see `scripts/scvi/README.md`), UMAP, multi-resolution clustering | `seurat_clustered.rds`, `umap_coordinates.csv`, `pca_variance.csv`, `cluster_summary.json` |
| 5 | `05_annotate.R` | Optional Azimuth, optional SingleR, marker-gene scoring, `FindAllMarkers`, build per-cluster annotation table | `seurat_final.rds`, `cluster_annotations.csv`, `cluster_markers.csv`, `annotation_summary.json` |

After all steps complete, `run.sh` collects every `*_summary.json` plus
the config used into a single `outputs/pipeline_report.json` — the merged
state of the entire run, in one file an LLM or report generator can read.

## Requirements

- R >= 4.2 (Seurat v5 recommended)
- Required R packages: `Seurat`, `yaml`, `jsonlite`, `ggplot2`
- Feature-gated packages — **required only when the feature is enabled**. If a config
  requests the feature but the package is missing, the step **errors out** (it does not
  silently skip):
  - `SoupX` — when `ambient.method: soupx` (the default)
  - `celda` — when `ambient.method: decontx`
  - `scDblFinder`, `SingleCellExperiment` — when `qc.doublets.remove: true`
  - `harmony` — when `integration.method: harmony`
  - `Azimuth` — when `annotation.azimuth.reference` is set
  - `SingleR`, `celldex` — when `annotation.singler.reference` is set
- Opt-in scVI subsystem — for `integration.method: scvi` only. A separate Python env,
  set up once with `cd scripts/scvi && uv sync` (installs scvi-tools + torch). Nothing
  Python loads unless you select scVI. See `scripts/scvi/README.md`.
- Python 3 with `pyyaml` — used by `run.sh` to merge JSON summaries
  (the pipeline still completes without it; only `pipeline_report.json` is skipped)

## Quickstart

```bash
git clone https://github.com/sanghoonio/seurat-pipeline.git my-analysis
cd my-analysis

# 1. Copy and edit the config for your data
cp config.template.yaml config.yaml
$EDITOR config.yaml

# 2. Run the full pipeline
./run.sh config.yaml
```

Outputs land in `outputs/`, logs in `logs/`. To re-run a single step
after editing the config:

```bash
Rscript scripts/04_normalize_cluster.R --config config.yaml
Rscript scripts/04_normalize_cluster.R --config config.yaml --no-plots  # skip PNGs
```

## Pointing the pipeline at Cell Ranger output

The minimum information per sample is a unique id and the path to a Cell
Ranger H5 matrix. Ambient RNA correction (step 02) additionally needs the raw
H5 (skipped per-sample if missing).

**Standard `cellranger count` layout:**

```
data/sample1/outs/
├── filtered_feature_bc_matrix.h5
└── raw_feature_bc_matrix.h5
```

```yaml
samples:
  - id: "sample1"
    filtered: "data/sample1/outs/filtered_feature_bc_matrix.h5"
    raw: "data/sample1/outs/raw_feature_bc_matrix.h5"
```

**`cellranger multi` / Flex per-sample-outs layout:**

```
data/run1/per_sample_outs/sample1/count/
├── sample_filtered_feature_bc_matrix.h5
└── sample_raw_feature_bc_matrix.h5
```

```yaml
samples:
  - id: "sample1"
    filtered: "data/run1/per_sample_outs/sample1/count/sample_filtered_feature_bc_matrix.h5"
    raw:      "data/run1/per_sample_outs/sample1/count/sample_raw_feature_bc_matrix.h5"
```

### Sample metadata

Any extra fields you add under a sample become per-cell metadata columns
on the merged Seurat object — no code changes needed. For example:

```yaml
samples:
  - id: "pbmc_1k_v3"
    filtered: "data/pbmc_1k_v3/outs/filtered_feature_bc_matrix.h5"
    raw:      "data/pbmc_1k_v3/outs/raw_feature_bc_matrix.h5"
    chemistry: "v3"
    donor: "donor1"
    condition: "control"
```

becomes:

```r
obj <- readRDS("outputs/seurat_final.rds")
head(obj@meta.data[, c("sample_id", "chemistry", "donor", "condition")])
#>                                sample_id chemistry  donor condition
#> pbmc_1k_v3_AAACCCAAGGAGAGTA-1 pbmc_1k_v3        v3 donor1   control
#> ...

table(obj$chemistry, obj$cell_type)   # use them for downstream splits
DimPlot(obj, group.by = "condition")  # or as plot facets
```

The merged object's `sample_id` column is always set; everything else is
whatever you put in the YAML. Use this for `condition`, `batch`, `donor`,
`timepoint`, technical replicate IDs — anything you'd want available as a
Seurat metadata column or as a `group.by` / `split.by` argument later.

## Outputs

Everything lives under `outputs/`:

| File | What it contains |
|------|------------------|
| `seurat_raw.rds` | Merged, unfiltered object with raw counts and QC metrics |
| `seurat_decontaminated.rds` | After ambient RNA correction |
| `seurat_filtered.rds` | After QC threshold and doublet filtering |
| `seurat_clustered.rds` | After normalization, PCA, UMAP, clustering |
| `seurat_final.rds` | After cell type annotation (`obj$cell_type` set) |
| `qc_metrics.csv` / `qc_metrics_decontaminated.csv` | Per-cell QC metrics |
| `umap_coordinates.csv` | Per-cell UMAP, cluster, sample (for plotting outside R) |
| `pca_variance.csv` | Per-PC variance explained |
| `cluster_markers.csv` | `FindAllMarkers` output |
| `cluster_annotations.csv` | Per-cluster labels (Azimuth, SingleR, marker-based, suggested) |
| `*_summary.json` | One per step, structured summary of inputs/outputs/decisions |
| `pipeline_report.json` | Merged summary of every step + the config used |
| `plots/*.png` | All plots generated during the run |

Each `*_summary.json` shares this envelope:

```json
{
  "step": "01_load_qc_metrics",
  "timestamp": "2026-01-15T14:23:01-0500",
  "project": "my-scrnaseq-project",
  "seed": 42,
  "pipeline_commit": "1a2b3c4",
  "config": { "path": "/abs/path/config.yaml", "md5": "..." },
  "versions": { "R": "4.5.0", "Seurat": "...", "SeuratObject": "..." },
  "inputs": { "...": "files/matrices this step read" },
  "params": { "...": "config values this step actually used" },

  "...": "step-specific fields here",

  "outputs": { "rds": "seurat_raw.rds", "csv": "qc_metrics.csv" },
  "warnings": [],
  "status": "success"
}
```

`status` is `success`, `warning` (succeeded with entries in `warnings`),
`skipped` (step disabled) or `failed`. Each step deletes its previous summary
once its config is read; any later error writes a `failed` summary (with
`error` and `details`) in its place. A job killed outright (e.g. out of
memory) leaves no summary, shown as `missing` in the report's `_index`. Output
paths are relative to `paths.output_dir`. After running steps individually,
rebuild `pipeline_report.json` with `./run.sh config.yaml --report-only`.

## Letting Claude Code drive

The JSON summaries exist so an LLM can read them, judge whether the
pipeline made good decisions, and propose config edits *without ever
loading the Seurat object into memory*. Each one carries the numbers a
reviewer would actually want:

- **`cell_calling_summary.json`** (step 00, opt-in) — emptyDrops params,
  per-sample barcode set sizes, and stats of the recovered-only cells.
- **`qc_summary.json`** — samples configured vs loaded; per sample and merged:
  cell and feature counts, quantiles (min/q01/q05/…/q99/max) of nCount,
  nFeature, percent.mt, percent.ribo and complexity, and
  `threshold_preview_pct` (what each step 03 threshold would remove).
- **`ambient_summary.json`** — method (SoupX/DecontX) and its params;
  per sample: estimated contamination (and its basis), UMIs removed, top
  removed genes, post-correction quantiles; merged threshold preview on the
  corrected counts.
- **`filter_summary.json`** — per sample and total: cells before, removed by
  each threshold (and by any), doublets removed, cells after. Catches "we
  silently dropped 70% of sample X" before it pollutes downstream results.
- **`cluster_summary.json`** — normalization, integration status (applied /
  disabled / skipped, with reason), the embedding used, PCA variance at the
  dims used, per-resolution cluster counts, and per cluster the dominant
  sample's share against its expected share.
- **`annotation_summary.json`** — per-method status and params; per cluster:
  suggested label and its source, Azimuth label share / median score /
  low-confidence %, marker score and margin, top markers.

This repo ships with a `CLAUDE.md` documenting the QC iteration checklist
and operating rules. Pick one of two modes:

### Assisted mode

Claude reads the JSONs and proposes edits; you approve each re-run.
Useful when you want to learn what the pipeline is doing.

```
1. ./run.sh
2. "Read outputs/pipeline_report.json and tell me which QC thresholds
   need adjusting and why."
3. Approve the suggested config edits.
4. Re-run only the affected steps.
5. Repeat until clean.
```

### Autonomous mode

Hand the whole loop to Claude in one instruction:

> *"Run this pipeline on the Cell Ranger output at `/path/to/run` in
> autonomous mode. Iterate the config on your own using the JSON
> summaries until the QC looks clean, log every parameter you tried in
> `logs/iteration_log.md`, and just report to me at the end."*

Claude will wire up the config, run the pipeline, read the JSONs, edit
`config.yaml`, re-run only the affected steps, append a structured entry
to `logs/iteration_log.md` (iteration number, what was wrong, what
changed, resulting metrics), and repeat — usually 4–5 rounds — until the
checklist in `CLAUDE.md` is satisfied. It then writes a short summary
covering final cell counts, cluster count, integration status,
annotation confidence, and anything still suspicious.

The iteration log is a normal markdown file. Read it at the end to audit
every decision and roll back any change you don't agree with — every
config and JSON from the run is on disk.

Other useful prompts:

- *"My data is mouse PBMCs — update the config and marker genes accordingly."*
- *"sample_03 lost 60% of cells in filtering — figure out which
  threshold caused it and propose a fix."*
- *"Add a new sample to the config and re-run only the steps that need it."*

Autonomous mode is opt-in: Claude only skips per-step confirmation if
you say so explicitly. Otherwise it falls back to the assisted loop.

## Configuration

`config.template.yaml` is the canonical reference and is heavily
commented inline, including tuning notes for mouse data, 10x Flex
(probe-based) data, and tissue-specific Azimuth references. Any key left out
takes the template's value; keys from the previous layout are rejected with a
message saying where each moved. Method-specific options are nested under the
method's name (e.g. `ambient.decontx`, `normalization.sctransform`). Top-level
sections:

| Section | Purpose |
|---------|---------|
| `project` | Name, description, random seed |
| `paths` | Where to write outputs and logs |
| `samples` | List of samples with `filtered` / `raw` matrix paths and arbitrary metadata |
| `cell_calling` | Which cells to use: Cell Ranger filtered, emptyDrops, or their union |
| `qc` | `metrics` (MT/ribo patterns), `thresholds` (step 03 filters), `doublets` |
| `ambient` | Ambient RNA correction settings (SoupX or DecontX) |
| `normalization` | SCTransform vs LogNormalize, variables to regress, per-method options |
| `reduction` | PCs computed/used, `umap` parameters |
| `clustering` | Algorithm, resolutions, default resolution |
| `integration` | Batch correction method (Harmony / RPCA / CCA / scVI) |
| `annotation` | `azimuth` (reference, primary/plot levels), `singler`, `markers` |

## Layout

```
seurat-pipeline/
├── README.md
├── CLAUDE.md
├── config.template.yaml
├── run.sh
├── scripts/
│   ├── utils.R
│   ├── 00_cell_calling.R
│   ├── 01_load_qc_metrics.R
│   ├── 02_ambient.R
│   ├── 03_filter.R
│   ├── 04_normalize_cluster.R
│   ├── 05_annotate.R
│   └── scvi/        # opt-in scVI/sysVI Python subsystem (see scripts/scvi/README.md)
├── outputs/         # generated, gitignored
└── logs/            # generated, gitignored
```
