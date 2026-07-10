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
| 2 | `02_soupx.R` | Estimate and remove ambient RNA contamination per sample using SoupX | `seurat_decontaminated.rds`, `soupx_summary.json` |
| 3 | `03_filter.R` | Apply QC thresholds, optional doublet detection (scDblFinder) | `seurat_filtered.rds`, `filter_summary.json` |
| 4 | `04_normalize_cluster.R` | SCTransform / LogNormalize, PCA, optional batch integration (Harmony / RPCA / CCA, or opt-in scVI for large cross-technology gaps — see `scripts/scvi/README.md`), UMAP, multi-resolution clustering | `seurat_clustered.rds`, `umap_coordinates.csv`, `pca_variance.csv`, `cluster_summary.json` |
| 5 | `05_annotate.R` | Optional Azimuth, optional SingleR, marker-gene scoring, `FindAllMarkers`, build per-cluster annotation table | `seurat_final.rds`, `cluster_annotations.csv`, `cluster_markers.csv`, `annotation_summary.json` |

After all steps complete, `run.sh` collects every `*_summary.json` plus
the config used into a single `outputs/pipeline_report.json` — the merged
state of the entire run, in one file an LLM or report generator can read.

## Requirements

- R >= 4.2 (Seurat v5 recommended)
- Required R packages: `Seurat`, `SoupX`, `yaml`, `jsonlite`, `ggplot2`
- Feature-gated packages — **required only when the feature is enabled**. If a config
  requests the feature but the package is missing, the step **errors out** (it does not
  silently skip):
  - `scDblFinder`, `SingleCellExperiment` — when `qc.doublet_detection: true`
  - `harmony` — when `integration.method: harmony`
  - `Azimuth` — when `annotation.azimuth_reference` is set
  - `SingleR`, `celldex` — when `annotation.singler_reference` is set
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
Ranger H5 matrix. SoupX additionally needs the raw H5 (skipped per-sample
if missing).

**Standard `cellranger count` layout:**

```
data/sample1/outs/
├── filtered_feature_bc_matrix.h5
└── raw_feature_bc_matrix.h5
```

```yaml
samples:
  - id: "sample1"
    filtered_h5: "data/sample1/outs/filtered_feature_bc_matrix.h5"
    raw_h5: "data/sample1/outs/raw_feature_bc_matrix.h5"
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
    filtered_h5: "data/run1/per_sample_outs/sample1/count/sample_filtered_feature_bc_matrix.h5"
    raw_h5:      "data/run1/per_sample_outs/sample1/count/sample_raw_feature_bc_matrix.h5"
```

### Sample metadata

Any extra fields you add under a sample become per-cell metadata columns
on the merged Seurat object — no code changes needed. For example:

```yaml
samples:
  - id: "pbmc_1k_v3"
    filtered_h5: "data/pbmc_1k_v3/outs/filtered_feature_bc_matrix.h5"
    raw_h5:      "data/pbmc_1k_v3/outs/raw_feature_bc_matrix.h5"
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
| `seurat_decontaminated.rds` | After SoupX correction |
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

  "...": "step-specific fields here",

  "outputs": { "rds": "...", "csv": "..." },
  "status": "success"
}
```

Errors are emitted to stderr in the same shape:

```json
{
  "error": "Filtered H5 file not found",
  "details": { "sample": "sample1", "path": "..." },
  "timestamp": "..."
}
```

## Letting Claude Code drive

The JSON summaries exist so an LLM can read them, judge whether the
pipeline made good decisions, and propose config edits *without ever
loading the Seurat object into memory*. Each one carries the numbers a
reviewer would actually want:

- **`qc_summary.json`** — per-sample cells loaded, genes detected, median
  `nFeature`/`nCount`/`percent.mt`, and `q25/q50/q75/q95/q99` quantiles of
  `percent.mt`. Spot a too-strict (or too-loose) MT cutoff in one read.
- **`soupx_summary.json`** — per-sample estimated contamination, cells
  affected, samples skipped.
- **`filter_summary.json`** — cells before/after each criterion, per-sample
  retention rates, doublet counts. Catches "we silently dropped 70% of
  sample X" before it pollutes downstream results.
- **`cluster_summary.json`** — PCA variance curve, suggested PC count,
  per-resolution cluster counts, per-sample cluster composition. Tells
  you whether `dims_use` looks right and whether one sample is driving
  an entire cluster.
- **`annotation_summary.json`** — per-method status (Azimuth, SingleR),
  median Azimuth prediction score, low-confidence cell %, top markers
  per cluster, suggested label table.

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
- *"Sample MAD-03 lost 60% of cells in filtering — figure out which
  threshold caused it and propose a fix."*
- *"Add a new sample to the config and re-run only the steps that need it."*

Autonomous mode is opt-in: Claude only skips per-step confirmation if
you say so explicitly. Otherwise it falls back to the assisted loop.

## Configuration

`config.template.yaml` is the canonical reference and is heavily
commented inline, including tuning notes for mouse data, 10x Flex
(probe-based) data, and tissue-specific Azimuth references. Top-level
sections:

| Section | Purpose |
|---------|---------|
| `project` | Name, description, random seed |
| `paths` | Where to write outputs and logs |
| `samples` | List of samples with H5 paths and arbitrary metadata |
| `qc` | QC thresholds, MT/ribo patterns, doublet detection toggle |
| `soupx` | Ambient RNA correction settings |
| `normalization` | SCTransform vs LogNormalize, variables to regress |
| `reduction` | PCA + UMAP parameters |
| `clustering` | Algorithm and resolutions |
| `integration` | Batch correction method (Harmony / RPCA / CCA / scVI) |
| `annotation` | Azimuth ref, SingleR ref, custom marker gene dict |

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
│   ├── 02_soupx.R
│   ├── 03_filter.R
│   ├── 04_normalize_cluster.R
│   ├── 05_annotate.R
│   └── scvi/        # opt-in scVI/sysVI Python subsystem (see scripts/scvi/README.md)
├── outputs/         # generated, gitignored
└── logs/            # generated, gitignored
```
