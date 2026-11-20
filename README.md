# seurat-pipeline

A config-driven Seurat v5 pipeline for single-cell RNA-seq. Point one YAML file
at 10x Cell Ranger output, run the steps one at a time, and get a clustered,
annotated dataset plus a JSON summary per run that downstream tools (or LLMs)
can read directly.

Repo: https://github.com/sanghoonio/seurat-pipeline

> **Using [Claude Code](https://claude.com/claude-code)?** This repo ships with a
> `CLAUDE.md`. Say: *"Set this pipeline up for the Cell Ranger output at
> `/path/to/run` and iterate the QC for me."* See
> [Letting Claude Code drive](#letting-claude-code-drive).

## Pipeline

Eight steps, each a standalone R script. A step reads the output of one run of
the previous step and writes its own folder.

| Step | Script | What it does | Key outputs |
|------|--------|--------------|-------------|
| 00 | `00_cell_calling.R` | Optional. emptyDrops on the raw matrices; records which barcodes emptyDrops and/or Cell Ranger called | `cell_calling/<id>_barcodes.csv` |
| 01 | `01_load_qc_metrics.R` | Load the matrices (Cell Ranger cells, or a step 00 run's emptyDrops / union set), merge, per-cell QC metrics and their distributions | `seurat_raw.rds`, `qc_metrics.csv` |
| 02 | `02_ambient.R` | Ambient RNA correction per sample, SoupX or DecontX; per-gene removal table | `seurat_decontaminated.rds`, `ambient_removed_by_gene.csv` |
| 03 | `03_doublets.R` | scDblFinder per sample on every cell, before any threshold; calls and scores only | `doublet_calls.csv` |
| 04 | `04_filter.R` | QC thresholds on the corrected counts plus removal of called doublets; removed cells kept with their reason | `seurat_filtered.rds`, `seurat_filtered_out.rds` |
| 05 | `05_normalize.R` | SCTransform or LogNormalize, variable features, PCA | `seurat_normalized.rds`, `pca_variance.csv` |
| 06 | `06_cluster.R` | Integration (Harmony, RPCA, CCA, or scVI), neighbors, UMAP, clustering at several resolutions, `FindAllMarkers` | `seurat_clustered.rds`, `cluster_markers.csv`, `umap_coordinates.csv` |
| 07 | `07_annotate.R` | Azimuth, SingleR, CellTypist, a custom Seurat reference and marker-panel scores, each optional and independent; per-cluster table with a suggested label | `seurat_final.rds`, `cluster_annotations.csv`, `annotation/*.csv` |

Expensive work is split so variants of a later step reuse it: doublets are called
once per ambient run, normalization runs once per filter run, and markers once
per clustering.

## Runs, labels and inputs

Every step can run as several variants side by side, for example a SoupX and a
DecontX ambient run. Each variant is a **run block** in the config:

```yaml
02_ambient/decontx:
  ambient: {method: "decontx", round_to_int: true, decontx: {use_background: true}}

02_ambient/soupx:
  ambient: {method: "soupx", round_to_int: true, soupx: {contamination_range: [0.01, 0.8]}}

03_doublets/dx:
  input: "decontx"                 # reads outputs/02_ambient/decontx/
  doublets: {enabled: true, expected_rate: 0.05}
```

- **Keys:** `<step>` is the unlabeled run (label `default`); `<step>/<label>` is a
  labeled run. Labels use letters, digits, `_` and `-`.
- **Self-contained:** a block holds every setting its run uses. Keys left out take
  the built-in defaults (`config.template.yaml`); nothing is inherited from other
  blocks. A block may contain only the keys its step uses.
- **`input`:** the label of the previous step's run to read; absent means that
  step's unlabeled run. Settings a step needs from earlier steps (the cell set,
  the MT pattern, the normalization) come from the runs it reads.
- **Output folders:** each run writes to `outputs/<step>/<label>/` (RDS, CSVs,
  summary JSON, `plots/`). Rerunning a run clears and replaces only its own
  folder.
- **Run IDs:** every finished run gets a `run_id`, stamped into its RDS and
  summary. A step refuses to start unless its input run finished (`success` or
  `warning`), and refuses an RDS whose stamp does not match its run's summary.
- **Stale runs:** a run is stale when any run upstream of it was rerun, or its own
  block changed in the config since it ran. The report flags it; the next step
  warns when its input is stale.

## Requirements

- R >= 4.2 with Seurat >= 5
- Required R packages: `Seurat`, `yaml`, `jsonlite`, `ggplot2`, `Matrix`
- Packages needed only when a feature is configured. If a config asks for the
  feature and the package is missing, the step stops with an error rather than
  skipping it:
  - `SoupX` (`ambient.method: soupx`), `celda` (`decontx`)
  - `DropletUtils` (step 00)
  - `scDblFinder`, `SingleCellExperiment` (`doublets.enabled`)
  - `glmGamPoi` (SCTransform), `harmony` (`integration.method: harmony`)
  - `Azimuth`, `SeuratData` (`annotation.azimuth.reference`)
  - `SingleR`, `celldex` (`annotation.singler.reference`)
- Bundled Python env, set up once with `cd scripts/scvi && uv sync`, for scVI
  integration (step 06), CellTypist annotation (step 07) and the reference
  builder in `scripts/refs/`. Nothing Python loads unless one of these is used.
  See `scripts/scvi/README.md`.

## Quickstart

```bash
git clone https://github.com/sanghoonio/seurat-pipeline.git
cd seurat-pipeline
cp config.template.yaml config.yaml
$EDITOR config.yaml                      # samples, then one block per step

Rscript scripts/01_load_qc_metrics.R --config config.yaml
Rscript scripts/02_ambient.R         --config config.yaml
# ... read each run's summary and plots, adjust the config, then the next step
Rscript scripts/04_filter.R --config config.yaml --label dx_floor100 --no-plots

./run.sh config.yaml --report-only       # every run and its status
```

`--label` is needed when a step has several runs; with one run it can be
omitted. `--no-plots` skips PNGs. Unknown arguments are errors.

## Pointing the pipeline at Cell Ranger output

Each sample needs a unique `id` and its filtered matrix. Ambient correction and
cell calling also need the raw matrix. Either an `.h5` file or a 10x mtx
directory works.

**`cellranger count`:**

```yaml
samples:
  - id: "sample1"
    filtered: "data/sample1/outs/filtered_feature_bc_matrix.h5"
    raw: "data/sample1/outs/raw_feature_bc_matrix.h5"
```

**`cellranger multi` / Flex per-sample outputs:**

```yaml
samples:
  - id: "sample1"
    filtered: "data/run1/per_sample_outs/sample1/count/sample_filtered_feature_bc_matrix.h5"
    raw:      "data/run1/per_sample_outs/sample1/count/sample_raw_feature_bc_matrix.h5"
```

For Flex, `raw` is the gene-level `sample_raw_feature_bc_matrix.h5`, not the probe
matrix. Raw matrices contribute barcodes only; their features are restricted to
the filtered matrix's.

### Sample metadata

Any extra field under a sample becomes a per-cell metadata column:

```yaml
samples:
  - id: "pbmc_1k_v3"
    filtered: "data/pbmc_1k_v3/outs/filtered_feature_bc_matrix.h5"
    chemistry: "v3"
    donor: "donor1"
```

With many per-sample fields, keep them in a CSV and point the config at it:

```yaml
sample_metadata:
  file: "sample_metadata.csv"   # relative to base_dir
  id_column: "sample_id"        # matches samples[].id
  columns: ["condition"]        # null = all columns
```

Every sample needs exactly one row; a column that is also an inline field is an
error. Summaries record the file's path, md5 and column names, never its values.

## Outputs

Each run's folder, `outputs/<step>/<label>/`, holds:

| Step | Files |
|------|-------|
| 00 | `cell_calling/<id>_barcodes.csv`, `cell_calling/<id>_barcode_ranks.csv` |
| 01 | `seurat_raw.rds`, `qc_metrics.csv` |
| 02 | `seurat_decontaminated.rds`, `qc_metrics_decontaminated.csv` (with pre-correction totals), `ambient_removed_by_gene.csv` (per sample and gene: UMIs before, removed, share of the ambient profile) |
| 03 | `doublet_calls.csv` |
| 04 | `seurat_filtered.rds`; `seurat_filtered_out.rds` and `filtered_out_cells.csv` (removed cells with `removal_reason` and `failed_<threshold>`) |
| 05 | `seurat_normalized.rds`, `pca_variance.csv` |
| 06 | `seurat_clustered.rds`, `cluster_markers.csv`, `umap_coordinates.csv` |
| 07 | `seurat_final.rds` (`obj$cell_type` set), `cluster_annotations.csv`, `annotation/*.csv` (per-cell labels per SingleR reference, CellTypist model, custom reference), `filtered_out_annotation.csv` |

Every run also writes its summary JSON and `plots/*.png`.
`./run.sh config.yaml --report-only` writes `outputs/pipeline_report.json`: an
index of every run (defined in the config or found on disk) with its status,
plus every run's summary.

Every summary shares this envelope:

```json
{
  "step": "04_filter",
  "timestamp": "2026-11-20T14:23:01-0500",
  "project": "my-scrnaseq-project",
  "pipeline_commit": "1a2b3c4",
  "config": { "path": "/abs/path/config.yaml", "md5": "..." },
  "versions": { "R": "4.5.0", "Seurat": "...", "SeuratObject": "..." },
  "run": { "step": "04_filter", "label": "dx_floor100", "run_id": "...", "settings_hash": "...",
           "settings": { "...": "the run's block with defaults filled in" } },
  "input": { "step": "03_doublets", "label": "dx", "run_id": "...", "dir": "03_doublets/dx" },
  "lineage": [ { "step": "01_load_qc_metrics", "label": "default", "run_id": "...", "settings": {} } ],
  "inputs": { "...": "files this run read" },
  "params": { "...": "settings this run used" },

  "...": "step-specific fields",

  "outputs": { "rds": "04_filter/dx_floor100/seurat_filtered.rds" },
  "warnings": [],
  "status": "success"
}
```

`status` is `success`, `warning` (finished with entries in `warnings`) or
`failed` (with `error` and `details`). A run killed outright (out of memory, time limit) leaves no summary
and shows as `missing` in the report. Output paths are relative to
`paths.output_dir`.

In the report, each run's status is one of `success`, `warning`, `failed`,
`missing` (defined but never finished), `stale` (with the reasons), or
`not_in_config` (a folder on disk whose block was removed).

## Letting Claude Code drive

The summaries exist so an LLM can judge each step and propose config edits
without loading a Seurat object. What each carries:

- **00 `cell_calling_summary.json`**: emptyDrops params, per-sample barcode set
  sizes, stats of the recovered-only cells, knee and inflection.
- **01 `qc_summary.json`**: samples loaded; per sample and merged, cell and
  feature counts and quantiles of nCount, nFeature, percent.mt, percent.ribo and
  complexity.
- **02 `ambient_summary.json`**: method and params; per sample, contamination
  (and its basis), UMIs removed, top removed genes, quantiles after correction.
- **03 `doublets_summary.json`**: doublets called per sample and overall.
- **04 `filter_summary.json`**: per sample and total, cells removed by each
  threshold, removed only because ambient correction pushed them under a bound,
  doublets removed, cells kept.
- **05 `normalize_summary.json`**: method, split, variable features, PCA variance.
- **06 `cluster_summary.json`**: integration status, embedding used, PCA variance
  at the dims used, clusters per resolution, per cluster the dominant sample's
  share against its expected share, number of markers.
- **07 `annotation_summary.json`**: per method, status, gene coverage and labels
  found; per cluster, the suggested label and its source, Azimuth label share,
  score and low-confidence %, marker score and margin, top markers; optionally
  the labels of removed cells, diagnostics per level of `compare_by` columns,
  and how much of each marker set ambient correction removed.

`CLAUDE.md` holds the review checklist and operating rules. In the default
assisted mode Claude reads the summaries and proposes edits, and you approve
each run. In autonomous mode, which you ask for explicitly, it runs the steps,
edits the config, reruns only what a change affects, and logs every iteration
in `logs/iteration_log.md`.

## Configuration

`config.template.yaml` is the reference for every key, with comments. Top level
holds only shared settings; everything else is in run blocks:

| Key | Purpose |
|-----|---------|
| `project` | Name, description, random seed |
| `paths` | Base, output and log directories |
| `samples` | Samples with `filtered` / `raw` matrix paths and any metadata fields |
| `sample_metadata` | Optional per-sample metadata CSV |
| `00_cell_calling` | `cell_calling.emptydrops` |
| `01_load_qc_metrics` | `cell_calling.cells` (cellranger / emptydrops / union), `qc.metrics` |
| `02_ambient` | `ambient` (SoupX or DecontX) |
| `03_doublets` | `doublets` |
| `04_filter` | `qc.thresholds`, `qc.prefilter_umap`, `filter.remove_doublets` |
| `05_normalize` | `normalization` (incl. `split_by`), `reduction.npcs` |
| `06_cluster` | `integration`, `reduction.dims_use`, `reduction.umap`, `clustering`, `markers` |
| `07_annotate` | `annotation` (`azimuth`, `singler`, `celltypist`, `custom_reference`, `markers`, `removed_cells`, `compare_by`) |

The config is validated as a whole every time a step starts. Unknown keys,
unknown steps, keys a step does not use, old-layout top-level sections, bad
labels and an `input` that names no defined run are all errors.

## Layout

```
seurat-pipeline/
├── README.md
├── CLAUDE.md
├── config.template.yaml
├── run.sh                 # ./run.sh <config> --report-only
├── scripts/
│   ├── utils.R            # config, run blocks, lineage, summaries, shared helpers
│   ├── 00_cell_calling.R … 07_annotate.R
│   ├── report.R           # run index + staleness (used by run.sh)
│   ├── scvi/              # bundled Python env + scVI/sysVI integration (see README)
│   ├── celltypist/        # CellTypist annotation, run by step 07
│   └── refs/              # build a custom Seurat reference from an .h5ad (see README)
├── outputs/               # generated, gitignored
└── logs/                  # generated, gitignored
```
