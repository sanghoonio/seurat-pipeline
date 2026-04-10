# seurat-pipeline

A reproducible, config-driven Seurat v5 pipeline for single-cell RNA-seq.
Drop in any 10x Cell Ranger output, edit one YAML file, run one script —
get a clustered, annotated dataset plus structured JSON summaries that
downstream tools (or LLMs) can read directly.

> **Using [Claude Code](https://claude.com/claude-code)?** This repo ships
> with a `CLAUDE.md`. Skip the README and just ask:
> *"Set this pipeline up for the Cell Ranger output at `/path/to/run` and
> iterate the QC for me."* See [Just ask Claude Code](#just-ask-claude-code)
> below.

## What it does

Five sequential, modular steps. Each step is a standalone R script that
reads from a single config file, writes its own RDS checkpoint and a JSON
summary, and can be re-run independently.

| Step | Script | What it does | Key outputs |
|------|--------|--------------|-------------|
| 1 | `01_load_qc_metrics.R` | Load Cell Ranger H5 files, merge samples, compute per-cell QC metrics (`nFeature`, `nCount`, `percent.mt`, `percent.ribo`, complexity) | `seurat_raw.rds`, `qc_metrics.csv`, `qc_summary.json` |
| 2 | `02_soupx.R` | Estimate and remove ambient RNA contamination per sample using SoupX | `seurat_decontaminated.rds`, `soupx_summary.json` |
| 3 | `03_filter.R` | Apply QC thresholds, optional doublet detection (scDblFinder) | `seurat_filtered.rds`, `filter_summary.json` |
| 4 | `04_normalize_cluster.R` | SCTransform / LogNormalize, PCA, optional batch integration (Harmony / RPCA / CCA), UMAP, multi-resolution clustering | `seurat_clustered.rds`, `umap_coordinates.csv`, `pca_variance.csv`, `cluster_summary.json` |
| 5 | `05_annotate.R` | Optional Azimuth, optional SingleR, marker-gene scoring, `FindAllMarkers`, build per-cluster annotation table | `seurat_final.rds`, `cluster_annotations.csv`, `cluster_markers.csv`, `annotation_summary.json` |

After all steps complete, `run.sh` collects every `*_summary.json` plus the
config used into a single `pipeline_report.json` for downstream consumption.

## Why this exists

- **Drop-in for Cell Ranger output.** No glue code per project — just edit
  paths in the config.
- **Single source of truth.** Every parameter (QC thresholds, normalization
  method, integration strategy, marker genes, references) lives in one YAML.
- **AI/LLM-friendly.** Each step writes a JSON summary with a consistent
  envelope, and `run.sh` merges them into one `pipeline_report.json`. An
  LLM (or any downstream report generator) can read the entire run state
  from a single file without parsing R logs.
- **Resumable.** Each step writes an RDS checkpoint. Re-run from any step
  by invoking `Rscript scripts/0X_*.R --config config.yaml` directly.
- **Multi-sample aware.** Handles 1..N samples; integration is auto-skipped
  with a single sample.

## Requirements

- R >= 4.2 (Seurat v5 recommended)
- Required R packages: `Seurat`, `SoupX`, `yaml`, `jsonlite`, `ggplot2`
- Optional R packages (warnings only if missing):
  - `scDblFinder`, `SingleCellExperiment` — doublet detection
  - `harmony` — batch integration
  - `Azimuth` — reference-based cell type prediction
  - `SingleR`, `celldex` — alternative reference annotation
- Python 3 with `pyyaml` — used by `run.sh` to merge JSON summaries
  (the pipeline still completes without it; only `pipeline_report.json` is skipped)

## Quickstart

```bash
git clone <this-repo> my-analysis
cd my-analysis

# 1. Copy and edit the config for your data
cp config.template.yaml config.yaml
$EDITOR config.yaml

# 2. Run the full pipeline
./run.sh config.yaml
```

That's it. Outputs land in `outputs/`, logs in `logs/`.

## Pointing the pipeline at Cell Ranger output

The minimum information per sample is a unique id and the path to a Cell
Ranger H5 matrix. SoupX additionally needs the raw H5 (skipped per-sample
if missing).

**Standard `cellranger count` layout:**

```
data/
└── sample1/
    └── outs/
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
data/run1/
└── per_sample_outs/
    └── sample1/
        └── count/
            ├── sample_filtered_feature_bc_matrix.h5
            └── sample_raw_feature_bc_matrix.h5
```

```yaml
samples:
  - id: "sample1"
    filtered_h5: "data/run1/per_sample_outs/sample1/count/sample_filtered_feature_bc_matrix.h5"
    raw_h5:      "data/run1/per_sample_outs/sample1/count/sample_raw_feature_bc_matrix.h5"
```

You can add **any extra fields** under a sample (`condition`, `batch`,
`donor`, `timepoint`, ...) and they will be automatically attached as
per-cell metadata columns on the merged Seurat object.

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

## AI-friendly JSON outputs

Every step writes a `<step>_summary.json` with the same envelope:

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

Errors are also emitted as JSON to stderr:

```json
{
  "error": "Filtered H5 file not found",
  "details": { "sample": "sample1", "path": "..." },
  "timestamp": "..."
}
```

After `run.sh` completes, **`outputs/pipeline_report.json`** is the single
file an LLM (or downstream report generator) can read to understand the
entire run — every per-step summary plus the exact config used. The schema
is intentionally flat and stable.

### Letting Claude Code (or another LLM) iterate the pipeline for you

The JSON summaries aren't just for reporting — they exist so an LLM can
read them, judge whether the pipeline made good decisions, and propose
config edits *without ever loading the Seurat object into memory*. Each
summary deliberately includes the numbers a reviewer would actually want:

- **`qc_summary.json`** — per-sample cells loaded, genes detected, median
  `nFeature`/`nCount`/`percent.mt`, and the `q25/q50/q75/q95/q99` quantiles
  of `percent.mt`. An LLM can spot a too-strict (or too-loose) MT cutoff
  in one read.
- **`soupx_summary.json`** — per-sample estimated contamination fraction,
  number of cells affected, and which samples were skipped.
- **`filter_summary.json`** — cells before/after for each criterion
  (`nFeature_min`, `nFeature_max`, `nCount_min`, `nCount_max`, `percent_mt`),
  per-sample retention rates, and total doublet counts. Catches "we
  silently dropped 70% of sample X" before it pollutes downstream results.
- **`cluster_summary.json`** — PCA variance curve, suggested PC count,
  per-resolution cluster counts, per-sample cluster composition. Lets an
  LLM tell you whether `dims_use` looks right and whether one sample is
  driving an entire cluster.
- **`annotation_summary.json`** — per-method status (Azimuth, SingleR),
  median Azimuth prediction score, low-confidence cell %, top markers per
  cluster, and the suggested label table. The honest "is this annotation
  trustworthy?" signal.

A typical **assisted** iteration looks like:

```
1. ./run.sh                                      # initial run with defaults
2. Ask Claude Code: "read outputs/pipeline_report.json and tell me which
   QC thresholds need adjusting and why"
3. Apply Claude's suggested edits to config.yaml
4. Re-run only the affected steps:
     Rscript scripts/03_filter.R --config config.yaml
     Rscript scripts/04_normalize_cluster.R --config config.yaml
     Rscript scripts/05_annotate.R --config config.yaml
5. Repeat until the LLM agrees the run looks clean.
```

Because the config is a single YAML and every step is resumable, the
iteration loop is fast and the LLM has full context for *why* each
parameter was set the way it was — the config used is embedded in
`pipeline_report.json` itself.

### Full autonomous mode

If you'd rather not babysit the loop at all, hand the whole thing to
Claude Code as a single instruction:

> *"Run this pipeline on the Cell Ranger output at `/path/to/run` in
> autonomous mode. Iterate the config on your own using the JSON
> summaries until the QC looks clean, log every parameter you tried in
> `logs/iteration_log.md`, and just report to me at the end."*

Claude will:

1. Wire `config.yaml` to your data and run the pipeline.
2. Read `outputs/pipeline_report.json`, judge each step against the QC
   checklist in `CLAUDE.md`, and decide what to change.
3. Edit `config.yaml`, re-run only the affected steps, and append a new
   entry to `logs/iteration_log.md` with: iteration number, what was
   wrong, what changed, and the resulting key metrics.
4. Repeat until the run converges (or until it hits the iteration cap
   it sets for itself — usually 4–5 rounds).
5. Stop and write a short summary report covering: final cell counts,
   cluster count, integration status, annotation confidence, anything
   that still looks suspicious, and the full iteration trail.

The iteration log is a normal markdown file you can read at the end to
audit every decision Claude made and roll back any change you don't
agree with — the previous configs and JSON summaries are all on disk.
Treat it like a lab notebook written by a research assistant who walked
the data through QC for you.

Autonomous mode is opt-in: Claude will only run the pipeline without
asking each time if you say so explicitly in the first prompt. Otherwise
it falls back to the assisted loop above and confirms each re-run.

## Just ask Claude Code

This repo ships with a `CLAUDE.md` at the root, so if you're using
[Claude Code](https://claude.com/claude-code) (or any agent that respects
`CLAUDE.md`) you can skip most of this README and just ask. Reasonable
opening prompts:

- *"Set up this pipeline for the Cell Ranger output at `/path/to/run`."*
- *"My data is mouse PBMCs — update the config and marker genes accordingly."*
- *"Run the pipeline, then read `outputs/pipeline_report.json` and tell me
  what to adjust."*
- *"Sample MAD-03 lost 60% of cells in filtering — figure out which
  threshold caused it and propose a fix."*
- *"Add a new sample to the config and re-run only the steps that need
  it."*

Claude will read `CLAUDE.md`, the config, and any existing outputs, and
take it from there.

## Configuration reference

`config.template.yaml` is the canonical reference and is heavily commented.
Top-level sections:

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
| `integration` | Batch correction method (Harmony / RPCA / CCA) |
| `annotation` | Azimuth ref, SingleR ref, custom marker gene dict |

### Customizing for tissue / species

- **Mouse:** set `qc.mt_pattern: "^mt-"`.
- **Tissue-specific Azimuth references:** see
  https://azimuth.hubmapconsortium.org/references/ — set
  `annotation.azimuth_reference` to e.g. `lungref`, `pbmcref`, `kidneyref`,
  `heartref`, `bonemarrowref`, `tonsilref`, `mousecortexref`.
- **Custom markers:** edit `annotation.marker_genes` — top-level groups are
  just for organization, the leaf entries are `cell_type: [gene, gene, ...]`.
  Each cell type gets an `AddModuleScore` and the highest-scoring type per
  cluster becomes that cluster's marker-based label.
- **10x Flex (probe-based) data:** loosen QC drastically — Flex produces
  much sparser counts. Try `nFeature_min: 20`, `nFeature_max: 5000`,
  `nCount_min: 30`, `nCount_max: 25000`, `percent_mt_max: 10`. Only ~12 MT
  genes are in the probe panel, so `percent.mt` is a weak signal.

## Resuming / re-running individual steps

Each script reads its input from the previous step's RDS, so you can re-run
any step alone after editing the config:

```bash
# Re-run only normalization + clustering after tweaking parameters
Rscript scripts/04_normalize_cluster.R --config config.yaml

# Re-run only annotation after adding new marker genes
Rscript scripts/05_annotate.R --config config.yaml
```

Pass `--no-plots` to skip PNG generation for faster iteration.

## Project layout

```
seurat-pipeline/
├── README.md
├── config.template.yaml
├── run.sh
├── scripts/
│   ├── utils.R
│   ├── 01_load_qc_metrics.R
│   ├── 02_soupx.R
│   ├── 03_filter.R
│   ├── 04_normalize_cluster.R
│   └── 05_annotate.R
├── outputs/         # generated, gitignored
└── logs/            # generated, gitignored
```

## License

MIT.
