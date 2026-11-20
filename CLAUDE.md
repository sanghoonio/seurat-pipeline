# CLAUDE.md

Guidance for Claude Code (and other LLM agents) working in this repo.

## What this repo is

A config-driven, 8-step Seurat v5 pipeline for single-cell RNA-seq. One
`config.yaml` holds shared settings (project, paths, samples) and one **run
block** per run of each step. Each step is a standalone R script that reads the
output of one run of the previous step and writes its own folder,
`outputs/<step>/<label>/`, with an RDS checkpoint, CSVs, plots and a JSON summary.

The pipeline is designed to be **drop-in for 10x Cell Ranger output** and
**driven by an LLM in a QC iteration loop**. Most requests will be:

1. Wire the config up to a dataset.
2. Run a step (one run at a time).
3. Read its summary and plots and recommend config edits.
4. Rerun that step, or add a labeled variant, and re-evaluate.

## Repo layout

```
seurat-pipeline/
├── README.md                  # user-facing docs
├── CLAUDE.md                  # this file
├── config.template.yaml       # reference for every key, in run-block form
├── run.sh                     # ./run.sh <config> --report-only -> pipeline_report.json
├── scripts/
│   ├── utils.R                # STEPS registry, CONFIG_DEFAULTS, run-block parsing/validation,
│   │                          #   load_config() (run selection, lineage, effective settings),
│   │                          #   ancestor_dir(), read_ancestor_rds(), save_run_rds(), summaries,
│   │                          #   matrix readers, emptyDrops, scVI and CellTypist bridges
│   ├── 00_cell_calling.R      # optional: emptyDrops on raw -> called barcodes
│   ├── 01_load_qc_metrics.R   # load + merge + QC metrics
│   ├── 02_ambient.R           # SoupX or DecontX per sample
│   ├── 03_doublets.R          # scDblFinder per sample, calls only
│   ├── 04_filter.R            # thresholds + doublet removal; removed cells kept
│   ├── 05_normalize.R         # SCTransform / LogNormalize + PCA
│   ├── 06_cluster.R           # integration, UMAP, clustering, FindAllMarkers
│   ├── 07_annotate.R          # Azimuth / SingleR / CellTypist / custom reference / markers
│   ├── report.R               # run index + staleness; --list-labels <step>
│   ├── scvi/                  # bundled Python env + scVI/sysVI (see scripts/scvi/README.md)
│   ├── celltypist/            # CellTypist annotation, called by step 07
│   └── refs/                  # one-time helpers to build a custom Seurat reference
├── outputs/                   # gitignored: outputs/<step>/<label>/, pipeline_report.json
└── logs/                      # gitignored
```

No package, no build step. Scripts run with `Rscript` from the repo root.

## Runs, labels and inputs

- **Run blocks:** `<step>:` is the step's unlabeled run (label `default`);
  `<step>/<label>:` is a labeled run. A step can have any number of runs.
- **Self-contained:** a block holds every setting its run uses. Absent keys take
  `CONFIG_DEFAULTS` (explicit nulls stay null); nothing is inherited from other
  blocks. A block may contain only the keys its step owns:

  | Step | Owns | `input` names a run of |
  |---|---|---|
  | 00_cell_calling | `cell_calling.emptydrops` | none |
  | 01_load_qc_metrics | `cell_calling.cells`, `qc.metrics` | 00 (only when cells is emptydrops/union) |
  | 02_ambient | `ambient` | 01 |
  | 03_doublets | `doublets` | 02 |
  | 04_filter | `qc.thresholds`, `qc.prefilter_umap`, `filter` | 03 |
  | 05_normalize | `normalization`, `reduction.npcs` | 04 |
  | 06_cluster | `integration`, `reduction.dims_use`, `reduction.umap`, `clustering`, `markers` | 05 |
  | 07_annotate | `annotation` | 06 |

- **`input`:** a quoted label; absent = the previous step's unlabeled run, which
  must then exist. A step only runs if it has a block.
- **Effective settings:** a run uses built-in defaults, then every ancestor run's
  *recorded* settings (from their summaries, in chain order), then its own block.
  So step 02 uses the MT pattern its 01 run used, and step 06 knows how its 05 run
  normalized, regardless of what other blocks say.
- **Ancestor files:** steps reach older runs through the lineage, not config:
  04 reads 02's object and 03's calls; 07 reads 06's object and markers, 04's
  removed cells and 02's per-gene removal table.
- **Safety:** a step refuses an input run without a `success`/`warning` summary
  (never run, still running, or failed), and refuses an RDS whose `run_id` stamp
  does not match its run's summary. Starting a run clears its folder.
- **Staleness:** a run is stale when an ancestor was rerun after it was made, or
  its block's settings or `input` changed since. `./run.sh <config> --report-only`
  (or `Rscript scripts/report.R --config <config>`) lists every run as `success`,
  `warning`, `failed`, `missing`, `stale` (with reasons) or `not_in_config`, and
  writes `outputs/pipeline_report.json`. A step's own `warnings` say when its
  input is stale.

## Standard workflow for a new Cell Ranger output

1. **Find the matrices.** `cellranger count`: `<sample>/outs/filtered_feature_bc_matrix.h5`
   and `raw_feature_bc_matrix.h5`. `cellranger multi` / Flex:
   `per_sample_outs/<id>/[count/]sample_filtered_feature_bc_matrix.h5` (and `_raw_`).
   Glob `**/*filtered_feature_bc_matrix*`. Older runs may only have mtx
   directories; `filtered`/`raw` accept a directory.
2. **Copy `config.template.yaml` to `config.yaml`**, fill `samples`, set
   `project.name`. Leave the step blocks at their defaults for a first pass.
3. **Run step 01**, `Rscript scripts/01_load_qc_metrics.R --config config.yaml`,
   only after telling the user (steps write GB-scale RDS files).
4. Read the run's summary (`outputs/01_load_qc_metrics/default/qc_summary.json`)
   and plots, walk the checklist below, propose edits, then the next step.

## QC review checklist

Read each run's summary, in its folder. Every summary shares one envelope:
`step`, `timestamp`, `pipeline_commit` (`-dirty` = uncommitted changes),
`config.path`/`md5`, `versions`, `run` (step, label, run_id, settings_hash,
settings), `input` (the run it read), `lineage` (every ancestor run with its
run_id and settings), `inputs`, `params`, `outputs` (paths relative to
`paths.output_dir`), `warnings`, and `status`: `success`, `warning` (read the
warnings) or `failed` (`error`/`details` say why). A run killed outright (out of
memory, time limit) leaves no summary; the report shows it as `missing`.

### 00 cell calling (`cell_calling_summary.json`)
- `per_sample[sid].n_emptydrops_only`: cells emptyDrops recovered beyond Cell
  Ranger's set. Zero means nothing was recovered (raise `fdr`, lower `lower`).
- `emptydrops_only_stats` (median nCount/nFeature) says whether the recovered
  cells look like real low-RNA cells or junk. `barcode_rank` has knee/inflection.
- Raw matrices contribute barcodes only; their features are restricted to the
  filtered matrix's (Flex raw matrices keep excluded probes).
- **Coupling:** recovered cells are low-UMI, and 04's `nCount_min`/`nFeature_min`
  cut on the same axis. When 01 uses emptydrops/union, lower them and confirm in
  `filter_summary.json` `per_sample[sid].n_cells_by_source` that the recovered
  cells survive. A cell recovered in 00 and removed in 04 is still lost.

### 01 loading (`qc_summary.json`)
- `n_samples_loaded` equals `n_samples_configured`; flag samples much smaller or
  larger than the rest (`per_sample[sid].n_cells`).
- `merged.quantiles` / `per_sample[sid].quantiles` (min, q01 … q99, max of nCount,
  nFeature, percent.mt, percent.ribo, complexity) inform 04's thresholds. These
  are uncorrected counts; 04 filters the corrected counts, whose quantiles are in
  `ambient_summary.json`.
- On 10x Flex, few MT genes are in the probe panel, so MT quantiles are low; do
  not recommend loosening the MT cutoff for that reason, document it instead.

### 02 ambient (`ambient_summary.json`)
- `per_sample[sid].estimated_contamination`: above ~0.2 is high, above ~0.4 suspect.
  `contamination_basis` says what it is (DecontX: median of per-cell values;
  `contamination_mean`/`_q90` give the spread; each cell's value is
  `decontX_contamination`).
- Per-sample `status` `skipped`/`failed`: find out why (missing `raw`, listed in
  `ambient.skip_samples`, estimation failed).
- SoupX `at_contamination_bound` other than `none`: `upper` is usually an
  unreliable fit, `lower` a near no-op. On Flex autoEstCont can fail or pin;
  skip those samples or use DecontX.
- `top_removed_genes` should look like ambient transcripts (highly expressed,
  cell-type restricted), not housekeeping genes; `umis_removed_pct` should track
  contamination. `ambient_removed_by_gene.csv` has every gene with its
  `ambient_profile_frac`.
- `02_decontx_contamination_vs_nFeature.png` should show contamination rising as
  genes per cell fall.
- `merged.quantiles` are the corrected counts 04 filters.
- `corrected: false` means a passthrough run (`ambient.enabled: false`).

### 03 doublets (`doublets_summary.json`)
- `per_sample[sid].pct_doublets` against `params.expected_rate`.
  `doublet_detection` is `partial`/`failed` when samples were skipped (too few
  cells) or failed; their cells are all recorded as singlets.

### 04 filter (`filter_summary.json`)
- `total.pct_removed` over 30% is worth a comment; over 50% adds a warning.
- `removed_by_threshold` (keys match `qc.thresholds`; a cell can fail several,
  `any` is the union) says which threshold removed the most. A huge
  `nFeature_min` count usually means the bound is too high for the assay.
- `per_sample[sid].pct_removed` far above the others: thresholds wrong for that
  sample, or a real quality problem.
- `n_removed_only_after_correction`: cells failing on corrected counts that pass
  on uncorrected totals. A large number means ambient correction, not quality,
  removes them; check what they are before tightening anything.
- `n_removed_doublet_only`: doublets that passed every threshold.
- Removed cells are in `seurat_filtered_out.rds` / `filtered_out_cells.csv` with
  `removal_reason` (`thresholds`, `doublet`, `thresholds+doublet`) and
  `failed_<threshold>`.
- Plots: violins and scatter with this run's thresholds, cells removed per
  sample, and with `qc.prefilter_umap` all input cells coloured by kept / reason.

### 05 normalize (`normalize_summary.json`)
- `pca.cumulative_variance_pct`: shares of total variance across the variable
  features, so well under 100% is normal. `05_elbow_plot.png` informs 06's
  `dims_use`.

### 06 cluster (`cluster_summary.json`)
- `integration.status`: `applied`, `disabled` or `skipped` (single-level
  `batch_var`); `reduction_used` is the embedding behind neighbors/UMAP/clusters.
  `batch_var` is any metadata column (e.g. `library`, `chemistry`).
- RPCA/CCA need the 05 run to have `normalization.split_by` equal to `batch_var`.
- scVI is for large cross-technology gaps (e.g. probe panel vs whole
  transcriptome) that Harmony leaves split; judge it by whether shared cell types
  mix across batches while rare or confounded populations survive, not by a
  uniformly blended UMAP. See `scripts/scvi/README.md`. Its cache lives in the 05
  run's folder, keyed on the 05 run_id and scVI params.
- `clustering.resolutions`: clusters per resolution; raise `default_resolution`
  if far fewer clusters than expected cell types.
- `per_cluster[cl].top_sample_pct` far above `top_sample_expected_pct`: a
  donor-dominated cluster (`cluster_counts_by_sample` has the table).
- `06_elbow_plot.png` marks `dims_use`; `cluster_markers.n_markers` should be > 0.

### 07 annotate (`annotation_summary.json`)
- `methods.azimuth.median_prediction_score` below ~0.6 means a poor reference
  match; `low_confidence_pct` high means Azimuth is guessing.
- Each other method has its own block: `methods.singler.references[ref]`,
  `methods.celltypist.models[model]`, `methods.custom_reference`. Each reports
  `gene_coverage_pct` (reference genes present in the data); low coverage (probe
  panels miss genes) weakens that method. Labels are metadata columns
  (`singler_<ref>`, `celltypist_<model>`, `customref_<name>_<col>`) and CSVs in
  `annotation/`. The pipeline does not compare methods; that is for the user.
- `cluster_annotations[]`: `suggested_label` and `label_source` (azimuth /
  markers / markers (override) / singler / none), `azimuth_label_fraction` (low =
  mixed cluster), Azimuth score and low-confidence %, `marker_score` /
  `marker_margin` (near-zero margin = arbitrary pick). Many `Unknown` clusters
  mean the references and marker panel miss relevant types.
- Reference levels that say "None" (lungref smooth muscle at `ann_level_3`) are
  filled from the nearest coarser level (`n_cells_label_from_coarser_level`).
- `removed_cells` (`annotation.removed_cells: true`): per removal reason, Azimuth
  labels, score and best marker type of the 04 run's removed cells (under 10
  UMIs listed but not annotated). Removed cells that confidently map to a real
  type mean the thresholds cut that type.
- `compare_by[col]`: per level, cells, samples, score, low-confidence %, label
  composition, per-label score; `labels_with_score_gap` flags labels whose score
  differs by >= `params.compare_score_gap` between levels. With few samples per
  level, differences can be donor effects.
- `marker_genes_removed_by_ambient_pct` (needs a marker panel): per sample and
  marker type, the share of marker UMIs the 02 run removed. High = correction is
  taking that type's own signal.

When a problem shows up, **propose the smallest config edit that fixes it** and
rerun only the affected run. To compare settings without losing the current
result, add a labeled block instead of editing the existing one. Downstream runs
of a rerun become stale (the report says which) and must be rerun to pick it up.

## Running steps

```bash
Rscript scripts/04_filter.R --config config.yaml                       # the step's only run
Rscript scripts/04_filter.R --config config.yaml --label dx_floor100   # a chosen run
Rscript scripts/04_filter.R --config config.yaml --label dx_floor100 --no-plots
./run.sh config.yaml --report-only
```

`--label` is required when a step has several runs (the step stops and lists
them). Unknown arguments are errors.

## Autonomous mode

If the user explicitly asks you to run the pipeline "autonomously", "end to
end", "without asking each time" or "iterate it until it's clean", you may run
step scripts without confirming each one. Otherwise, ask before each run.

In autonomous mode:

1. Wire the config (samples, `project.name`, blocks at defaults).
2. Run the steps in order, one run each, reading each summary before the next.
3. Walk the checklist; decide the smallest fix.
4. Edit the config. Prefer a new labeled block when the change is an
   experiment; edit in place when it is a correction.
5. Rerun only the affected runs and their descendants (the report's `stale`
   entries list them). Use `--no-plots` on intermediate iterations.
6. Append to `logs/iteration_log.md` (create it on the first iteration):

   ```markdown
   ## Iteration N — YYYY-MM-DD HH:MM
   **Trigger:** what looked wrong (e.g. "sample_03 lost 64% of cells, almost all to nFeature_min")

   **Change:**
   - `04_filter/default qc.thresholds.nFeature_min`: 200 -> 80

   **Runs:** 04_filter/default, 05_normalize/default, 06_cluster/default, 07_annotate/default

   **Result:**
   - Total cells: 12,403 -> 14,891 (+20%)
   - sample_03 retention: 36% -> 78%
   - Clusters at res 0.5: 14 -> 16
   ```

7. Stop when every checklist item is fine, after 5 iterations without
   meaningful improvement, before any change you can't justify from the
   summaries, or when a step fails for a reason that needs the user.
8. Report: final cell counts (total and per sample), clusters and resolution,
   integration status, annotation confidence, anything still suspicious, and a
   pointer to `logs/iteration_log.md`.

Hard rules: never delete earlier blocks or configs to "clean up" (the log and the
run folders are the record); never widen thresholds just to silence a warning;
never run on data you cannot find.

## Customizing the config

- **Mouse:** `01_load_qc_metrics` block, `qc.metrics.mt_pattern: "^mt-"`.
- **Thresholds:** `04_filter` block; any bound can be null. Tune from 01's and
  02's quantiles.
- **10x Flex / probe-based:** 3'/5' defaults may not transfer; percent.mt is a
  weak signal; some genes are absent from the panel (check marker genes).
- **New tissue:** `07_annotate` block, `annotation.markers.genes` as compartment ->
  cell type -> genes (compartments are for readability only).
- **Cross-technology integration:** `06_cluster` block, `integration.method: scvi`
  and `batch_var` (e.g. `chemistry`); one-time `cd scripts/scvi && uv sync`.
- **Azimuth:** `annotation.azimuth.reference` (`lungref`, `pbmcref`, ...) and
  `primary_level` (e.g. `ann_level_3`); a level the reference lacks is an error.
- **Other annotation:** `annotation.singler.reference` (one or more celldex
  references), `annotation.celltypist.models` (Python env),
  `annotation.custom_reference` (an RDS built with `scripts/refs/`).
- **Sample metadata:** extra fields under a sample become metadata columns; for
  many fields use `sample_metadata.file` (a CSV joined on sample id). Summaries
  record its columns, never values.
- **Defaults:** `CONFIG_DEFAULTS` in `utils.R` must stay in sync with
  `config.template.yaml`. Ownership lives in `STEPS`.

## Things not to do

- **Do not delete `outputs/` or `logs/`** without asking; run folders are the
  record of every variant.
- **Do not hand-edit RDS files or mutate checkpoints in ad-hoc R sessions.** Change
  the config and rerun the run.
- **Do not invent new steps.** The chain is a fixed contract; propose any change
  to it first. Sanctioned helpers that are not steps: `scripts/scvi/` (called by
  06), `scripts/celltypist/` (called by 07), `scripts/refs/` (run by hand),
  `scripts/report.R` (the report).
- **Do not add R dependencies to `utils.R`** beyond `yaml`, `jsonlite`, `ggplot2`,
  `Matrix`, `Seurat` and the optional packages above. The Python bridges shell
  out via `system2` and link no R package.
- **Missing deps for a requested feature hard-stop, they don't skip.** Don't
  "fix" an error by turning the feature off; install the dependency or confirm
  with the user.
- **Do not commit `config.yaml`, `outputs/` or `logs/`.**
- **Do not run a step without telling the user first** (outside autonomous mode).

## Known issues

- **Only RDS files and summaries are written atomically.** CSVs are written in
  place. A step that starts while its input run is being rerun is refused (no
  summary yet), but a long step that reads an ancestor's files late (07 reading
  04's removed cells) while that ancestor is rerunning can read a partial file.
  Don't rerun an ancestor while its descendants are running.
- **Combining samples measured with different probe panels** (e.g. Flex v1 and
  v2) merges their gene sets in step 01; genes only in one panel read as zero in
  the other. There is no option yet to keep only shared genes.
