# CLAUDE.md

Guidance for Claude Code (and other LLM agents) working in this repo.

## What this repo is

A config-driven, 5-step Seurat v5 pipeline for single-cell RNA-seq. A
single `config.yaml` controls every parameter; each step is a standalone
R script that reads a checkpoint RDS from the previous step, writes its
own RDS, and emits a JSON summary. After all steps run, `run.sh` merges
every `*_summary.json` plus the config used into one `outputs/pipeline_report.json`.

The pipeline is designed to be **drop-in for 10x Cell Ranger output** and
**driven by an LLM in a tight QC iteration loop** — that's the primary
use case, not "Claude writes new R code." Most user requests will be:

1. Wire the config up to a new dataset.
2. Run the pipeline.
3. Read the JSON summaries and recommend config edits.
4. Re-run the affected steps and re-evaluate.

## Repo layout

```
seurat-pipeline/
├── README.md                  # user-facing docs (read this first if you need depth)
├── CLAUDE.md                  # this file
├── config.template.yaml       # canonical reference for every config field
├── run.sh                     # runs all steps + builds pipeline_report.json
├── scripts/
│   ├── utils.R                # config loading, JSON writing, %||%, attach_sample_metadata(),
│   │                          #   read_10x_any(), call_cells_emptydrops(), read_called_barcodes(),
│   │                          #   run_scvi_integration() (shells out to scripts/scvi/)
│   ├── 00_cell_calling.R      # (opt-in) emptyDrops cell-calling from RAW → called barcodes
│   ├── 01_load_qc_metrics.R   # load matrix (h5 or mtx dir) → merge → QC metrics
│   ├── 02_ambient.R           # ambient RNA correction, SoupX or DecontX (per-sample)
│   ├── 03_filter.R            # QC thresholds + scDblFinder
│   ├── 04_normalize_cluster.R # SCTransform/PCA/UMAP/clustering; integration: harmony|rpca|cca|scvi
│   ├── 05_annotate.R          # Azimuth/SingleR/marker scoring/FindAllMarkers
│   └── scvi/                  # OPT-IN scVI/sysVI Python subsystem (see scripts/scvi/README.md)
├── outputs/                   # gitignored — RDS checkpoints, CSVs, JSONs, plots/
└── logs/                      # gitignored — per-step stdout/stderr per run
```

There is no `src/`, no package, no build step. Scripts are run directly
with `Rscript` from the repo root.

## Standard workflow when a user gives you a Cell Ranger output

1. **Find the matrices.** A standard `cellranger count` run has them at
   `<sample>/outs/filtered_feature_bc_matrix.h5` and
   `<sample>/outs/raw_feature_bc_matrix.h5`. A `cellranger multi` / Flex
   run nests them under `per_sample_outs/<id>/count/sample_filtered_feature_bc_matrix.h5`
   (and `_raw_`). Use `Glob` with patterns like `**/*filtered_feature_bc_matrix.h5`.
   Older runs may have only the mtx **directories** (`filtered_feature_bc_matrix/`
   with barcodes/features/matrix.mtx.gz) and no `.h5` — that's fine, the
   `filtered`/`raw` config fields accept a directory path too.
2. **Copy `config.template.yaml` → `config.yaml`** if it doesn't exist
   yet, then edit the `samples:` block to point at the discovered H5
   files. Set `project.name` to something descriptive. Leave QC defaults
   alone on the first pass — you'll tune them after the first run.
3. **Run `./run.sh config.yaml`.** Don't run without confirming with the
   user first; this can take many minutes and writes large files.
4. **Read `outputs/pipeline_report.json`** and walk through the QC iteration
   checklist below.

## QC iteration checklist (the AI loop this repo is designed for)

After every run, read `outputs/pipeline_report.json` (or the individual
`*_summary.json` files) and check, in order. After running steps one at a time,
rebuild the report with `./run.sh <config> --report-only`; its `_index` lists
each step's `status`, `timestamp`, `pipeline_commit` and `n_warnings`, so a
stale or failed step is obvious.

Every summary shares one envelope: `step`, `timestamp`, `pipeline_commit`
(`-dirty` = uncommitted changes), `config.path`/`config.md5`, `versions`,
`inputs`, `params` (the values actually used), `outputs` (paths relative to
`output_dir`), `warnings`, and `status`:
- `success`, `warning` (succeeded but `warnings` is non-empty: read them),
  `skipped` (step disabled), or `failed` (the step errored; `error`/`details`
  say why). A step deletes its previous summary as soon as its config is
  read, and any error after that (config validation, `emit_error`, or an
  unhandled R error) writes a `failed` summary instead. A job killed outright
  (out of memory, time limit) leaves no summary for that step, which the
  report's `_index` shows as `missing`.

### 0. Cell calling (`cell_calling` section, from `00_cell_calling`)
- Only present when `cell_calling.cells` is `emptydrops` or `union`. Check
  `per_sample[sid].n_emptydrops_only` — how many cells emptyDrops recovered
  beyond the Cell Ranger filtered set. Zero means the params recovered nothing
  (raise `fdr`, lower `lower`).
- `per_sample[sid].emptydrops_only_stats` (median nCount/nFeature/percent.mt of the
  recovered fraction) tells you whether you caught real low-RNA cells or junk.
  Use `<id>_barcode_ranks.csv` + the `barcode_rank` knee/inflection for a
  barcode-rank plot.
- Raw matrices contribute **barcodes only**: steps 00/01/02 restrict raw features to
  the filtered matrix's feature set (`restrict_to_filtered_features()`), because
  Cell Ranger's raw matrix keeps excluded probes (Flex: ~21k, incl. `DEPRECATED_*`
  duplicates of valid genes). The log line `features: raw=… -> keeping …` records it.
- **CRITICAL coupling:** emptyDrops recovers low-UMI cells; the step-03
  `nFeature_min`/`nCount_min` thresholds cut on the same axis and will silently
  delete them. When cell-calling is on, lower those thresholds (e.g.
  `nFeature_min ~50-100`, `nCount_min ~100`) and confirm via
  `filter_summary.json` `per_sample[sid].n_cells_by_source` (before/after)
  that the recovered fraction survives. Do NOT judge
  recovery by step 00 alone — a cell recovered in 00 and deleted in 03 is
  still lost.

### 1. Loading (`qc` section, from `01_load_qc_metrics`)
- Did every sample load? `n_samples_loaded` should equal
  `n_samples_configured`; `samples_loaded` lists them.
- Are any samples drastically smaller/larger than the rest
  (`per_sample[sid].n_cells`)? Note it.
- Look at `merged.quantiles` (and `per_sample[sid].quantiles`): min/q01/q05/
  .../q99/max of nCount_RNA, nFeature_RNA, percent.mt, percent.ribo,
  complexity. If `percent.mt` `q95` is, say, 25% but the
  config has `percent_mt_max: 10`, you're about to lose the entire
  high-MT tail — flag it before step 03 runs. For 10x Flex data, MT
  quantiles are usually low because few MT genes are in the panel; do
  NOT recommend loosening the MT cutoff in that case, recommend
  documenting it instead.
- `threshold_preview_pct` (per sample and merged) shows what each
  `qc.thresholds` bound would remove in step 03 (`any` = union). Flag any
  single threshold removing a large share of a sample. Step 02 repeats the
  preview on corrected counts, which is what step 03 actually filters.

### 2. Ambient RNA (`ambient` section, from `02_ambient`)
- `per_sample[sid].estimated_contamination` — anything above ~0.2 is
  high; above ~0.4 is suspect. Mention it.
- Any samples with `status: "skipped"` or `"failed"`? Find out why
  (usually missing `raw`, or the sample is listed in `ambient.skip_samples`).
- `at_contamination_bound` other than `"none"` means the estimate sits on a
  `contamination_range` bound: `upper` is usually an unreliable fit, `lower`
  a near no-op. On Flex, autoEstCont sometimes finds no soup markers and
  fails or pins at a bound; it can also work fine. Skip those samples.
- `top_removed_genes` should look like ambient transcripts (highly expressed,
  cell-type-restricted genes), not housekeeping genes across the board.
  `umis_removed_pct` should track `estimated_contamination`.
- `params.method` is `soupx` or `decontx`; `contamination_basis` says what
  `estimated_contamination` is. With DecontX it's the median per-cell fraction
  (`contamination_mean`/`_q90` give the spread)
  and each cell's value is kept as `decontX_contamination` in the metadata.
  `02_decontx_contamination_vs_nFeature.png` should show contamination rising
  as genes per cell fall; a flat cloud means the per-cell model added little
  over SoupX.

### 3. Filtering (`filter` section)
- `total.pct_removed` — over 30% is worth a comment, over 50% adds a
  warning automatically; investigate.
- `total.removed_by_threshold` / `per_sample[sid].removed_by_threshold` say
  which threshold dropped the most cells (keys match `qc.thresholds`; a cell
  can fail several, `any` is the union). If `nFeature_min` is dropping huge
  numbers, the threshold is probably too high for the assay (especially Flex).
- `per_sample[sid].pct_removed` — if one sample lost 60% and the others
  lost 5%, the thresholds are wrong for that sample's distribution OR
  that sample has a real quality problem.
- `per_sample[sid].doublets.pct_of_tested` — sanity-check against
  `params.doublets.expected_rate`. `doublets.status` is `partial`/`failed`
  when samples were skipped or failed (their cells were all kept).

### 4. Clustering (`cluster` section)
- PCA variance shares are of total variance across the variable features,
  so `pca.cumulative_variance_at_dims_used` well under 100% is normal.
  Use `04_elbow_plot.png` to check the curve has levelled off before
  `dims_use`.
- `clustering.resolutions` shows the cluster count at each resolution.
  If 0.5 gives 4 clusters and the user expects ~15 cell types, bump
  `default_resolution` up.
- `per_cluster[cl]` — `top_sample_pct` far above `top_sample_expected_pct`
  means a donor-dominated cluster (`cluster_counts_by_sample` has the full
  table). With integration off, that may be real biology or a batch effect;
  annotation usually tells which.
- `integration.status` is `applied`, `disabled` or `skipped` (with `reason`,
  e.g. a single-level `batch_var`); `reduction_used` is the embedding that
  neighbors/UMAP/clustering used.
- `integration.batch_var` controls what Harmony/RPCA/CCA/scVI correct over (default
  `sample_id`). For mixed chemistries/assays in one object, set it to a metadata
  column like `chemistry` to correct the assay effect.
- `integration.method: scvi` is the OPT-IN escalation for LARGE cross-technology gaps
  (e.g. Flex probe-panel vs 3′ whole-transcriptome) where harmony/rpca leave clusters
  split by batch. It needs the bundled Python env (`cd scripts/scvi && uv sync`); see
  `scripts/scvi/README.md`. Judge it by whether SHARED cell types mix across batch while
  rare/confounded populations survive — not by a uniformly-blended UMAP.

### 5. Annotation (`annotation` section)
- `methods.azimuth.median_prediction_score` — below ~0.6 means the
  reference is a poor match; consider a different `annotation.azimuth.reference`.
- `methods.azimuth.low_confidence_pct` — high values mean Azimuth is
  guessing; the marker-based labels may be more trustworthy.
- `cluster_annotations[]` — per cluster: `suggested_label` and its
  `label_source` (azimuth / markers / markers (override) / singler / none),
  plus `azimuth_label_fraction` (share of the cluster with that label; low =
  mixed cluster), `azimuth_median_score`, `azimuth_low_conf_pct`, and for
  markers `marker_score` / `marker_margin` (near-zero margin = arbitrary
  pick). Many `Unknown` clusters mean the references and
  `annotation.markers.genes` are missing relevant cell types.
- Reference levels that say "None" (e.g. lungref smooth muscle at
  `ann_level_3`) are filled from the nearest coarser level;
  `methods.azimuth.n_cells_label_from_coarser_level` counts them.

When you find a problem, **propose the smallest config edit that fixes
it** and **only re-run the steps it affects** — don't rerun the whole
pipeline if you only changed clustering parameters. The dependency chain
is: 00 → 01 → 02 → 03 → 04 → 05. Step 00 (cell calling) is opt-in and writes
barcode sets read by 01/02; 01→05 each read the previous step's RDS checkpoint.
Re-run 00 only when you change `cell_calling` params (it's the expensive
emptyDrops step); switching `cell_calling.cells` between `emptydrops` and
`union` alone needs 01/02 re-run but not 00.

## Autonomous mode

If the user explicitly asks you to run the pipeline "autonomously,"
"end-to-end," "without asking each time," or "iterate it for me until
it's clean," you have standing permission to run `./run.sh` and the
individual `Rscript scripts/0X_*.R` commands without confirming each
invocation. Otherwise, default to the assisted loop and ask before each
run.

In autonomous mode, operate as a research assistant doing QC iteration:

1. **Wire the config** to the user's data (find H5 files, populate
   `samples:`, set `project.name`).
2. **Initial run.** `./run.sh config.yaml`.
3. **Read `outputs/pipeline_report.json`** and walk the QC checklist
   above. Decide what's wrong and what the smallest fix is.
4. **Edit `config.yaml`** with the proposed fix. Edit the file, do not
   regenerate it from scratch.
5. **Re-run only the affected steps.** Use the dependency chain to
   minimize work — if you only changed `clustering.resolutions`, just
   rerun `04_normalize_cluster.R` and `05_annotate.R`. Pass `--no-plots`
   on intermediate iterations to save time.
6. **Append to `logs/iteration_log.md`** (create it on the first
   iteration). One entry per iteration, in this format:

   ```markdown
   ## Iteration N — YYYY-MM-DD HH:MM
   **Trigger:** what looked wrong in the previous run
   (e.g. "sample_03 lost 64% of cells, almost all to nFeature_min")

   **Change:**
   - `qc.thresholds.nFeature_min`: 200 → 80
   - `qc.thresholds.percent_mt_max`: 15 → 20

   **Steps re-run:** 03, 04, 05

   **Result:**
   - Total cells: 12,403 → 14,891 (+20%)
   - sample_03 retention: 36% → 78%
   - Clusters at res 0.5: 14 → 16
   - Azimuth median score: 0.71 → 0.73
   ```

   Keep entries terse but specific — numbers, not adjectives.
7. **Stop conditions.** Halt iteration when ANY of:
   - All checklist items are green (no `pct_removed > 30%`, integration
     succeeded, Azimuth median score > 0.6, no all-`Unknown` clusters,
     PCA cumulative variance reasonable, etc.).
   - You've made 5 iterations without meaningfully improving any
     metric (don't churn).
   - You're about to make a change you can't justify from the JSON
     alone (escalate to the user instead of guessing).
   - A step errors out for a reason that needs human input
     (missing reference, missing package, ambiguous sample metadata).
8. **Write a final report** as your message to the user. Cover, in this
   order: final cell counts (total + per sample), final cluster count
   and resolution, integration status, annotation confidence summary,
   anything that still looks suspicious, and a one-line pointer to
   `logs/iteration_log.md` for the full trail. Do not paste the
   iteration log into the message — link to it.

Hard rules in autonomous mode:

- **Never delete or rewrite previous configs.** If you want to compare,
  the iteration log records what each iteration used.
- **Never iterate past step 5 without a passing reason.** Random walks
  in parameter space are not iteration.
- **Never silently widen thresholds to "make the warning go away."**
  If `pct_removed > 50%` is fixed by setting `nFeature_min: 1`, that's
  a regression, not a fix. Investigate why first.
- **Never run the pipeline on data you cannot find.** If H5 files
  aren't where the user said they were, stop and ask — do not invent
  paths.

## How to re-run a single step

```bash
Rscript scripts/04_normalize_cluster.R --config config.yaml
Rscript scripts/04_normalize_cluster.R --config config.yaml --no-plots  # faster
```

The `--no-plots` flag skips PNG generation for tighter iteration loops.

## Customizing the config — common asks

- **Mouse data:** set `qc.metrics.mt_pattern: "^mt-"` (case matters).
- **QC thresholds** (`qc.thresholds`): any bound can be `null` for no filter on that side.
  Tune from step 01's violins and `threshold_preview_pct`.
- **10x Flex / probe-based:** the 3'/5' defaults may not transfer; set
  thresholds from step 01. `percent.mt` is a weaker signal here.
- **New tissue:** set `annotation.markers.genes` to
  markers organized by compartment. The top-level group names are just
  for readability; only the leaf `cell_type: [genes]` mapping matters.
- **Cross-technology integration (e.g. Flex + 3′ in one object):** harmony/rpca
  correct an additive shift and can't bridge a modality gap where the two batches
  measure partly different genes. Set `integration.method: scvi` (a count-based deep
  model) and `integration.batch_var` to the batch column (e.g. `chemistry`). One-time
  setup: `cd scripts/scvi && uv sync`. Tune via the `integration.scvi` block; judge by
  whether SHARED cell types mix while rare/confounded populations survive (not by a
  uniformly-blended UMAP). See `scripts/scvi/README.md`.
- **Azimuth references:** set `annotation.azimuth.reference` to one of
  the names from https://azimuth.hubmapconsortium.org/references/
  (`pbmcref`, `lungref`, `kidneyref`, `heartref`, `bonemarrowref`,
  `tonsilref`, `mousecortexref`, etc.). Set to `null` to skip. Set
  `annotation.azimuth.primary_level` to that reference's level name for the
  per-cluster label (e.g. `ann_level_3` for lungref, `celltype.l2` for
  pbmcref); a level the reference lacks is an error.
- **Config layout changes:** keys from the previous layout are rejected by
  `load_config()` with a message mapping each to its new location
  (`LEGACY_KEYS` in `utils.R`); absent keys take `CONFIG_DEFAULTS`, which
  must stay in sync with `config.template.yaml`.
- **Adding sample metadata columns:** any extra field under a sample
  (e.g. `condition: "treated"`, `donor: "D01"`, `timepoint: 24`) is
  automatically attached to every cell from that sample as a Seurat
  metadata column — `attach_sample_metadata()` in `utils.R` handles it.
  No code changes needed.

## Things not to do

- **Do not delete `outputs/` or `logs/`** without asking. Each `seurat_*.rds`
  is a checkpoint the user may want; `logs/` carries failure history.
- **Do not hand-edit RDS files or run ad-hoc R sessions** that mutate the
  pipeline checkpoints. If something needs to change, change the config
  and rerun the affected step — that's the whole point of the design.
- **Do not invent new scripts or steps.** If you think a sixth step is
  needed, propose it to the user first; this repo's value is being a
  small, fixed contract. (`scripts/scvi/` is the one sanctioned exception:
  an opt-in Python subsystem invoked *by* step 04, not a new step in the chain.)
- **Do not rewrite `utils.R` to take dependencies on tidyverse, dplyr,
  data.table, etc.** It is intentionally limited to `yaml`, `jsonlite`,
  `ggplot2`, `Seurat`, and the optional packages already documented.
  Adding deps breaks the "drop-in" promise. (The scVI helper `run_scvi_integration`
  is fine because it *shells out* to the Python env via `system2` — it links no new
  R package, and nothing Python loads unless `method: scvi`.)
- **Missing deps for a REQUESTED feature hard-stop, they don't skip.** If a config
  enables doublet detection / harmony / scVI / Azimuth / SingleR but the package or
  scVI env is absent, the step `emit_error`s rather than silently degrading. Don't
  "fix" that by disabling the feature to make the error go away — install the dep, or
  confirm with the user that the feature really should be off.
- **Do not commit `config.yaml`.** It's gitignored on purpose — only
  `config.template.yaml` is checked in. Users edit their own copy.
- **Do not commit anything in `outputs/` or `logs/`.** Same reason.
- **Do not run the pipeline without telling the user first.** It writes
  GB-scale RDS files and can take a long time on large datasets.

## Known issues (deferred)

Found in review; not fixed yet. Keep them in mind when interpreting results.

- **scVI cache ignores counts.** `run_scvi_integration()` keys its cache on
  parameters plus cell and gene names, not count values. Changing
  `ambient.method`, `round_to_int` or `skip_samples` with the same surviving
  cells can reuse a latent trained on the old counts. Workaround: change
  `integration.scvi.run_label` (or delete `outputs/scvi/<run_label>/`) after
  any upstream change. Fix: add a count checksum to the cache key.
- **Steps don't check that their input RDS is current.** Each step reads
  whatever upstream RDS is on disk; running step 04 after a failed step 03
  silently uses the previous `seurat_filtered.rds`. `inputs.rds` records the
  path only. Check the upstream summary's `status`/`timestamp` (or the
  report `_index`) before running a step on its own. Fix: record the input's
  mtime/md5 and warn when the upstream summary isn't `success`/`warning`.
- **`run.sh` needs python3 + PyYAML** to resolve `paths.output_dir` /
  `log_dir`; without them it silently falls back to `./outputs` and `./logs`
  while R writes to the configured dirs, so logs and the report land in the
  wrong place. Fix: fail loudly when the lookup fails.
- **Doublets are called after the QC thresholds** (step 03 filters, then runs
  scDblFinder on the survivors). The usual recommendation is to call doublets
  on lightly filtered cells, then apply QC. This is a methods choice to
  decide, not a crash; `doublets.pct_of_tested` is relative to post-threshold
  cells.

## Where to read more

- `README.md` — full user-facing docs, including the AI iteration story.
- `config.template.yaml` — every config field with inline comments and
  tuning notes for different assays.
- The `*_summary.json` files in `outputs/` after a run — the actual
  data you'll reason about during QC iteration.
