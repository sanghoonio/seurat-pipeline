suppressPackageStartupMessages({
  library(yaml)
  library(jsonlite)
})

# Null-coalescing operator (defined for portability across R/Seurat versions)
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Parse command-line arguments
#' Expects: --config path/to/config.yaml [--no-plots]
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  config_path <- "config.yaml"
  make_plots <- TRUE

  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--config" && i < length(args)) {
      config_path <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--no-plots") {
      make_plots <- FALSE
      i <- i + 1
    } else {
      i <- i + 1
    }
  }

  if (!file.exists(config_path)) {
    emit_error("Config file not found", list(path = config_path))
  }

  list(config_path = config_path, make_plots = make_plots)
}

# Defaults applied to any key the config leaves out. They match
# config.template.yaml, so an absent key behaves like the documented one. A key
# set to null in the config stays null (e.g. a QC bound of null = no filter).
CONFIG_DEFAULTS <- list(
  project = list(seed = 42),
  paths = list(base_dir = ".", output_dir = "outputs", log_dir = "logs"),
  cell_calling = list(
    cells = "cellranger",
    emptydrops = list(fdr = 0.01, lower = 100, niters = 10000, retain = NULL,
                      alpha = NULL, test_ambient = FALSE)
  ),
  qc = list(
    metrics = list(mt_pattern = "^MT-", ribo_pattern = "^RP[SL]"),
    thresholds = list(nFeature_min = 200, nFeature_max = 7500, nCount_min = 500,
                      nCount_max = 50000, percent_mt_max = 15),
    doublets = list(remove = TRUE, expected_rate = 0.05),
    prefilter_umap = FALSE
  ),
  ambient = list(
    enabled = TRUE, method = "soupx", round_to_int = TRUE, skip_samples = list(),
    soupx = list(contamination_range = c(0.01, 0.8)),
    decontx = list(use_background = TRUE, max_iter = 500, delta = c(10, 10),
                   estimate_delta = TRUE)
  ),
  normalization = list(
    method = "SCTransform", vars_to_regress = "percent.mt",
    sctransform = list(variable_features = 3000, vst_flavor = "v2"),
    lognormalize = list(variable_features = 2000, scale_factor = 10000)
  ),
  reduction = list(npcs = 30, dims_use = 30, umap = list(n_neighbors = 30, min_dist = 0.3)),
  clustering = list(algorithm = "louvain", resolutions = c(0.3, 0.5, 0.8, 1.0),
                    default_resolution = 0.5),
  integration = list(enabled = FALSE, method = "harmony", batch_var = "sample_id",
                     scvi = list()),
  annotation = list(
    azimuth = list(reference = NULL, primary_level = NULL, plot_levels = list()),
    singler = list(reference = NULL),
    markers = list(genes = NULL, override_azimuth = list(), min_score = 0)
  )
)

# Keys from the previous config layout -> where they live now. Any of these in a
# config is an error (silently ignoring them would quietly change behavior).
LEGACY_KEYS <- list(
  "soupx"                              = "ambient (method: soupx; contamination_range under ambient.soupx)",
  "cell_calling.method"                = "cell_calling.cells (cellranger | emptydrops | union)",
  "cell_calling.use"                   = "cell_calling.cells (cellranger | emptydrops | union)",
  "cell_calling.fdr"                   = "cell_calling.emptydrops.fdr",
  "cell_calling.lower"                 = "cell_calling.emptydrops.lower",
  "cell_calling.niters"                = "cell_calling.emptydrops.niters",
  "cell_calling.retain"                = "cell_calling.emptydrops.retain",
  "cell_calling.alpha"                 = "cell_calling.emptydrops.alpha",
  "cell_calling.test_ambient"          = "cell_calling.emptydrops.test_ambient",
  "qc.nFeature_min"                    = "qc.thresholds.nFeature_min",
  "qc.nFeature_max"                    = "qc.thresholds.nFeature_max",
  "qc.nCount_min"                      = "qc.thresholds.nCount_min",
  "qc.nCount_max"                      = "qc.thresholds.nCount_max",
  "qc.percent_mt_max"                  = "qc.thresholds.percent_mt_max",
  "qc.mt_pattern"                      = "qc.metrics.mt_pattern",
  "qc.ribo_pattern"                    = "qc.metrics.ribo_pattern",
  "qc.doublet_detection"               = "qc.doublets.remove",
  "qc.expected_doublet_rate"           = "qc.doublets.expected_rate",
  "normalization.sct_vars_to_regress"  = "normalization.vars_to_regress (used by both methods)",
  "normalization.sct_variable_features"= "normalization.sctransform.variable_features",
  "normalization.sct_vst_flavor"       = "normalization.sctransform.vst_flavor",
  "normalization.scale_factor"         = "normalization.lognormalize.scale_factor",
  "normalization.log_variable_features"= "normalization.lognormalize.variable_features",
  "reduction.umap_n_neighbors"         = "reduction.umap.n_neighbors",
  "reduction.umap_min_dist"            = "reduction.umap.min_dist",
  "annotation.azimuth_reference"       = "annotation.azimuth.reference",
  "annotation.azimuth_levels"          = "annotation.azimuth.primary_level (first entry) and annotation.azimuth.plot_levels",
  "annotation.singler_reference"       = "annotation.singler.reference",
  "annotation.marker_genes"            = "annotation.markers.genes",
  "annotation.marker_override"         = "annotation.markers.override_azimuth",
  "integration.reference_sample"       = "(removed: never used)"
)

has_key <- function(cfg, dotted) {
  parts <- strsplit(dotted, ".", fixed = TRUE)[[1]]
  node <- cfg
  for (p in parts) {
    if (!is.list(node) || !(p %in% names(node))) return(FALSE)
    node <- node[[p]]
  }
  TRUE
}

#' Recursively fill keys absent from `cfg` with `defaults`. Keys present in cfg
#' (including explicit nulls) win; samples are never merged.
fill_defaults <- function(cfg, defaults) {
  for (k in names(defaults)) {
    if (!(k %in% names(cfg))) {
      cfg[k] <- list(defaults[[k]])
    } else if (is.list(defaults[[k]]) && length(defaults[[k]]) && !is.null(names(defaults[[k]])) &&
               is.list(cfg[[k]])) {
      cfg[[k]] <- fill_defaults(cfg[[k]], defaults[[k]])
    }
  }
  cfg
}

#' Sections (keys whose default is a named list) set to null. Rejected rather
#' than filled with defaults: `ambient: null` or `thresholds: null` most likely
#' means "off", and silently running the defaults would do the opposite.
null_sections <- function(cfg, defaults, prefix = "") {
  out <- character(0)
  for (k in names(defaults)) {
    d <- defaults[[k]]
    if (!(is.list(d) && length(d) && !is.null(names(d))) || !(k %in% names(cfg))) next
    key <- paste0(prefix, k)
    if (is.null(cfg[[k]])) out <- c(out, key)
    else if (is.list(cfg[[k]])) out <- c(out, null_sections(cfg[[k]], d, paste0(key, ".")))
  }
  out
}

check_choice <- function(value, choices, key) {
  if (is.null(value) || length(value) != 1 || !(value %in% choices))
    emit_error(sprintf("%s must be one of: %s", key, paste(choices, collapse = ", ")),
               list(got = value %||% "null"))
  invisible(value)
}

#' Load and validate config
load_config <- function(path, step = NULL, summary_file = NULL) {
  cfg <- yaml::read_yaml(path)
  attr(cfg, "config_path") <- normalizePath(path)
  # Start the step as soon as paths are known, so a config error below still
  # replaces the previous summary with a "failed" one.
  if (!is.null(step) && !is.null(cfg$paths)) begin_step(step, summary_file, cfg)

  required <- c("project", "samples", "paths")
  missing <- setdiff(required, names(cfg))
  if (length(missing) > 0) {
    emit_error("Missing required config sections", list(missing = missing))
  }
  if (length(cfg$samples) == 0) {
    emit_error("No samples defined in config")
  }

  legacy <- Filter(function(k) has_key(cfg, k), names(LEGACY_KEYS))
  if (length(legacy)) {
    emit_error("Config uses keys from the old layout; move them as listed (see config.template.yaml)",
               setNames(lapply(legacy, function(k) LEGACY_KEYS[[k]]), legacy))
  }

  nulls <- null_sections(cfg, CONFIG_DEFAULTS)
  if (length(nulls))
    emit_error("Config sections set to null; delete the section to use the defaults, or set its keys (e.g. enabled: false)",
               list(sections = nulls))

  cfg <- fill_defaults(cfg, CONFIG_DEFAULTS)
  for (k in c("reduction.npcs", "reduction.dims_use", "clustering.default_resolution")) {
    parts <- strsplit(k, ".", fixed = TRUE)[[1]]
    v <- cfg[[parts[1]]][[parts[2]]]
    if (is.null(v) || length(v) != 1 || !is.numeric(v))
      emit_error(sprintf("%s must be a single number", k), list(got = v %||% "null"))
  }

  check_choice(cfg$cell_calling$cells, c("cellranger", "emptydrops", "union"), "cell_calling.cells")
  check_choice(cfg$ambient$method, c("soupx", "decontx"), "ambient.method")
  check_choice(cfg$normalization$method, c("SCTransform", "LogNormalize"), "normalization.method")
  check_choice(cfg$integration$method, c("harmony", "rpca", "cca", "scvi"), "integration.method")
  check_choice(cfg$clustering$algorithm, c("louvain", "leiden"), "clustering.algorithm")

  if (isTRUE(cfg$integration$enabled) && cfg$normalization$method == "LogNormalize")
    emit_error("integration.enabled requires normalization.method: SCTransform",
               list(fix = "use SCTransform, or set integration.enabled: false"))
  res <- as.numeric(unlist(cfg$clustering$resolutions))
  if (!length(res) || !(as.numeric(cfg$clustering$default_resolution) %in% res))
    emit_error("clustering.default_resolution must be one of clustering.resolutions",
               list(default_resolution = cfg$clustering$default_resolution, resolutions = res))
  if (cfg$reduction$dims_use > cfg$reduction$npcs)
    emit_error("reduction.dims_use cannot exceed reduction.npcs",
               list(dims_use = cfg$reduction$dims_use, npcs = cfg$reduction$npcs))

  # Per-sample validation: must have id and a filtered source. emptyDrops
  # cell-calling additionally requires a raw source.
  emptydrops_on <- cfg$cell_calling$cells != "cellranger"
  for (i in seq_along(cfg$samples)) {
    s <- cfg$samples[[i]]
    if (is.null(s$id)) {
      emit_error("Sample missing 'id' field", list(index = i))
    }
    if (is.null(sample_filtered_path(s))) {
      emit_error("Sample missing 'filtered' field", list(sample = s$id))
    }
    if (emptydrops_on && is.null(sample_raw_path(s))) {
      emit_error(sprintf("cell_calling.cells = %s requires 'raw' for every sample", cfg$cell_calling$cells),
                 list(sample = s$id))
    }
  }

  attr(cfg, "config_path") <- normalizePath(path)
  if (!is.null(step)) .step$cfg <- cfg
  cfg
}

#' Per-sample matrix-source accessors: `filtered` / `raw`, each an .h5 file or a
#' 10x mtx directory. `filtered_h5` / `raw_h5` are still accepted as aliases.
sample_filtered_path <- function(s) s$filtered %||% s$filtered_h5
sample_raw_path <- function(s) s$raw %||% s$raw_h5

#' Write JSON summary to file
write_json_summary <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  # Attach this step's warnings; "success" with warnings becomes "warning".
  if (!is.null(.step$name)) {
    data$warnings <- I(.step$warnings)
    if (identical(data$status, "success") && length(.step$warnings)) data$status <- "warning"
  }
  jsonlite::write_json(data, path, pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null")
  log_msg(paste("Wrote summary:", path))
}

#' Emit structured error to stderr and exit
emit_error <- function(msg, details = list()) {
  err <- list(
    error = msg,
    details = details,
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
  write(jsonlite::toJSON(err, pretty = TRUE, auto_unbox = TRUE), stderr())
  # Replace the step summary with a failure record (see begin_step()).
  if (!is.null(.step$summary_path)) {
    try(jsonlite::write_json(
      c(summary_header(.step$name, .step$cfg),
        list(status = "failed", error = msg, details = details, warnings = I(.step$warnings))),
      .step$summary_path, pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null"), silent = TRUE)
  }
  quit(status = 1, save = "no")
}

#' Log message to stderr with timestamp
log_msg <- function(msg) {
  ts <- format(Sys.time(), "%H:%M:%S")
  message(sprintf("[%s] %s", ts, msg))
}

#' Resolve path relative to config base_dir (absolute paths pass through)
resolve_path <- function(path, cfg) {
  if (startsWith(path, "/")) return(path)
  base <- cfg$paths$base_dir %||% "."
  file.path(base, path)
}

#' Get output directory from config, create if needed
output_dir <- function(cfg) {
  d <- resolve_path(cfg$paths$output_dir %||% "outputs", cfg)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

#' Get log directory from config, create if needed
log_dir <- function(cfg) {
  d <- resolve_path(cfg$paths$log_dir %||% "logs", cfg)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

#' Get plot directory, create if needed
plot_dir <- function(cfg) {
  d <- file.path(resolve_path(cfg$paths$output_dir %||% "outputs", cfg), "plots")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

#' Save a ggplot as PNG
#' Distribution of per-cell QC metrics for JSON summaries: min, 1/5/25/50/75/
#' 95/99th percentiles and max of each metric present in `df`.
QC_QUANTILE_PROBS <- c(min = 0, q01 = 0.01, q05 = 0.05, q25 = 0.25, q50 = 0.5,
                       q75 = 0.75, q95 = 0.95, q99 = 0.99, max = 1)
qc_quantiles <- function(df, metrics = c("nCount_RNA", "nFeature_RNA", "percent.mt",
                                         "percent.ribo", "complexity")) {
  metrics <- intersect(metrics, colnames(df))
  setNames(lapply(metrics, function(m) {
    q <- stats::quantile(df[[m]], QC_QUANTILE_PROBS, names = FALSE, na.rm = TRUE)
    as.list(setNames(round(q, 3), names(QC_QUANTILE_PROBS)))
  }), metrics)
}

#' Threshold preview for step summaries:
#' % of cells each step-03 threshold in cfg$qc$thresholds WOULD remove (nothing is removed
#' here). Cells can fail several thresholds, so `any` is not the column sum.
#' Doublet removal (step 03) is not included.
threshold_preview <- function(df, qc) {
  fails <- list(
    nFeature_min   = df$nFeature_RNA < (qc$nFeature_min %||% -Inf),
    nFeature_max   = df$nFeature_RNA > (qc$nFeature_max %||% Inf),
    nCount_min     = df$nCount_RNA   < (qc$nCount_min   %||% -Inf),
    nCount_max     = df$nCount_RNA   > (qc$nCount_max   %||% Inf),
    percent_mt_max = if (is.null(qc$percent_mt_max)) rep(FALSE, nrow(df))
                     else is.na(df$percent.mt) | df$percent.mt > qc$percent_mt_max   # NaN for 0-UMI cells
  )
  fails$any <- Reduce(`|`, fails)
  lapply(fails, function(x) round(100 * mean(x, na.rm = TRUE), 2))
}

#' QC thresholds are optional: a null bound in cfg$qc$thresholds means "no filter on
#' that side". qc_pass() applies them with that rule.
qc_pass <- function(nFeature, nCount, percent_mt, qc) {
  nFeature >= (qc$nFeature_min %||% -Inf) & nFeature <= (qc$nFeature_max %||% Inf) &
    nCount >= (qc$nCount_min %||% -Inf) & nCount <= (qc$nCount_max %||% Inf) &
    (if (is.null(qc$percent_mt_max)) TRUE else !is.na(percent_mt) & percent_mt <= qc$percent_mt_max)
}

#' Dashed threshold lines for QC plots; draws nothing when every bound is null.
qc_threshold_lines <- function(axis, values, ...) {
  v <- unlist(values)
  if (!length(v)) return(NULL)
  if (axis == "h") geom_hline(yintercept = v, linetype = "dashed", ...)
  else geom_vline(xintercept = v, linetype = "dashed", ...)
}

save_plot <- function(p, filename, cfg, width = 10, height = 7, dpi = 200) {
  path <- file.path(plot_dir(cfg), filename)
  # bg = "white": themes without a plot background (theme_minimal, some Seurat
  # themes) otherwise save as transparent PNGs.
  ggplot2::ggsave(path, plot = p, width = width, height = height, dpi = dpi, bg = "white")
  log_msg(paste("Saved plot:", path))
}

#' Build standard summary header
summary_header <- function(step_name, cfg, packages = character(0)) {
  cfg_path <- attr(cfg, "config_path")
  list(
    step = step_name,
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    project = cfg$project$name,
    seed = cfg$project$seed,
    pipeline_commit = pipeline_commit(),
    config = list(path = cfg_path %||% NA_character_,
                  md5 = if (!is.null(cfg_path) && file.exists(cfg_path)) unname(tools::md5sum(cfg_path)) else NA_character_),
    versions = pkg_versions(unique(c("Seurat", "SeuratObject", packages)))
  )
}

#' Short git commit of the pipeline checkout ("-dirty" if it has uncommitted
#' changes); NA when not a git checkout.
pipeline_commit <- function() {
  d <- if (exists(".utils_dir")) get(".utils_dir") else "scripts"
  run <- function(...) tryCatch(suppressWarnings(system2("git", c("-C", shQuote(d), ...), stdout = TRUE, stderr = FALSE)),
                                error = function(e) character(0))
  sha <- run("rev-parse", "--short", "HEAD")
  if (!length(sha) || !nzchar(sha[1])) return(NA_character_)
  dirty <- length(run("status", "--porcelain", "--untracked-files=no")) > 0
  paste0(sha[1], if (dirty) "-dirty" else "")
}

pkg_versions <- function(pkgs) {
  v <- lapply(pkgs, function(p) if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else NA_character_)
  c(list(R = paste(R.version$major, R.version$minor, sep = ".")), setNames(v, pkgs))
}

#' Path relative to output_dir when it lives under it (so summaries stay valid
#' if the outputs folder is moved); otherwise unchanged.
rel_out <- function(p, cfg) {
  if (is.null(p) || is.na(p)) return(NA_character_)
  # Normalize through the deepest existing ancestor, so not-yet-written files
  # resolve symlinks (e.g. /var -> /private/var) the same way as the folder.
  norm <- function(x) {
    x <- path.expand(x)
    if (file.exists(x) || dirname(x) == x) return(normalizePath(x, mustWork = FALSE))
    file.path(norm(dirname(x)), basename(x))
  }
  od <- norm(output_dir(cfg))
  pp <- norm(p)
  if (startsWith(pp, paste0(od, "/"))) substring(pp, nchar(od) + 2) else p
}

#' Per-sample matrix inputs, for the `inputs` block of step summaries.
sample_inputs <- function(cfg) {
  setNames(lapply(cfg$samples, function(s) list(
    filtered = resolve_path(sample_filtered_path(s), cfg),
    raw = if (!is.null(sample_raw_path(s))) resolve_path(sample_raw_path(s), cfg) else NA_character_
  )), vapply(cfg$samples, function(s) s$id, ""))
}

# --- Per-step state: summary path and collected warnings ---
.step <- new.env()
.step$name <- NULL; .step$summary_path <- NULL; .step$cfg <- NULL; .step$warnings <- character(0)

#' Call at the start of every step. Deletes the step's previous summary so a
#' failed run can't leave a stale "success" behind; emit_error() then writes a
#' status: "failed" summary in its place.
begin_step <- function(step_name, summary_file, cfg) {
  path <- file.path(output_dir(cfg), summary_file)
  if (file.exists(path)) file.remove(path)
  .step$name <- step_name; .step$summary_path <- path; .step$cfg <- cfg
  .step$warnings <- character(0)
  # Any uncaught R error (not just emit_error) also writes a failed summary.
  if (!interactive()) options(error = function() emit_error("Unhandled R error", list(message = trimws(geterrmessage()))))
  log_msg(sprintf("Step %s, pipeline %s", step_name, pipeline_commit()))
}

#' Log a warning and record it in the step summary's `warnings` list.
warn_msg <- function(msg) {
  log_msg(paste("WARNING:", msg))
  .step$warnings <- c(.step$warnings, msg)
  invisible(msg)
}

#' Attach optional sample metadata fields to a Seurat object.
#' Only sets columns that are present in the sample entry, so the
#' template config can omit unused fields.
attach_sample_metadata <- function(obj, s) {
  obj$sample_id <- s$id
  reserved <- c("id", "filtered", "raw", "filtered_h5", "raw_h5")
  for (field in setdiff(names(s), reserved)) {
    val <- s[[field]]
    if (!is.null(val)) obj[[field]] <- val
  }
  obj
}

#' Read a 10x matrix from either an .h5 file or a 10x mtx directory.
#' Dispatches on the path: *.h5 -> Read10X_h5; directory -> Read10X.
#' A directory must contain the Cell Ranger v3 triplet
#' (barcodes.tsv.gz / features.tsv.gz / matrix.mtx.gz).
read_10x_any <- function(path, use.names = TRUE, unique.features = TRUE) {
  if (grepl("\\.h5$", path, ignore.case = TRUE)) {
    if (!file.exists(path)) {
      emit_error("Matrix .h5 file not found", list(path = path))
    }
    return(Read10X_h5(path, use.names = use.names, unique.features = unique.features))
  }
  if (dir.exists(path)) {
    needed <- c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
    have <- file.exists(file.path(path, needed))
    if (!all(have)) {
      emit_error("mtx directory missing required files",
                 list(path = path, missing = needed[!have]))
    }
    # Read10X returns column-2 symbols and de-dupes (== use.names/unique.features)
    return(Read10X(path))
  }
  emit_error("Matrix path is neither an .h5 file nor a directory", list(path = path))
}

#' md5 of a character vector (base R only: written to a temp file, tools::md5sum).
md5_of_strings <- function(x) {
  f <- tempfile(); on.exit(unlink(f))
  writeLines(as.character(x), f)
  unname(tools::md5sum(f))
}

#' Restrict a RAW matrix to the feature set of the Cell Ranger FILTERED matrix.
#' Cell Ranger's raw matrix carries every feature; the filtered one drops the
#' probes 10x excludes (for Flex: ~21k unvalidated / DEPRECATED_* probes, most of
#' which duplicate a valid probe's gene and so split its counts across two rows
#' and act as chemistry-only features). Raw-derived objects therefore keep only
#' the filtered feature set. Barcodes (the emptyDrops recovery axis) are untouched.
restrict_to_filtered_features <- function(raw_counts, filt_features, sid = "") {
  keep <- rownames(raw_counts) %in% filt_features
  n_drop <- sum(!keep)
  log_msg(sprintf("  [%s] features: raw=%d filtered=%d -> keeping %d (dropped %d raw-only)",
                  sid, nrow(raw_counts), length(filt_features), sum(keep), n_drop))
  raw_counts[keep, , drop = FALSE]
}

#' Call cells from a raw count matrix with DropletUtils::emptyDrops.
#' Returns the significant-cell barcodes plus the full emptyDrops DataFrame
#' (so the caller can record per-barcode FDR for diagnostics). Set operations
#' against the Cell Ranger filtered set are done by the caller (step 00).
call_cells_emptydrops <- function(raw_counts, fdr = 0.01, lower = 100,
                                  niters = 10000, retain = NULL, alpha = NULL,
                                  test_ambient = FALSE, seed = 42) {
  if (!requireNamespace("DropletUtils", quietly = TRUE)) {
    emit_error("cell_calling.cells = emptydrops/union requires DropletUtils",
               list(install = "BiocManager::install('DropletUtils')"))
  }
  set.seed(seed)
  ed_args <- list(m = raw_counts, lower = lower, niters = niters,
                  test.ambient = test_ambient)
  if (!is.null(retain)) ed_args$retain <- retain
  if (!is.null(alpha)) ed_args$alpha <- alpha
  ed <- do.call(DropletUtils::emptyDrops, ed_args)

  # Barcodes with total <= lower are never tested (FDR = NA) -> not cells.
  is_cell <- !is.na(ed$FDR) & ed$FDR <= fdr
  list(emptydrops_barcodes = rownames(ed)[is_cell], ed = ed)
}

#' Attach cell-calling provenance to a merged object as metadata, if step 00 ran:
#'   cell_source     : "both" | "emptydrops_only" | "filtered_only" per cell
#'   emptydrops_FDR  : emptyDrops FDR per cell (NA for barcodes <= lower)
#' Barcodes in the cell_calling CSVs are un-prefixed; the object's are prefixed
#' with "<sample_id>_". No-op (returns obj unchanged) when no cell_calling output.
attach_cell_calling_meta <- function(obj, cfg) {
  if ((cfg$cell_calling$cells %||% "cellranger") == "cellranger") return(obj)
  cc_dir <- file.path(output_dir(cfg), "cell_calling")
  if (!dir.exists(cc_dir)) return(obj)
  src <- setNames(rep(NA_character_, ncol(obj)), colnames(obj))
  fdr <- setNames(rep(NA_real_, ncol(obj)), colnames(obj))
  for (s in cfg$samples) {
    bcf <- file.path(cc_dir, paste0(s$id, "_barcodes.csv"))
    if (!file.exists(bcf)) next
    b <- read.csv(bcf, stringsAsFactors = FALSE)
    pref <- paste0(s$id, "_", b$barcode)
    sv <- ifelse(b$in_filtered & b$in_emptydrops, "both",
          ifelse(b$in_emptydrops, "emptydrops_only", "filtered_only"))
    names(sv) <- pref
    common <- intersect(pref, colnames(obj)); src[common] <- sv[common]
    rkf <- file.path(cc_dir, paste0(s$id, "_barcode_ranks.csv"))
    if (file.exists(rkf)) {
      r <- read.csv(rkf, stringsAsFactors = FALSE)
      fv <- setNames(r$FDR, paste0(s$id, "_", r$barcode))
      common2 <- intersect(names(fv), colnames(obj)); fdr[common2] <- fv[common2]
    }
  }
  if (any(!is.na(src))) {
    obj$cell_source <- src
    obj$emptydrops_FDR <- fdr
  }
  obj
}

#' Cross-technology integration with an external scVI/sysVI model (opt-in).
#'
#' Bridges to a standalone Python script (scvi-tools) via a plain-file contract:
#' writes the RAW RNA counts (Matrix Market) + gene/barcode/metadata sidecars,
#' shells out to the script, and reads back a barcode-keyed latent CSV which it
#' aligns to the object's cells. Returns the embedding for the caller to install
#' as a DimReduc; the model itself is trained on counts, so the caller still runs
#' SCT/PCA separately for the variance report and downstream markers.
#'
#' Opt-in and guarded like harmony/DropletUtils: needs `scvi_cfg$python` and
#' `scvi_cfg$script` to point at an existing interpreter + script, else it
#' stops; a runtime failure returns {success=FALSE} and the caller stops too
#' (requested integration never falls back to un-integrated PCA). Adds NO R
#' package dependency — it only shells out.
#'
#' Caches on the call parameters plus an input fingerprint (cell and gene name
#' md5s): if a prior latent for this run_label exists and both are unchanged, the
#' (expensive) training is skipped and the cached latent reused. Any param change
#' or any change to the input cells/genes retrains.
run_scvi_integration <- function(obj, batch_var, scvi_cfg, cfg,
                                 pipeline_dir = getOption("seurat_pipeline_dir")) {
  run_label <- scvi_cfg$run_label %||% "default"
  scvi_dir <- file.path(output_dir(cfg), "scvi", run_label)
  dir.create(scvi_dir, recursive = TRUE, showWarnings = FALSE)

  # Defaults resolve to the pipeline-bundled env under scripts/scvi/ (pipeline_dir is the
  # scripts dir). Config integration.scvi.python/.script override. This keeps downstream
  # configs path-free: they need only `method: scvi` + tunables.
  bundled <- if (!is.null(pipeline_dir)) file.path(pipeline_dir, "scvi") else NULL
  script <- scvi_cfg$script %||% (if (!is.null(bundled)) file.path(bundled, "scvi_integrate.py") else NULL)
  python <- scvi_cfg$python %||% (if (!is.null(bundled)) file.path(bundled, ".venv", "bin", "python") else NULL)
  # Missing dependency -> HARD STOP. An explicitly requested integration method must not
  # silently degrade to un-integrated PCA. (Runtime failures below return success=FALSE,
  # and step 04 hard-stops on that too.)
  if (is.null(python) || is.null(script) || !file.exists(python %||% "") || !file.exists(script %||% "")) {
    emit_error("Integration method 'scvi' dependencies not found",
               list(python = python %||% "<unresolved>",
                    script = script %||% "<unresolved>",
                    fix = "Install the bundled env: cd scripts/scvi && uv sync (or set integration.scvi.python/.script)"))
  }
  if (!(batch_var %in% colnames(obj@meta.data))) {
    return(list(success = FALSE, error = sprintf("batch_var '%s' not in object metadata", batch_var)))
  }

  include_genes <- as.character(scvi_cfg$include_genes %||% character(0))
  params <- list(
    model               = scvi_cfg$model %||% "scvi",
    batch_key           = batch_var,
    n_latent            = as.integer(scvi_cfg$n_latent %||% 30),
    n_layers            = as.integer(scvi_cfg$n_layers %||% 2),
    gene_likelihood     = scvi_cfg$gene_likelihood %||% "nb",
    dispersion          = scvi_cfg$dispersion %||% "gene",
    n_hvg               = as.integer(scvi_cfg$n_hvg %||% 2000),
    min_cells_per_batch = as.integer(scvi_cfg$min_cells_per_batch %||% 3),
    include_genes       = include_genes,
    max_epochs          = scvi_cfg$max_epochs,   # may stay NULL -> scvi heuristic
    seed                = as.integer(scvi_cfg$seed %||% cfg$project$seed %||% 42),
    # Input fingerprint: a cached latent is only valid for the exact cells and
    # gene space it was trained on, so the key includes both (md5 of the names).
    n_cells             = ncol(obj),
    n_genes             = nrow(obj[["RNA"]]),
    cells_md5           = md5_of_strings(colnames(obj)),
    genes_md5           = md5_of_strings(rownames(obj[["RNA"]]))
  )
  latent_csv  <- file.path(scvi_dir, "scvi_latent.csv")
  params_json <- file.path(scvi_dir, "scvi_call_params.json")
  run_summary <- file.path(scvi_dir, "scvi_run_summary.json")

  # Cache: reuse latent iff a prior run had identical params.
  cache_hit <- FALSE
  if (file.exists(latent_csv) && file.exists(params_json)) {
    prev <- tryCatch(jsonlite::fromJSON(params_json), error = function(e) NULL)
    key <- function(p) jsonlite::toJSON(p, auto_unbox = TRUE, null = "null")
    if (!is.null(prev) && identical(key(prev), key(params))) {
      cache_hit <- TRUE
      log_msg(sprintf("scVI: reusing cached latent for run_label '%s' (params unchanged)", run_label))
    }
  }

  if (!cache_hit) {
    counts <- tryCatch(GetAssayData(obj, assay = "RNA", layer = "counts"),
                       error = function(e) LayerData(obj, assay = "RNA", layer = "counts"))
    if (!inherits(counts, "dgCMatrix")) counts <- as(counts, "dgCMatrix")
    Matrix::writeMM(counts, file.path(scvi_dir, "counts.mtx"))
    writeLines(rownames(counts), file.path(scvi_dir, "genes.txt"))
    writeLines(colnames(counts), file.path(scvi_dir, "barcodes.txt"))
    meta <- data.frame(barcode = colnames(obj), check.names = FALSE)
    meta[[batch_var]] <- as.character(obj@meta.data[[batch_var]])
    if ("sample_id" %in% colnames(obj@meta.data)) meta$sample_id <- as.character(obj$sample_id)
    write.csv(meta, file.path(scvi_dir, "metadata.csv"), row.names = FALSE)

    args <- c(script,
      "--counts", file.path(scvi_dir, "counts.mtx"),
      "--genes", file.path(scvi_dir, "genes.txt"),
      "--barcodes", file.path(scvi_dir, "barcodes.txt"),
      "--metadata", file.path(scvi_dir, "metadata.csv"),
      "--out", latent_csv, "--summary", run_summary,
      "--model", params$model, "--batch-key", params$batch_key,
      "--n-latent", params$n_latent, "--n-layers", params$n_layers,
      "--gene-likelihood", params$gene_likelihood, "--dispersion", params$dispersion,
      "--n-hvg", params$n_hvg, "--min-cells-per-batch", params$min_cells_per_batch,
      "--seed", params$seed)
    if (length(include_genes) > 0) args <- c(args, "--include-genes", paste(include_genes, collapse = ","))
    if (!is.null(params$max_epochs)) args <- c(args, "--max-epochs", params$max_epochs)

    log_msg(sprintf("scVI: training %s (n_latent=%d, n_layers=%d) via %s",
                    params$model, params$n_latent, params$n_layers, basename(python)))
    status <- system2(python, args = as.character(args), stdout = "", stderr = "")
    if (!identical(as.integer(status), 0L)) {
      return(list(success = FALSE, error = sprintf("scvi_integrate.py exited with status %s", status)))
    }
    jsonlite::write_json(params, params_json, auto_unbox = TRUE, null = "null", pretty = TRUE)
  }

  if (!file.exists(latent_csv)) return(list(success = FALSE, error = "scVI latent CSV not produced"))
  lat <- read.csv(latent_csv, row.names = 1, check.names = FALSE)
  missing_cells <- setdiff(colnames(obj), rownames(lat))
  if (length(missing_cells) > 0)
    return(list(success = FALSE, error = sprintf("%d cells missing from scVI latent", length(missing_cells))))
  emb <- as.matrix(lat[colnames(obj), , drop = FALSE])   # align by barcode, never by position
  colnames(emb) <- paste0("SCVI_", seq_len(ncol(emb)))
  summ <- tryCatch(jsonlite::fromJSON(run_summary), error = function(e) NULL)
  list(success = TRUE, embedding = emb, summary = summ, cache_hit = cache_hit, run_label = run_label)
}

#' Read the called-cell barcode set written by step 00, choosing the subset
#' named by cell_calling.cells (union / emptydrops). Returns NULL when
#' cells = cellranger (callers then use the Cell Ranger filtered matrix).
read_called_barcodes <- function(sid, cfg) {
  # Gate on the config, not just on file existence: step 00 output left over
  # from an earlier emptyDrops run must not be used once cells = cellranger.
  cells <- cfg$cell_calling$cells %||% "cellranger"
  if (cells == "cellranger") return(NULL)
  f <- file.path(output_dir(cfg), "cell_calling", paste0(sid, "_barcodes.csv"))
  if (!file.exists(f))
    emit_error(sprintf("cell_calling.cells = %s but step 00 output is missing; run step 00 first", cells),
               list(sample = sid, expected = f))
  df <- read.csv(f, stringsAsFactors = FALSE)
  sel <- switch(cells,
    "union"      = as.logical(df$in_filtered) | as.logical(df$in_emptydrops),
    "emptydrops" = as.logical(df$in_emptydrops)
  )
  df$barcode[sel]
}
