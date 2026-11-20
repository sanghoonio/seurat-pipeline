suppressPackageStartupMessages({
  library(yaml)
  library(jsonlite)
})

# Null-coalescing operator (defined for portability across R/Seurat versions)
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Parse command-line arguments: --config <path> [--label <label>] [--no-plots].
#' Unknown arguments are errors (a mistyped flag must not silently run something else).
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(config_path = "config.yaml", label = NULL, make_plots = TRUE)
  i <- 1
  while (i <= length(args)) {
    a <- args[i]
    if (a == "--config" && i < length(args)) { out$config_path <- args[i + 1]; i <- i + 2 }
    else if (a == "--label" && i < length(args)) { out$label <- args[i + 1]; i <- i + 2 }
    else if (a == "--no-plots") { out$make_plots <- FALSE; i <- i + 1 }
    else emit_error("Unknown command-line argument", list(argument = a,
                    usage = "--config <path> [--label <label>] [--no-plots]"))
  }
  if (!file.exists(out$config_path)) emit_error("Config file not found", list(path = out$config_path))
  out
}

# --- Steps ---------------------------------------------------------------------
# Each step owns a fixed set of config keys; a run block may contain only those.
# `parent` is the step whose run a block's `input` names.
STEPS <- list(
  "00_cell_calling"      = list(owns = "cell_calling.emptydrops", parent = NULL,
                                summary = "cell_calling_summary.json"),
  "01_load_qc_metrics"   = list(owns = c("cell_calling.cells", "qc.metrics"), parent = "00_cell_calling",
                                summary = "qc_summary.json"),
  "02_ambient"           = list(owns = "ambient", parent = "01_load_qc_metrics",
                                summary = "ambient_summary.json"),
  "03_doublets"          = list(owns = "doublets", parent = "02_ambient",
                                summary = "doublets_summary.json"),
  "04_filter"            = list(owns = c("qc.thresholds", "qc.prefilter_umap", "filter"), parent = "03_doublets",
                                summary = "filter_summary.json"),
  "05_normalize"         = list(owns = c("normalization", "reduction.npcs"), parent = "04_filter",
                                summary = "normalize_summary.json"),
  "06_cluster"           = list(owns = c("integration", "reduction.dims_use", "reduction.umap", "clustering", "markers"),
                                parent = "05_normalize", summary = "cluster_summary.json"),
  "07_annotate"          = list(owns = "annotation", parent = "06_cluster",
                                summary = "annotation_summary.json")
)
SHARED_KEYS <- c("project", "paths", "samples", "sample_metadata")
LABEL_RE <- "^[A-Za-z0-9_-]+$"

# Built-in defaults for every key a run block leaves out (explicit nulls are kept:
# a QC bound of null = no filter). They match config.template.yaml.
CONFIG_DEFAULTS <- list(
  project = list(seed = 42),
  paths = list(base_dir = ".", output_dir = "outputs", log_dir = "logs"),
  sample_metadata = list(file = NULL, id_column = "sample_id", columns = NULL),
  cell_calling = list(
    cells = "cellranger",
    emptydrops = list(fdr = 0.01, lower = 100, niters = 10000, retain = NULL,
                      alpha = NULL, test_ambient = FALSE)
  ),
  qc = list(
    metrics = list(mt_pattern = "^MT-", ribo_pattern = "^RP[SL]"),
    thresholds = list(nFeature_min = 200, nFeature_max = 7500, nCount_min = 500,
                      nCount_max = 50000, percent_mt_max = 15),
    prefilter_umap = FALSE
  ),
  ambient = list(
    enabled = TRUE, method = "soupx", round_to_int = TRUE, skip_samples = list(),
    soupx = list(contamination_range = c(0.01, 0.8)),
    decontx = list(use_background = TRUE, max_iter = 500, delta = c(10, 10),
                   estimate_delta = TRUE)
  ),
  doublets = list(enabled = TRUE, expected_rate = 0.05),
  filter = list(remove_doublets = TRUE),
  normalization = list(
    method = "SCTransform", vars_to_regress = "percent.mt", split_by = NULL,
    sctransform = list(variable_features = 3000, vst_flavor = "v2"),
    lognormalize = list(variable_features = 2000, scale_factor = 10000)
  ),
  reduction = list(npcs = 30, dims_use = 30, umap = list(n_neighbors = 30, min_dist = 0.3)),
  clustering = list(algorithm = "louvain", resolutions = c(0.3, 0.5, 0.8, 1.0),
                    default_resolution = 0.5),
  integration = list(enabled = FALSE, method = "harmony", batch_var = "sample_id",
                     scvi = list(python = NULL, script = NULL, model = "scvi", n_latent = 30,
                                 n_layers = 2, gene_likelihood = "nb", dispersion = "gene",
                                 n_hvg = 2000, min_cells_per_batch = 3, include_genes = list(),
                                 max_epochs = NULL, seed = NULL)),
  markers = list(only_pos = TRUE, min_pct = 0.25, logfc_threshold = 0.25, test = "wilcox"),
  annotation = list(
    azimuth = list(reference = NULL, primary_level = NULL, plot_levels = list()),
    singler = list(reference = NULL, labels = "main"),
    celltypist = list(models = list(), majority_voting = FALSE, python = NULL, script = NULL),
    custom_reference = list(path = NULL, label_columns = list(), dims = 50),
    markers = list(genes = NULL, override_azimuth = list(), min_score = 0),
    removed_cells = FALSE,
    compare_by = list()
  )
)

is_named_list <- function(x) is.list(x) && length(x) > 0 && !is.null(names(x))

has_key <- function(cfg, dotted) {
  node <- cfg
  for (p in strsplit(dotted, ".", fixed = TRUE)[[1]]) {
    if (!is.list(node) || !(p %in% names(node))) return(FALSE)
    node <- node[[p]]
  }
  TRUE
}
get_key <- function(cfg, dotted) {
  node <- cfg
  for (p in strsplit(dotted, ".", fixed = TRUE)[[1]]) node <- node[[p]]
  node
}
set_key <- function(cfg, dotted, value) {
  parts <- strsplit(dotted, ".", fixed = TRUE)[[1]]
  if (length(parts) == 1) { cfg[parts] <- list(value); return(cfg) }
  sub <- cfg[[parts[1]]]
  if (!is.list(sub)) sub <- list()
  cfg[parts[1]] <- list(set_key(sub, paste(parts[-1], collapse = "."), value))
  cfg
}

#' Recursively fill keys absent from `cfg` with `defaults`. Keys present in cfg
#' (including explicit nulls) win.
fill_defaults <- function(cfg, defaults) {
  for (k in names(defaults)) {
    if (!(k %in% names(cfg))) {
      cfg[k] <- list(defaults[[k]])
    } else if (is_named_list(defaults[[k]]) && is.list(cfg[[k]])) {
      cfg[[k]] <- fill_defaults(cfg[[k]], defaults[[k]])
    }
  }
  cfg
}

#' Unknown keys anywhere under `node` (compared with the defaults' structure), and
#' sections set to null whose default is a section. Leaves (default NULL, a
#' scalar, or an unnamed list) accept any value.
check_keys <- function(node, defaults, where, prefix = "") {
  for (k in names(node)) {
    key <- paste0(prefix, k)
    if (!(k %in% names(defaults)))
      emit_error(sprintf("%s: unknown config key '%s'", where, key),
                 list(allowed_here = names(defaults)))
    d <- defaults[[k]]
    if (is_named_list(d)) {
      if (is.null(node[[k]]))
        emit_error(sprintf("%s: '%s' is set to null; remove it to use the defaults, or set its keys", where, key))
      if (is.list(node[[k]]) && length(node[[k]])) check_keys(node[[k]], d, where, paste0(key, "."))
    }
  }
}

#' Keys in a run block that its step does not own.
check_owned <- function(node, owns, where, step, prefix = "") {
  for (k in names(node)) {
    key <- paste0(prefix, k)
    if (any(key == owns | startsWith(key, paste0(owns, ".")))) next   # inside an owned subtree
    if (any(startsWith(owns, paste0(key, "."))) && is.list(node[[k]])) {
      check_owned(node[[k]], owns, where, step, paste0(key, "."))      # parent of an owned key
      next
    }
    emit_error(sprintf("%s: step %s does not use '%s'", where, step, key),
               list(step_owns = owns))
  }
}

check_choice <- function(value, choices, key) {
  if (is.null(value) || length(value) != 1 || !(value %in% choices))
    emit_error(sprintf("%s must be one of: %s", key, paste(choices, collapse = ", ")),
               list(got = value %||% "null"))
  invisible(value)
}

#' A block's own settings: its owned keys, filled with built-in defaults.
block_settings <- function(block, step) {
  out <- list()
  for (o in STEPS[[step]]$owns) {
    d <- get_key(CONFIG_DEFAULTS, o)
    v <- if (has_key(block, o)) get_key(block, o) else d
    if (is_named_list(d) && is.list(v)) v <- fill_defaults(v, d)
    out <- set_key(out, o, v)
  }
  out
}

settings_hash <- function(settings)
  md5_of_strings(jsonlite::toJSON(settings, auto_unbox = TRUE, null = "null", digits = NA))

#' Value checks on one block's settings (everything not set is a valid default).
validate_settings <- function(s, where) {
  chk <- function(key, choices) if (has_key(s, key)) check_choice(get_key(s, key), choices, paste0(where, ": ", key))
  chk("cell_calling.cells", c("cellranger", "emptydrops", "union"))
  chk("ambient.method", c("soupx", "decontx"))
  chk("normalization.method", c("SCTransform", "LogNormalize"))
  chk("integration.method", c("harmony", "rpca", "cca", "scvi"))
  chk("clustering.algorithm", c("louvain", "leiden"))
  chk("annotation.singler.labels", c("main", "fine"))
  if (has_key(s, "clustering")) {
    res <- as.numeric(unlist(s$clustering$resolutions))
    dr <- s$clustering$default_resolution
    if (is.null(dr) || length(dr) != 1 || !is.numeric(dr) || !(dr %in% res))
      emit_error(sprintf("%s: clustering.default_resolution must be one of clustering.resolutions", where),
                 list(default_resolution = dr %||% "null", resolutions = res))
  }
  for (k in c("reduction.npcs", "reduction.dims_use")) if (has_key(s, k)) {
    v <- get_key(s, k)
    if (is.null(v) || length(v) != 1 || !is.numeric(v))
      emit_error(sprintf("%s: %s must be a single number", where, k), list(got = v %||% "null"))
  }
  if (has_key(s, "annotation")) {
    if (!is.logical(s$annotation$removed_cells) || length(s$annotation$removed_cells) != 1)
      emit_error(sprintf("%s: annotation.removed_cells must be true or false", where))
    cr <- s$annotation$custom_reference
    if (!is.null(cr$path) && !length(unlist(cr$label_columns)))
      emit_error(sprintf("%s: annotation.custom_reference.path is set but label_columns is empty", where))
  }
}

#' Parse and validate every run block in the config. Returns
#' list(<step> = list(<label> = list(block, settings, input))).
parse_run_blocks <- function(raw) {
  runs <- list()
  for (key in names(raw)) {
    if (key %in% SHARED_KEYS) next
    m <- regmatches(key, regexec("^([0-9]{2}_[a-z_]+)(/(.*))?$", key))[[1]]
    if (!length(m)) {
      if (key %in% names(CONFIG_DEFAULTS))
        emit_error(sprintf("top-level '%s' is the old layout; step settings now live in run blocks", key),
                   list(example = "02_ambient/decontx:\n  ambient: {method: decontx}", see = "config.template.yaml"))
      emit_error(sprintf("unknown top-level config key '%s'", key),
                 list(allowed = c(SHARED_KEYS, "<step>", "<step>/<label>"), steps = names(STEPS)))
    }
    step <- m[2]
    if (!step %in% names(STEPS)) emit_error(sprintf("unknown step '%s' in key '%s'", step, key), list(steps = names(STEPS)))
    label <- if (nzchar(m[4])) m[4] else "default"
    if (nzchar(m[4]) && tolower(label) == "default")
      emit_error(sprintf("'%s': the label 'default' is reserved; use the bare key '%s:'", key, step))
    if (!grepl(LABEL_RE, label)) emit_error(sprintf("'%s': labels may contain only letters, digits, _ and -", key))
    if (tolower(label) %in% tolower(names(runs[[step]])))
      emit_error(sprintf("'%s': another %s run has the same label (labels are compared ignoring case)", key, step))
    block <- raw[[key]] %||% list()
    if (!is.list(block) || (length(block) && is.null(names(block))))
      emit_error(sprintf("'%s' must be a mapping of settings", key))
    input <- block$input
    block$input <- NULL
    if (!is.null(input) && !(is.character(input) && length(input) == 1))
      emit_error(sprintf("'%s': input must be a quoted label string", key), list(got = input))
    check_keys(block, CONFIG_DEFAULTS, key)
    check_owned(block, STEPS[[step]]$owns, key, step)
    settings <- block_settings(block, step)
    validate_settings(settings, key)
    runs[[step]][[label]] <- list(settings = settings, input = input, key = key)
  }
  # Inputs must name a defined run of the parent step.
  for (step in names(runs)) for (label in names(runs[[step]])) {
    r <- runs[[step]][[label]]
    parent <- STEPS[[step]]$parent
    needs_parent <- !is.null(parent) &&
      !(step == "01_load_qc_metrics" && r$settings$cell_calling$cells == "cellranger")
    if (!needs_parent) {
      if (!is.null(r$input))
        emit_error(sprintf("'%s': input is not used here", r$key),
                   list(reason = if (is.null(parent)) "step 00 has no input" else "cell_calling.cells is cellranger, so no step 00 run is read"))
      next
    }
    want <- r$input %||% "default"
    if (!want %in% names(runs[[parent]]))
      emit_error(sprintf("'%s': input '%s' is not a defined %s run", r$key, want, parent),
                 list(defined = names(runs[[parent]]) %||% character(0)))
    runs[[step]][[label]]$input <- want
  }
  runs
}

#' Read a run's summary (or NULL).
read_summary <- function(dir, step) {
  f <- file.path(dir, STEPS[[step]]$summary)
  if (!file.exists(f)) return(NULL)
  tryCatch(jsonlite::fromJSON(f, simplifyVector = FALSE), error = function(e) NULL)
}

run_dir_of <- function(root, step, label) file.path(root, step, label)

#' Load the config for one run of `step`. Validates the whole config, selects the
#' run (`label`, or the only run defined), checks its parent run finished, and
#' builds the effective settings: defaults, then each ancestor's recorded
#' settings in chain order, then this run's own. Clears and starts the run folder.
load_config <- function(path, step, label = NULL) {
  raw <- tryCatch(yaml::read_yaml(path), error = function(e)
    emit_error("Config could not be parsed", list(path = path, error = conditionMessage(e))))
  for (k in c("project", "samples", "paths"))
    if (is.null(raw[[k]])) emit_error(sprintf("Missing required config section '%s'", k))
  if (length(raw$samples) == 0) emit_error("No samples defined in config")
  check_keys(raw[intersect(names(raw), c("paths", "sample_metadata"))], CONFIG_DEFAULTS, "config")
  runs <- parse_run_blocks(raw)

  cfg <- fill_defaults(raw[intersect(names(raw), SHARED_KEYS)], CONFIG_DEFAULTS[SHARED_KEYS])
  attr(cfg, "config_path") <- normalizePath(path)
  for (i in seq_along(cfg$samples)) {
    s <- cfg$samples[[i]]
    if (is.null(s$id)) emit_error("Sample missing 'id' field", list(index = i))
    if (is.null(sample_filtered_path(s))) emit_error("Sample missing 'filtered' field", list(sample = s$id))
  }
  cfg <- load_sample_metadata(cfg)

  # --- Select the run ---
  defined <- names(runs[[step]])
  if (!length(defined)) emit_error(sprintf("No run block for step %s in the config", step),
                                   list(add = sprintf("%s: {}  (default settings) or %s/<label>: {...}", step, step)))
  if (is.null(label)) {
    if (length(defined) > 1)
      emit_error(sprintf("Step %s has several runs; choose one with --label", step), list(labels = defined))
    label <- defined
  }
  if (!label %in% defined) emit_error(sprintf("No %s run labelled '%s'", step, label), list(defined = defined))
  run <- runs[[step]][[label]]
  root <- output_dir(cfg)

  # --- Parent run and lineage ---
  lineage <- list()
  parent <- NULL
  if (!is.null(run$input)) {
    pstep <- STEPS[[step]]$parent
    pdir <- run_dir_of(root, pstep, run$input)
    ps <- read_summary(pdir, pstep)
    if (is.null(ps))
      emit_error(sprintf("Input run %s/%s has no summary: it has not run, or is still running", pstep, run$input),
                 list(expected = file.path(pdir, STEPS[[pstep]]$summary)))
    if (!(ps$status %in% c("success", "warning")))
      emit_error(sprintf("Input run %s/%s did not finish successfully (status: %s)", pstep, run$input, ps$status))
    parent <- list(step = pstep, label = run$input, run_id = ps$run$run_id,
                   dir = rel_out(pdir, cfg))
    lineage <- c(ps$lineage %||% list(),
                 list(c(parent, list(settings = ps$run$settings))))
  }
  # Ancestors with no recorded run (01 without a 00 parent) leave defaults in place.
  for (k in names(CONFIG_DEFAULTS)) if (!k %in% SHARED_KEYS) cfg[k] <- list(CONFIG_DEFAULTS[[k]])
  overlay <- function(cfg, st, settings) {
    for (o in STEPS[[st]]$owns) if (has_key(settings, o)) cfg <- set_key(cfg, o, get_key(settings, o))
    cfg
  }
  for (a in lineage) cfg <- overlay(cfg, a$step, a$settings)
  cfg <- overlay(cfg, step, run$settings)
  validate_settings(cfg, sprintf("%s/%s (effective)", step, label))
  if (isTRUE(cfg$integration$enabled) && cfg$normalization$method == "LogNormalize")
    emit_error("integration.enabled requires normalization.method: SCTransform in the 05_normalize run it uses",
               list(fix = "use an SCTransform normalize run, or set integration.enabled: false"))
  if (cfg$reduction$dims_use > cfg$reduction$npcs)
    emit_error("reduction.dims_use cannot exceed reduction.npcs of the 05_normalize run",
               list(dims_use = cfg$reduction$dims_use, npcs = cfg$reduction$npcs))

  h <- settings_hash(run$settings)
  cfg$run <- list(step = step, label = label, dir = run_dir_of(root, step, label),
                  run_id = paste0(format(Sys.time(), "%Y%m%dT%H%M%S"), "-", Sys.getpid(), "-", substr(h, 1, 8)),
                  settings = run$settings, settings_hash = h,
                  parent = parent, lineage = lineage)
  begin_step(cfg)
  # The input run finished, but an older ancestor may have been rerun since it was made.
  for (a in head(lineage, -1)) {
    cur <- read_summary(run_dir_of(root, a$step, a$label), a$step)
    if (is.null(cur) || !identical(cur$run$run_id, a$run_id))
      warn_msg(sprintf("input is stale: %s/%s was rerun or removed after the input run %s/%s was made",
                       a$step, a$label, parent$step, parent$label))
  }
  cfg
}

#' Folder of this run, or of an ancestor run on its lineage (NULL if the step is
#' not an ancestor, e.g. 00 when cells = cellranger).
run_dir <- function(cfg) cfg$run$dir
ancestor <- function(cfg, step) {
  for (a in cfg$run$lineage) if (identical(a$step, step)) return(a)
  NULL
}
ancestor_dir <- function(cfg, step, required = TRUE) {
  a <- ancestor(cfg, step)
  if (is.null(a)) {
    if (required) emit_error(sprintf("This run has no %s ancestor", step))
    return(NULL)
  }
  resolve_path(a$dir, list(paths = list(base_dir = output_dir(cfg))))
}

#' Read an RDS written by an ancestor run and check it belongs to the run recorded
#' in the lineage (not a leftover from an earlier run of that label).
read_ancestor_rds <- function(cfg, step, file) {
  path <- file.path(ancestor_dir(cfg, step), file)
  if (!file.exists(path)) emit_error(sprintf("%s output not found", step), list(expected = path))
  log_msg(sprintf("Loading %s", rel_out(path, cfg)))
  obj <- readRDS(path)
  stamp <- tryCatch(Misc(obj, slot = "pipeline_run"), error = function(e) NULL)
  want <- ancestor(cfg, step)$run_id
  if (is.null(stamp) || !identical(stamp$run_id, want))
    emit_error(sprintf("%s does not match its run's summary (run_id %s vs %s); rerun %s", basename(path),
                       stamp$run_id %||% "none", want, step), list(path = path))
  obj
}

#' Stamp an object with this run's identity, then save it atomically.
save_run_rds <- function(obj, cfg, file) {
  Misc(obj, slot = "pipeline_run") <- list(step = cfg$run$step, label = cfg$run$label, run_id = cfg$run$run_id)
  path <- file.path(run_dir(cfg), file)
  tmp <- paste0(path, ".tmp")
  saveRDS(obj, tmp)
  file.rename(tmp, path)
  log_msg(sprintf("Saved %s (%d cells)", rel_out(path, cfg), ncol(obj)))
  path
}

# Sample-entry fields that are not metadata, and column names the pipeline sets itself.
SAMPLE_RESERVED <- c("id", "filtered", "raw", "filtered_h5", "raw_h5")
METADATA_RESERVED <- c("sample_id", "orig.ident", "nCount_RNA", "nFeature_RNA")

#' Merge columns of the sample_metadata file into cfg$samples (in memory), so
#' attach_sample_metadata() adds them to every cell like inline sample fields.
#' Every configured sample needs exactly one row; extra rows are ignored. A
#' column that is also an inline field of a sample is an error, not an override.
#' Only the file path, md5 and column names are recorded (for summaries), never values.
load_sample_metadata <- function(cfg) {
  sm <- cfg$sample_metadata
  if (is.null(sm$file)) return(cfg)
  path <- resolve_path(sm$file, cfg)
  if (!file.exists(path)) emit_error("sample_metadata.file not found", list(path = path))
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA"))
  idc <- sm$id_column %||% "sample_id"
  if (!idc %in% names(df))
    emit_error("sample_metadata.id_column is not a column of the file", list(id_column = idc, columns = names(df)))
  ids <- vapply(cfg$samples, function(s) as.character(s$id), "")
  keys <- as.character(df[[idc]])
  if (length(miss <- setdiff(ids, keys)))
    emit_error("sample_metadata file has no row for some samples", list(samples = miss, file = path))
  if (length(dup <- intersect(ids, keys[duplicated(keys)])))
    emit_error("sample_metadata file has more than one row for some samples", list(samples = dup, file = path))
  cols <- if (is.null(sm$columns)) setdiff(names(df), idc) else as.character(unlist(sm$columns))
  if (length(absent <- setdiff(cols, names(df))))
    emit_error("sample_metadata.columns are not in the file", list(missing = absent, file = path))
  if (length(bad <- intersect(cols, c(SAMPLE_RESERVED, METADATA_RESERVED))))
    emit_error("sample_metadata columns use reserved names", list(columns = bad))
  clash <- unique(unlist(lapply(cfg$samples, function(s) intersect(cols, names(s)))))
  if (length(clash))
    emit_error("sample_metadata columns are also inline sample fields in the config; keep each in one place",
               list(columns = clash))
  rows <- match(ids, keys)
  for (i in seq_along(cfg$samples))
    for (cn in cols) cfg$samples[[i]][cn] <- list(df[[cn]][rows[i]])
  cfg$sample_metadata$resolved <- list(file = path, md5 = unname(tools::md5sum(path)), columns = I(cols))
  cfg
}

#' Marker panel (annotation.markers.genes) as cell type -> genes. Accepts
#' compartment -> cell type -> genes or a flat cell type -> genes map; a cell
#' type defined twice keeps its first definition.
flatten_markers <- function(genes_cfg, warn = TRUE) {
  out <- list()
  add <- function(ct, g) {
    if (ct %in% names(out)) {
      if (warn) warn_msg(sprintf("marker cell type '%s' is defined more than once; keeping the first", ct))
      return(invisible())
    }
    out[[ct]] <<- unique(as.character(unlist(g)))
  }
  for (key in names(genes_cfg)) {
    val <- genes_cfg[[key]]
    if (is.list(val)) for (ct in names(val)) add(ct, val[[ct]]) else add(key, val)
  }
  out
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
  tmp <- paste0(path, ".tmp")
  jsonlite::write_json(data, tmp, pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null", digits = NA)
  file.rename(tmp, path)
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
  d <- file.path(run_dir(cfg), "plots")
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
  hdr <- list(
    step = step_name,
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    project = cfg$project$name,
    seed = cfg$project$seed,
    pipeline_commit = pipeline_commit(),
    config = list(path = cfg_path %||% NA_character_,
                  md5 = if (!is.null(cfg_path) && file.exists(cfg_path)) unname(tools::md5sum(cfg_path)) else NA_character_),
    versions = pkg_versions(unique(c("Seurat", "SeuratObject", packages))),
    # This run: identity, the settings it used (block + defaults), the run it read,
    # and every ancestor run (with their run_ids and settings) back to the root.
    run = if (!is.null(cfg$run)) cfg$run[c("step", "label", "run_id", "settings_hash", "settings")],
    input = cfg$run$parent,
    lineage = cfg$run$lineage %||% list()
  )
  # Only when a sample_metadata file is used: its path, md5 and column names.
  if (!is.null(cfg$sample_metadata$resolved)) hdr$sample_metadata <- cfg$sample_metadata$resolved
  hdr
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
begin_step <- function(cfg) {
  # The run's folder is cleared, so nothing from an earlier run of this label survives.
  d <- run_dir(cfg)
  unlink(d, recursive = TRUE)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  step_name <- cfg$run$step
  path <- file.path(d, STEPS[[step_name]]$summary)
  .step$name <- step_name; .step$summary_path <- path; .step$cfg <- cfg
  .step$warnings <- character(0)
  # Any uncaught R error (not just emit_error) also writes a failed summary.
  if (!interactive()) options(error = function() emit_error("Unhandled R error", list(message = trimws(geterrmessage()))))
  log_msg(sprintf("Step %s, run %s (%s), input %s, pipeline %s", step_name, cfg$run$label, cfg$run$run_id,
                  if (is.null(cfg$run$parent)) "none" else paste(cfg$run$parent$step, cfg$run$parent$label, sep = "/"),
                  pipeline_commit()))
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
  for (field in setdiff(names(s), SAMPLE_RESERVED)) {
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
  cc_dir <- file.path(ancestor_dir(cfg, "00_cell_calling"), "cell_calling")
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
#' md5s): if a prior latent for this cache key exists and both are unchanged, the
#' (expensive) training is skipped and the cached latent reused. Any param change,
#' any change to the input cells/genes, or a rerun of the 05 run retrains.
run_scvi_integration <- function(obj, batch_var, scvi_cfg, cfg,
                                 pipeline_dir = getOption("seurat_pipeline_dir")) {

  # Defaults resolve to the pipeline-bundled env under scripts/scvi/ (pipeline_dir is the
  # scripts dir). Config integration.scvi.python/.script override. This keeps downstream
  # configs path-free: they need only `method: scvi` + tunables.
  bundled <- if (!is.null(pipeline_dir)) file.path(pipeline_dir, "scvi") else NULL
  script <- scvi_cfg$script %||% (if (!is.null(bundled)) file.path(bundled, "scvi_integrate.py") else NULL)
  python <- scvi_cfg$python %||% (if (!is.null(bundled)) file.path(bundled, ".venv", "bin", "python") else NULL)
  # Missing dependency -> HARD STOP. An explicitly requested integration method must not
  # silently degrade to un-integrated PCA. (Runtime failures below return success=FALSE,
  # and step 06 hard-stops on that too.)
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
    genes_md5           = md5_of_strings(rownames(obj[["RNA"]])),
    # The normalize run it came from: a rerun upstream gives a new run_id, so a
    # cached latent trained on older counts is never reused.
    input_run_id        = ancestor(cfg, "05_normalize")$run_id
  )
  # Cache lives with the 05_normalize run (06 variants differing only in UMAP or
  # clustering reuse one model; rerunning 05 clears it), keyed on the params.
  cache_key <- substr(md5_of_strings(jsonlite::toJSON(params, auto_unbox = TRUE, null = "null", digits = NA)), 1, 12)
  scvi_dir <- file.path(ancestor_dir(cfg, "05_normalize"), "scvi", cache_key)
  dir.create(scvi_dir, recursive = TRUE, showWarnings = FALSE)
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
      log_msg(sprintf("scVI: reusing cached latent %s (params unchanged)", cache_key))
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
  list(success = TRUE, embedding = emb, summary = summ, cache_hit = cache_hit, cache_key = cache_key)
}

#' Read the called-cell barcode set written by step 00, choosing the subset
#' named by cell_calling.cells (union / emptydrops). Returns NULL when
#' cells = cellranger (callers then use the Cell Ranger filtered matrix).
read_called_barcodes <- function(sid, cfg) {
  # Gate on the config, not just on file existence: step 00 output left over
  # from an earlier emptyDrops run must not be used once cells = cellranger.
  cells <- cfg$cell_calling$cells %||% "cellranger"
  if (cells == "cellranger") return(NULL)
  f <- file.path(ancestor_dir(cfg, "00_cell_calling"), "cell_calling", paste0(sid, "_barcodes.csv"))
  if (!file.exists(f))
    emit_error(sprintf("cell_calling.cells = %s but the step 00 run has no barcodes for this sample", cells),
               list(sample = sid, expected = f))
  df <- read.csv(f, stringsAsFactors = FALSE)
  sel <- switch(cells,
    "union"      = as.logical(df$in_filtered) | as.logical(df$in_emptydrops),
    "emptydrops" = as.logical(df$in_emptydrops)
  )
  df$barcode[sel]
}

#' Annotate cells with CellTypist models via the bundled Python env (the scVI env,
#' scripts/scvi/.venv, which also carries celltypist). Exports raw RNA counts;
#' the script normalizes them as CellTypist expects (10k per cell, log1p) and
#' writes one CSV per model into out_dir. Returns list(success, per_model =
#' list(<model> = data.frame aligned to the object's cells), summary) or
#' list(success = FALSE, error). A missing env or script is a hard stop.
run_celltypist <- function(obj, models, ct_cfg, out_dir, pipeline_dir = getOption("seurat_pipeline_dir")) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  python <- ct_cfg$python %||% (if (!is.null(pipeline_dir)) file.path(pipeline_dir, "scvi", ".venv", "bin", "python") else NULL)
  script <- ct_cfg$script %||% (if (!is.null(pipeline_dir)) file.path(pipeline_dir, "celltypist", "celltypist_annotate.py") else NULL)
  if (is.null(python) || is.null(script) || !file.exists(python) || !file.exists(script))
    emit_error("annotation.celltypist.models is set but the CellTypist env or script was not found",
               list(python = python %||% "<unresolved>", script = script %||% "<unresolved>",
                    fix = "cd scripts/scvi && uv sync (installs celltypist), or set annotation.celltypist.python/.script"))
  keys <- sub("\\.pkl$", "", models)
  model_csv <- setNames(file.path(out_dir, paste0("celltypist_", keys, ".csv")), keys)
  summary_json <- file.path(out_dir, "celltypist_run_summary.json")
  for (f in c(model_csv, summary_json)) if (file.exists(f)) file.remove(f)   # only this run's models
  tmp <- file.path(out_dir, "celltypist_input")
  dir.create(tmp, showWarnings = FALSE)
  counts <- LayerData(obj, assay = "RNA", layer = "counts")
  if (!inherits(counts, "dgCMatrix")) counts <- as(counts, "dgCMatrix")
  Matrix::writeMM(counts, file.path(tmp, "counts.mtx"))
  writeLines(rownames(counts), file.path(tmp, "genes.txt"))
  writeLines(colnames(counts), file.path(tmp, "barcodes.txt"))
  args <- c(script, "--counts", file.path(tmp, "counts.mtx"), "--genes", file.path(tmp, "genes.txt"),
            "--barcodes", file.path(tmp, "barcodes.txt"), "--models", paste(models, collapse = ","),
            "--out-dir", out_dir, "--summary", summary_json,
            if (isTRUE(ct_cfg$majority_voting)) "--majority-voting")
  log_msg(sprintf("CellTypist: %s via %s", paste(models, collapse = ", "), basename(python)))
  status <- system2(python, args = as.character(args), stdout = "", stderr = "")
  unlink(tmp, recursive = TRUE)   # large; re-exported each run
  if (!identical(as.integer(status), 0L))
    return(list(success = FALSE, error = sprintf("celltypist_annotate.py exited with status %s", status)))
  per_model <- list()
  for (k in keys) {
    if (!file.exists(model_csv[[k]])) return(list(success = FALSE, error = sprintf("no output for model %s", k)))
    d <- read.csv(model_csv[[k]], check.names = FALSE, stringsAsFactors = FALSE)
    if (length(setdiff(colnames(obj), d$barcode)))
      return(list(success = FALSE, error = sprintf("CellTypist output for %s is missing cells", k)))
    per_model[[k]] <- d[match(colnames(obj), d$barcode), , drop = FALSE]
  }
  list(success = TRUE, per_model = per_model, csv = model_csv,
       summary = tryCatch(jsonlite::fromJSON(summary_json, simplifyVector = FALSE), error = function(e) NULL))
}
