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

#' Load and validate config
load_config <- function(path) {
  cfg <- yaml::read_yaml(path)

  required <- c("project", "samples", "qc", "paths")
  missing <- setdiff(required, names(cfg))
  if (length(missing) > 0) {
    emit_error("Missing required config sections", list(missing = missing))
  }

  if (length(cfg$samples) == 0) {
    emit_error("No samples defined in config")
  }

  # Cell-calling config (optional). Default keeps current behavior.
  cc <- cfg$cell_calling
  if (!is.null(cc)) {
    m <- cc$method %||% "filtered"
    if (!m %in% c("filtered", "emptydrops")) {
      emit_error("cell_calling.method must be 'filtered' or 'emptydrops'", list(got = m))
    }
    u <- cc$use %||% "union"
    if (!u %in% c("union", "emptydrops", "filtered")) {
      emit_error("cell_calling.use must be 'union', 'emptydrops', or 'filtered'", list(got = u))
    }
  }
  emptydrops_on <- !is.null(cc) && (cc$method %||% "filtered") == "emptydrops"

  # Per-sample validation: must have id and a filtered source. emptyDrops
  # cell-calling additionally requires a raw source.
  for (i in seq_along(cfg$samples)) {
    s <- cfg$samples[[i]]
    if (is.null(s$id)) {
      emit_error("Sample missing 'id' field", list(index = i))
    }
    if (is.null(sample_filtered_path(s))) {
      emit_error("Sample missing 'filtered' (or 'filtered_h5') field", list(sample = s$id))
    }
    if (emptydrops_on && is.null(sample_raw_path(s))) {
      emit_error("cell_calling.method = emptydrops requires 'raw' (or 'raw_h5') for every sample",
                 list(sample = s$id))
    }
  }

  cfg
}

#' Per-sample matrix-source accessors. Accept clean names (`filtered`/`raw`)
#' or the legacy `filtered_h5`/`raw_h5` aliases. Each may be an .h5 file OR a
#' 10x mtx directory.
sample_filtered_path <- function(s) s$filtered %||% s$filtered_h5
sample_raw_path <- function(s) s$raw %||% s$raw_h5

#' Write JSON summary to file
write_json_summary <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(data, path, pretty = TRUE, auto_unbox = TRUE, null = "null")
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
save_plot <- function(p, filename, cfg, width = 10, height = 7, dpi = 200) {
  path <- file.path(plot_dir(cfg), filename)
  ggplot2::ggsave(path, plot = p, width = width, height = height, dpi = dpi)
  log_msg(paste("Saved plot:", path))
}

#' Build standard summary header
summary_header <- function(step_name, cfg) {
  list(
    step = step_name,
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    project = cfg$project$name,
    seed = cfg$project$seed
  )
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

#' Call cells from a raw count matrix with DropletUtils::emptyDrops.
#' Returns the significant-cell barcodes plus the full emptyDrops DataFrame
#' (so the caller can record per-barcode FDR for diagnostics). Set operations
#' against the Cell Ranger filtered set are done by the caller (step 00).
call_cells_emptydrops <- function(raw_counts, fdr = 0.01, lower = 100,
                                  niters = 10000, retain = NULL, alpha = NULL,
                                  test_ambient = FALSE, seed = 42) {
  if (!requireNamespace("DropletUtils", quietly = TRUE)) {
    emit_error("cell_calling.method = emptydrops requires DropletUtils",
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
#' `scvi_cfg$script` to point at an existing interpreter + script, else returns
#' {success=FALSE} and the caller falls back to un-integrated PCA. Adds NO R
#' package dependency — it only shells out.
#'
#' Caches on the call parameters: if a prior latent for this run_label exists and
#' the params are unchanged, the (expensive) training is skipped and the cached
#' latent reused. Any param change retrains.
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
  # silently degrade to un-integrated PCA. (Genuine runtime failures below still soft-fall.)
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
    seed                = as.integer(scvi_cfg$seed %||% cfg$project$seed %||% 42)
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
#' named by cell_calling.use (union / emptydrops / filtered). Returns NULL if
#' no cell-calling output exists for this sample (callers fall back to filtered).
read_called_barcodes <- function(sid, cfg) {
  f <- file.path(output_dir(cfg), "cell_calling", paste0(sid, "_barcodes.csv"))
  if (!file.exists(f)) return(NULL)
  use <- (cfg$cell_calling %||% list())$use %||% "union"
  df <- read.csv(f, stringsAsFactors = FALSE)
  sel <- switch(use,
    "union"      = as.logical(df$in_filtered) | as.logical(df$in_emptydrops),
    "emptydrops" = as.logical(df$in_emptydrops),
    "filtered"   = as.logical(df$in_filtered)
  )
  df$barcode[sel]
}
