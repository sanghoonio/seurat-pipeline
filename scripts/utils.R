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

  # Per-sample validation: must have id and filtered_h5
  for (i in seq_along(cfg$samples)) {
    s <- cfg$samples[[i]]
    if (is.null(s$id)) {
      emit_error("Sample missing 'id' field", list(index = i))
    }
    if (is.null(s$filtered_h5)) {
      emit_error("Sample missing 'filtered_h5' field", list(sample = s$id))
    }
  }

  cfg
}

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
  for (field in setdiff(names(s), c("id", "filtered_h5", "raw_h5"))) {
    val <- s[[field]]
    if (!is.null(val)) obj[[field]] <- val
  }
  obj
}
