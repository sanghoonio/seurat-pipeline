#!/usr/bin/env Rscript
# report.R
# Index every run of every step: the runs defined in the config plus any run
# folders on disk, each with its status, and write outputs/pipeline_report.json
# (index + every run's summary).
#
# Status: success | warning | failed | missing (defined, never finished) |
#         stale (an ancestor was rerun since, or the run's block changed in the
#         config since it ran) | not_in_config (folder on disk, block removed).
#
# Usage: Rscript scripts/report.R --config <path>
#        Rscript scripts/report.R --config <path> --list-labels <step>   (used by the SLURM runner)

suppressPackageStartupMessages({ library(yaml); library(jsonlite) })
.utils_dir <- local({
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else "scripts"
})
source(file.path(.utils_dir, "utils.R"))

a <- commandArgs(trailingOnly = TRUE)
arg <- function(flag) { i <- match(flag, a); if (is.na(i) || i == length(a)) NULL else a[i + 1] }
config_path <- arg("--config")
if (is.null(config_path) || !file.exists(config_path)) emit_error("usage: report.R --config <path> [--list-labels <step>]")
raw <- tryCatch(yaml::read_yaml(config_path), error = function(e)
  emit_error("Config could not be parsed", list(path = config_path, error = conditionMessage(e))))
runs <- parse_run_blocks(raw)

list_step <- arg("--list-labels")
if (!is.null(list_step)) {
  if (!list_step %in% names(STEPS)) emit_error(sprintf("unknown step '%s'", list_step), list(steps = names(STEPS)))
  cat(names(runs[[list_step]]), sep = "\n")
  quit(status = 0, save = "no")
}

cfg <- fill_defaults(raw[intersect(names(raw), SHARED_KEYS)], CONFIG_DEFAULTS[SHARED_KEYS])
root <- output_dir(cfg)

index <- list(); summaries <- list()
for (step in names(STEPS)) {
  on_disk <- if (dir.exists(file.path(root, step))) list.dirs(file.path(root, step), full.names = FALSE, recursive = FALSE) else character(0)
  for (label in union(names(runs[[step]]), on_disk)) {
    key <- paste(step, label, sep = "/")
    s <- read_summary(file.path(root, step, label), step)
    defined <- label %in% names(runs[[step]])
    stale <- character(0)
    if (!is.null(s)) {
      if (defined && !identical(s$run$settings_hash, settings_hash(runs[[step]][[label]]$settings)))
        stale <- c(stale, "settings changed in the config since this run")
      if (defined && !identical(s$input$label, runs[[step]][[label]]$input))
        stale <- c(stale, sprintf("input is now '%s' in the config", runs[[step]][[label]]$input %||% "none"))
      for (anc in s$lineage %||% list()) {
        cur <- read_summary(file.path(root, anc$dir), anc$step)
        if (is.null(cur) || !identical(cur$run$run_id, anc$run_id))
          stale <- c(stale, sprintf("%s/%s was rerun or removed since", anc$step, anc$label))
      }
      summaries[[key]] <- s
    }
    status <- if (!defined) "not_in_config" else if (is.null(s)) "missing" else
      if (identical(s$status, "failed")) "failed" else if (length(stale)) "stale" else s$status
    index[[key]] <- list(step = step, label = label, status = status,
                         summary_status = s$status, stale_reasons = I(stale),
                         run_id = s$run$run_id, timestamp = s$timestamp,
                         input = if (!is.null(s$input)) paste(s$input$step, s$input$label, sep = "/") else NULL,
                         pipeline_commit = s$pipeline_commit,
                         n_warnings = length(s$warnings %||% list()))
  }
}

out <- file.path(root, "pipeline_report.json")
jsonlite::write_json(list(config = normalizePath(config_path), generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
                          `_index` = index, runs = summaries),
                     out, pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null", digits = NA)

w <- max(nchar(names(index)), 10)
cat(sprintf("%-*s  %-13s  %s\n", w, "run", "status", "notes"))
for (k in names(index)) {
  r <- index[[k]]
  notes <- c(if (r$n_warnings > 0) sprintf("%d warning(s)", r$n_warnings), unlist(r$stale_reasons))
  cat(sprintf("%-*s  %-13s  %s\n", w, k, r$status, paste(notes, collapse = "; ")))
}
cat(sprintf("\nWrote %s\n", out))
