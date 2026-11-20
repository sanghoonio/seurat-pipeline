#!/usr/bin/env Rscript
# 03_doublets.R
# Call doublets with scDblFinder, per sample, on every cell of the 02_ambient run
# named by `input` (before any QC threshold). Records calls and scores only; the
# 04_filter run that reads this run removes them (filter.remove_doublets).
# Writes (run folder): doublet_calls.csv, doublets_summary.json

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})
.utils_dir <- local({
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else "scripts"
})
source(file.path(.utils_dir, "utils.R"))

args <- parse_args()
cfg <- load_config(args$config_path, step = "03_doublets", label = args$label)
set.seed(cfg$project$seed)

outdir <- run_dir(cfg)
dbl <- cfg$doublets
obj <- read_ancestor_rds(cfg, "02_ambient", "seurat_decontaminated.rds")
in_rds <- file.path(ancestor_dir(cfg, "02_ambient"), "seurat_decontaminated.rds")
sample_ids <- unique(obj$sample_id)

calls <- data.frame(barcode = colnames(obj), sample_id = obj$sample_id,
                    class = "singlet", score = NA_real_, stringsAsFactors = FALSE)
per_sample <- list()

if (!isTRUE(dbl$enabled)) {
  log_msg("doublets.enabled is false: every cell recorded as singlet")
  for (sid in sample_ids) per_sample[[sid]] <- list(status = "disabled", n_cells = sum(obj$sample_id == sid))
} else {
  # Requested but dep missing -> HARD STOP (don't silently ship un-de-doubleted data).
  if (!requireNamespace("scDblFinder", quietly = TRUE))
    emit_error("doublets.enabled is true but the 'scDblFinder' package is not installed",
               list(fix = "BiocManager::install('scDblFinder')  (or set doublets.enabled: false)"))
  suppressPackageStartupMessages({ library(scDblFinder); library(SingleCellExperiment) })
  log_msg("Running scDblFinder per sample")
  for (sid in sample_ids) {
    idx <- which(obj$sample_id == sid)
    if (length(idx) < 50) {
      warn_msg(sprintf("%s: too few cells (%d) for doublet detection; all recorded as singlets", sid, length(idx)))
      per_sample[[sid]] <- list(status = "skipped", reason = "too few cells", n_cells = length(idx), n_doublets = 0L)
      next
    }
    sce <- tryCatch(scDblFinder(as.SingleCellExperiment(obj[, idx]), dbr = dbl$expected_rate),
                    error = function(e) { warn_msg(sprintf("%s: scDblFinder failed (all recorded as singlets): %s", sid, e$message)); NULL })
    if (is.null(sce)) {
      per_sample[[sid]] <- list(status = "failed", n_cells = length(idx), n_doublets = 0L)
      next
    }
    calls$class[idx] <- as.character(sce$scDblFinder.class)
    calls$score[idx] <- round(sce$scDblFinder.score, 4)
    n_d <- sum(sce$scDblFinder.class == "doublet")
    per_sample[[sid]] <- list(status = "success", n_cells = length(idx), n_doublets = n_d,
                              pct_doublets = round(100 * n_d / length(idx), 2))
    log_msg(sprintf("  %s: %d doublets (%.1f%%)", sid, n_d, 100 * n_d / length(idx)))
    rm(sce)
  }
}

csv_path <- file.path(outdir, "doublet_calls.csv")
write.csv(calls, csv_path, row.names = FALSE)

st <- vapply(per_sample, function(x) x$status, "")
status_all <- if (!isTRUE(dbl$enabled)) "disabled" else if (all(st == "success")) "success" else
  if (any(st == "success")) "partial" else "failed"
summary_data <- c(
  summary_header("03_doublets", cfg, packages = if (isTRUE(dbl$enabled)) "scDblFinder" else character(0)),
  list(
    inputs = list(rds = rel_out(in_rds, cfg)),
    params = dbl,
    doublet_detection = status_all,   # partial/failed: some samples kept all cells as singlets
    n_cells = nrow(calls),
    n_doublets = sum(calls$class == "doublet"),
    pct_doublets = round(100 * mean(calls$class == "doublet"), 2),
    per_sample = per_sample,
    outputs = list(csv = rel_out(csv_path, cfg)),
    status = "success"
  )
)
write_json_summary(summary_data, file.path(outdir, STEPS[["03_doublets"]]$summary))

if (args$make_plots && isTRUE(dbl$enabled) && any(!is.na(calls$score))) {
  p <- ggplot(calls[!is.na(calls$score), ], aes(x = score, fill = class)) +
    geom_histogram(bins = 50) +
    facet_wrap(~sample_id, scales = "free_y") +
    scale_fill_manual(values = c(singlet = "grey60", doublet = "goldenrod")) +
    labs(title = "scDblFinder score per sample", x = "Doublet score", y = "Cells", fill = NULL) +
    theme_bw(base_size = 11)
  save_plot(p, "03_doublet_scores.png", cfg, width = 10, height = 7)

  rates <- data.frame(sample = names(per_sample),
                      pct = vapply(per_sample, function(x) x$pct_doublets %||% 0, numeric(1)))
  p <- ggplot(rates, aes(x = sample, y = pct)) +
    geom_col(fill = "goldenrod") +
    geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.3, size = 3) +
    geom_hline(yintercept = 100 * (dbl$expected_rate %||% 0), linetype = "dashed") +
    labs(title = "Doublets called per sample (dashed: expected_rate)", x = NULL, y = "% of cells") +
    theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p, "03_doublet_rate.png", cfg, width = 8, height = 5)
}

log_msg("Step 03 complete")
