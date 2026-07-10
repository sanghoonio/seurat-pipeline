#!/usr/bin/env Rscript
# 03_filter.R
# Apply QC thresholds from config. Optional doublet detection via scDblFinder.
# Writes: seurat_filtered.rds, filter_summary.json

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
cfg <- load_config(args$config_path)
set.seed(cfg$project$seed)

outdir <- output_dir(cfg)

# --- Load input ---
# Prefer decontaminated; fall back to raw
input_rds <- file.path(outdir, "seurat_decontaminated.rds")
if (!file.exists(input_rds)) {
  input_rds <- file.path(outdir, "seurat_raw.rds")
}
if (!file.exists(input_rds)) {
  emit_error("No input RDS found. Run steps 01/02 first.", list(checked = c("seurat_decontaminated.rds", "seurat_raw.rds")))
}

log_msg(sprintf("Loading %s", basename(input_rds)))
obj <- readRDS(input_rds)

total_before <- ncol(obj)
log_msg(sprintf("Cells before filtering: %d", total_before))

# --- Track per-sample counts ---
sample_ids <- unique(obj$sample_id)
before_counts <- table(obj$sample_id)

# --- Apply filters ---
qc <- cfg$qc

pass_nfeat_min <- obj$nFeature_RNA >= qc$nFeature_min
pass_nfeat_max <- obj$nFeature_RNA <= qc$nFeature_max
pass_ncount_min <- obj$nCount_RNA >= qc$nCount_min
pass_ncount_max <- obj$nCount_RNA <= qc$nCount_max
pass_mt <- obj$percent.mt <= qc$percent_mt_max

pass_all <- pass_nfeat_min & pass_nfeat_max & pass_ncount_min & pass_ncount_max & pass_mt

filter_counts <- list(
  nFeature_min = sum(!pass_nfeat_min),
  nFeature_max = sum(!pass_nfeat_max),
  nCount_min = sum(!pass_ncount_min),
  nCount_max = sum(!pass_ncount_max),
  percent_mt = sum(!pass_mt)
)

obj <- obj[, pass_all]
after_qc <- ncol(obj)
log_msg(sprintf("Cells after QC filtering: %d (removed %d)", after_qc, total_before - after_qc))

# --- Doublet detection ---
doublet_info <- list(enabled = FALSE)

if (isTRUE(qc$doublet_detection)) {
  # Requested but dep missing -> HARD STOP (don't silently ship un-de-doubleted data).
  if (!requireNamespace("scDblFinder", quietly = TRUE))
    emit_error("qc.doublet_detection is enabled but the 'scDblFinder' package is not installed",
               list(fix = "BiocManager::install('scDblFinder')  (or set qc.doublet_detection: false)"))
  {
    suppressPackageStartupMessages(library(scDblFinder))
    suppressPackageStartupMessages(library(SingleCellExperiment))

    log_msg("Running scDblFinder per sample")
    doublet_info <- list(enabled = TRUE, status = "success", per_sample = list())

    all_doublet_calls <- character(ncol(obj))
    names(all_doublet_calls) <- colnames(obj)
    all_doublet_scores <- numeric(ncol(obj))
    names(all_doublet_scores) <- colnames(obj)

    for (sid in sample_ids) {
      cells_idx <- which(obj$sample_id == sid)
      if (length(cells_idx) < 50) {
        log_msg(sprintf("  %s: too few cells (%d) for doublet detection, skipping", sid, length(cells_idx)))
        all_doublet_calls[cells_idx] <- "singlet"
        doublet_info$per_sample[[sid]] <- list(status = "skipped", reason = "too few cells")
        next
      }

      sub_obj <- obj[, cells_idx]
      sce <- as.SingleCellExperiment(sub_obj)

      sce <- tryCatch(
        scDblFinder(sce, dbr = qc$expected_doublet_rate %||% NULL),
        error = function(e) {
          log_msg(sprintf("  WARNING: scDblFinder failed for %s: %s", sid, e$message))
          NULL
        }
      )

      if (is.null(sce)) {
        all_doublet_calls[cells_idx] <- "singlet"
        doublet_info$per_sample[[sid]] <- list(status = "failed")
        next
      }

      calls <- sce$scDblFinder.class
      scores <- sce$scDblFinder.score
      all_doublet_calls[colnames(sub_obj)] <- as.character(calls)
      all_doublet_scores[colnames(sub_obj)] <- scores

      n_doublets <- sum(calls == "doublet")
      doublet_info$per_sample[[sid]] <- list(
        status = "success",
        cells_tested = length(cells_idx),
        doublets_found = n_doublets,
        doublet_pct = round(100 * n_doublets / length(cells_idx), 2)
      )

      log_msg(sprintf("  %s: %d doublets found (%.1f%%)", sid, n_doublets,
                       100 * n_doublets / length(cells_idx)))
    }

    obj$scDblFinder_class <- all_doublet_calls
    obj$scDblFinder_score <- all_doublet_scores

    before_doublet <- ncol(obj)
    obj <- subset(obj, scDblFinder_class == "singlet")
    doublets_removed <- before_doublet - ncol(obj)

    doublet_info$total_doublets_removed <- doublets_removed
    log_msg(sprintf("Removed %d doublets", doublets_removed))
  }
}

total_after <- ncol(obj)
after_counts <- table(obj$sample_id)

# --- Check for problems ---
pct_removed <- round(100 * (1 - total_after / total_before), 1)
if (pct_removed > 50) {
  log_msg(sprintf("WARNING: %.1f%% of cells filtered — review thresholds", pct_removed))
}

for (sid in sample_ids) {
  if (!(sid %in% names(after_counts)) || after_counts[sid] == 0) {
    emit_error("All cells filtered for a sample", list(
      sample = sid,
      thresholds = qc[c("nFeature_min", "nFeature_max", "nCount_min", "nCount_max", "percent_mt_max")]
    ))
  }
}

# --- Save ---
rds_path <- file.path(outdir, "seurat_filtered.rds")
log_msg(sprintf("Saving filtered object: %d cells", total_after))
saveRDS(obj, rds_path)

# --- JSON summary ---
per_sample_summary <- list()
for (sid in sample_ids) {
  b <- as.integer(before_counts[sid])
  a <- as.integer(after_counts[sid])
  per_sample_summary[[sid]] <- list(
    cells_before = b,
    cells_after = a,
    cells_removed = b - a,
    pct_removed = round(100 * (b - a) / b, 1)
  )
}

summary_data <- c(
  summary_header("03_filter", cfg),
  list(
    input = basename(input_rds),
    thresholds_applied = list(
      nFeature_min = qc$nFeature_min,
      nFeature_max = qc$nFeature_max,
      nCount_min = qc$nCount_min,
      nCount_max = qc$nCount_max,
      percent_mt_max = qc$percent_mt_max
    ),
    filter_counts = filter_counts,
    total = list(
      cells_before = total_before,
      cells_after_qc = after_qc,
      cells_after_doublet = total_after,
      cells_removed = total_before - total_after,
      pct_removed = pct_removed
    ),
    per_sample = per_sample_summary,
    doublet_detection = doublet_info,
    outputs = list(rds = rds_path),
    status = "success"
  )
)

write_json_summary(summary_data, file.path(outdir, "filter_summary.json"))

# --- Plots ---
if (args$make_plots) {
  log_msg("Generating filter plots")

  pre_csv <- file.path(outdir, "qc_metrics_decontaminated.csv")
  if (!file.exists(pre_csv)) pre_csv <- file.path(outdir, "qc_metrics.csv")

  if (file.exists(pre_csv)) {
    pre_df <- read.csv(pre_csv, stringsAsFactors = FALSE)
    pre_df$stage <- "Before filter"

    post_df <- data.frame(
      sample_id = obj$sample_id,
      nCount_RNA = obj$nCount_RNA,
      nFeature_RNA = obj$nFeature_RNA,
      percent.mt = obj$percent.mt,
      stringsAsFactors = FALSE
    )
    post_df$stage <- "After filter"

    comb <- rbind(
      pre_df[, c("sample_id", "nFeature_RNA", "stage")],
      post_df[, c("sample_id", "nFeature_RNA", "stage")]
    )
    comb$stage <- factor(comb$stage, levels = c("Before filter", "After filter"))

    p <- ggplot(comb, aes(x = sample_id, y = nFeature_RNA, fill = stage)) +
      geom_violin(alpha = 0.6, scale = "width") +
      geom_hline(yintercept = c(cfg$qc$nFeature_min, cfg$qc$nFeature_max),
                 linetype = "dashed", color = "red") +
      labs(title = "Genes per cell: before/after filtering", y = "nFeature_RNA", x = NULL, fill = NULL) +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    save_plot(p, "03_filter_nFeature_before_after.png", cfg, width = 8, height = 5)

    comb_nc <- rbind(
      pre_df[, c("sample_id", "nCount_RNA", "stage")],
      post_df[, c("sample_id", "nCount_RNA", "stage")]
    )
    comb_nc$stage <- factor(comb_nc$stage, levels = c("Before filter", "After filter"))

    p <- ggplot(comb_nc, aes(x = sample_id, y = nCount_RNA, fill = stage)) +
      geom_violin(alpha = 0.6, scale = "width") +
      scale_y_log10() +
      geom_hline(yintercept = c(cfg$qc$nCount_min, cfg$qc$nCount_max),
                 linetype = "dashed", color = "red") +
      labs(title = "UMI counts per cell: before/after filtering", y = "nCount_RNA (log10)", x = NULL, fill = NULL) +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    save_plot(p, "03_filter_nCount_before_after.png", cfg, width = 8, height = 5)

    comb_mt <- rbind(
      pre_df[, c("sample_id", "percent.mt", "stage")],
      post_df[, c("sample_id", "percent.mt", "stage")]
    )
    comb_mt$stage <- factor(comb_mt$stage, levels = c("Before filter", "After filter"))

    p <- ggplot(comb_mt, aes(x = sample_id, y = percent.mt, fill = stage)) +
      geom_violin(alpha = 0.6, scale = "width") +
      geom_hline(yintercept = cfg$qc$percent_mt_max, linetype = "dashed", color = "red") +
      labs(title = "Mitochondrial %: before/after filtering", y = "% MT", x = NULL, fill = NULL) +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    save_plot(p, "03_filter_pct_mt_before_after.png", cfg, width = 8, height = 5)
  }

  removal_df <- data.frame(
    sample = names(per_sample_summary),
    before = sapply(per_sample_summary, function(x) x$cells_before),
    after = sapply(per_sample_summary, function(x) x$cells_after),
    stringsAsFactors = FALSE
  )
  removal_df$removed <- removal_df$before - removal_df$after

  removal_long <- rbind(
    data.frame(sample = removal_df$sample, count = removal_df$after, category = "Kept", stringsAsFactors = FALSE),
    data.frame(sample = removal_df$sample, count = removal_df$removed, category = "Removed", stringsAsFactors = FALSE)
  )
  removal_long$category <- factor(removal_long$category, levels = c("Removed", "Kept"))

  p <- ggplot(removal_long, aes(x = sample, y = count, fill = category)) +
    geom_col() +
    scale_fill_manual(values = c("Kept" = "steelblue", "Removed" = "tomato")) +
    labs(title = "Cells kept vs removed per sample", y = "Cells", x = NULL, fill = NULL) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p, "03_filter_cells_removed.png", cfg, width = 8, height = 5)
}

log_msg("Step 03 complete")
