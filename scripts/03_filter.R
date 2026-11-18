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
cfg <- load_config(args$config_path, step = "03_filter", summary_file = "filter_summary.json")
set.seed(cfg$project$seed)

outdir <- output_dir(cfg)

# --- Load input ---
# Prefer decontaminated; fall back to raw
input_rds <- file.path(outdir, "seurat_decontaminated.rds")
if (!file.exists(input_rds)) {
  input_rds <- file.path(outdir, "seurat_raw.rds")
  if (file.exists(input_rds))
    warn_msg("seurat_decontaminated.rds not found; filtering seurat_raw.rds (step 02 not run?)")
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
qc <- cfg$qc$thresholds
dbl <- cfg$qc$doublets

# A null bound means no filter on that side.
pass_nfeat_min <- obj$nFeature_RNA >= (qc$nFeature_min %||% -Inf)
pass_nfeat_max <- obj$nFeature_RNA <= (qc$nFeature_max %||% Inf)
pass_ncount_min <- obj$nCount_RNA >= (qc$nCount_min %||% -Inf)
pass_ncount_max <- obj$nCount_RNA <= (qc$nCount_max %||% Inf)
# percent.mt is NaN for 0-UMI cells; with a bound set they fail it.
pass_mt <- if (is.null(qc$percent_mt_max)) rep(TRUE, ncol(obj)) else
  !is.na(obj$percent.mt) & obj$percent.mt <= qc$percent_mt_max

pass_all <- pass_nfeat_min & pass_nfeat_max & pass_ncount_min & pass_ncount_max & pass_mt

# Cells failing each threshold (a cell can fail several; `any` is the union),
# overall and per sample. Keys match qc.thresholds.
fails <- data.frame(nFeature_min = !pass_nfeat_min, nFeature_max = !pass_nfeat_max,
                    nCount_min = !pass_ncount_min, nCount_max = !pass_ncount_max,
                    percent_mt_max = !pass_mt, any = !pass_all)
count_fails <- function(rows) as.list(vapply(fails[rows, , drop = FALSE], sum, integer(1)))
removed_by_threshold <- count_fails(rep(TRUE, nrow(fails)))
removed_by_threshold_per_sample <- lapply(setNames(nm = sample_ids), function(sid) count_fails(obj$sample_id == sid))
has_source <- "cell_source" %in% colnames(obj@meta.data)
source_before <- if (has_source) table(obj$sample_id, obj$cell_source) else NULL

obj <- obj[, pass_all]
after_qc <- ncol(obj)
log_msg(sprintf("Cells after QC filtering: %d (removed %d)", after_qc, total_before - after_qc))

# --- Doublet detection ---
doublet_info <- list(enabled = FALSE, status = "disabled")

if (isTRUE(dbl$remove)) {
  # Requested but dep missing -> HARD STOP (don't silently ship un-de-doubleted data).
  if (!requireNamespace("scDblFinder", quietly = TRUE))
    emit_error("qc.doublets.remove is true but the 'scDblFinder' package is not installed",
               list(fix = "BiocManager::install('scDblFinder')  (or set qc.doublets.remove: false)"))
  {
    suppressPackageStartupMessages(library(scDblFinder))
    suppressPackageStartupMessages(library(SingleCellExperiment))

    log_msg("Running scDblFinder per sample")
    doublet_info <- list(enabled = TRUE, expected_rate = dbl$expected_rate %||% "scDblFinder default (scaled by cell count)",
                         per_sample = list())

    all_doublet_calls <- character(ncol(obj))
    names(all_doublet_calls) <- colnames(obj)
    all_doublet_scores <- numeric(ncol(obj))
    names(all_doublet_scores) <- colnames(obj)

    for (sid in sample_ids) {
      cells_idx <- which(obj$sample_id == sid)
      if (length(cells_idx) < 50) {
        warn_msg(sprintf("%s: too few cells (%d) for doublet detection; all kept as singlets", sid, length(cells_idx)))
        all_doublet_calls[cells_idx] <- "singlet"
        doublet_info$per_sample[[sid]] <- list(status = "skipped", reason = "too few cells",
                                               n_tested = 0L, n_removed = 0L)
        next
      }

      sub_obj <- obj[, cells_idx]
      sce <- as.SingleCellExperiment(sub_obj)

      sce <- tryCatch(
        scDblFinder(sce, dbr = dbl$expected_rate),
        error = function(e) {
          warn_msg(sprintf("%s: scDblFinder failed (all kept as singlets): %s", sid, e$message))
          NULL
        }
      )

      if (is.null(sce)) {
        all_doublet_calls[cells_idx] <- "singlet"
        doublet_info$per_sample[[sid]] <- list(status = "failed", n_tested = 0L, n_removed = 0L)
        next
      }

      calls <- sce$scDblFinder.class
      scores <- sce$scDblFinder.score
      all_doublet_calls[colnames(sub_obj)] <- as.character(calls)
      all_doublet_scores[colnames(sub_obj)] <- scores

      n_doublets <- sum(calls == "doublet")
      doublet_info$per_sample[[sid]] <- list(
        status = "success",
        n_tested = length(cells_idx),         # cells left after the QC thresholds
        n_removed = n_doublets,
        pct_of_tested = round(100 * n_doublets / length(cells_idx), 2)
      )

      log_msg(sprintf("  %s: %d doublets found (%.1f%%)", sid, n_doublets,
                       100 * n_doublets / length(cells_idx)))
    }

    obj$scDblFinder_class <- all_doublet_calls
    obj$scDblFinder_score <- all_doublet_scores

    before_doublet <- ncol(obj)
    obj <- subset(obj, scDblFinder_class == "singlet")
    doublets_removed <- before_doublet - ncol(obj)

    doublet_info$n_removed <- doublets_removed
    st <- vapply(doublet_info$per_sample, function(x) x$status, "")
    doublet_info$status <- if (all(st == "success")) "success" else if (any(st == "success")) "partial" else "failed"
    log_msg(sprintf("Removed %d doublets", doublets_removed))
  }
}

total_after <- ncol(obj)
after_counts <- table(obj$sample_id)

# --- Check for problems ---
pct_removed <- round(100 * (1 - total_after / total_before), 1)
if (pct_removed > 50) {
  warn_msg(sprintf("%.1f%% of cells filtered; review thresholds", pct_removed))
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
source_after <- if (has_source) table(obj$sample_id, obj$cell_source) else NULL
by_source <- function(tab, sid) if (is.null(tab)) NULL else as.list(setNames(as.integer(tab[sid, ]), colnames(tab)))
per_sample_summary <- list()
for (sid in sample_ids) {
  b <- as.integer(before_counts[sid])
  a <- as.integer(after_counts[sid])
  rt <- removed_by_threshold_per_sample[[sid]]
  dbl_s <- doublet_info$per_sample[[sid]]
  per_sample_summary[[sid]] <- c(
    list(n_cells_before = b,
         removed_by_threshold = rt,
         n_after_thresholds = b - rt$any,
         doublets = dbl_s,
         n_cells_after = a,
         pct_removed = round(100 * (b - a) / b, 1)),
    if (has_source) list(n_cells_by_source = list(before = by_source(source_before, sid),
                                                  after = by_source(source_after, sid)))
  )
}

summary_data <- c(
  summary_header("03_filter", cfg, packages = if (isTRUE(dbl$remove)) "scDblFinder" else character(0)),
  list(
    inputs = list(rds = rel_out(input_rds, cfg)),
    params = list(thresholds = qc, doublets = dbl),
    total = c(
      list(n_cells_before = total_before,
           removed_by_threshold = removed_by_threshold,
           n_after_thresholds = after_qc,
           n_doublets_removed = doublet_info$n_removed %||% 0L,
           n_cells_after = total_after,
           pct_removed = pct_removed),
      if (has_source) list(n_cells_by_source = list(before = as.list(colSums(source_before)),
                                                    after = as.list(colSums(source_after))))
    ),
    per_sample = per_sample_summary,
    doublets = doublet_info[c("enabled", "status", intersect("expected_rate", names(doublet_info)),
                              intersect("n_removed", names(doublet_info)))],
    outputs = list(rds = rel_out(rds_path, cfg)),
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
      qc_threshold_lines("h", c(cfg$qc$thresholds$nFeature_min, cfg$qc$thresholds$nFeature_max), color = "red") +
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
      qc_threshold_lines("h", c(cfg$qc$thresholds$nCount_min, cfg$qc$thresholds$nCount_max), color = "red") +
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
      qc_threshold_lines("h", cfg$qc$thresholds$percent_mt_max, color = "red") +
      labs(title = "Mitochondrial %: before/after filtering", y = "% MT", x = NULL, fill = NULL) +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    save_plot(p, "03_filter_pct_mt_before_after.png", cfg, width = 8, height = 5)
  }

  # Cells kept per sample, with what was removed split into QC thresholds and
  # doublets; kept count and % removed labelled on each bar.
  removal_df <- data.frame(
    sample = names(per_sample_summary),
    before = sapply(per_sample_summary, function(x) x$n_cells_before),
    thresholds = sapply(per_sample_summary, function(x) x$removed_by_threshold$any),
    doublets = sapply(per_sample_summary, function(x) x$doublets$n_removed %||% 0L),
    after = sapply(per_sample_summary, function(x) x$n_cells_after),
    stringsAsFactors = FALSE
  )
  removal_long <- rbind(
    data.frame(sample = removal_df$sample, count = removal_df$after, category = "Kept"),
    data.frame(sample = removal_df$sample, count = removal_df$thresholds, category = "Removed: QC thresholds"),
    data.frame(sample = removal_df$sample, count = removal_df$doublets, category = "Removed: doublets")
  )
  removal_long$category <- factor(removal_long$category,
                                  levels = c("Removed: doublets", "Removed: QC thresholds", "Kept"))
  removal_df$label <- sprintf("%s kept\n(%.1f%% removed)", format(removal_df$after, big.mark = ","),
                              100 * (removal_df$before - removal_df$after) / removal_df$before)

  p <- ggplot(removal_long, aes(x = sample, y = count, fill = category)) +
    geom_col() +
    geom_text(data = removal_df, aes(x = sample, y = before, label = label), inherit.aes = FALSE,
              vjust = -0.25, size = 3, lineheight = 0.9) +
    scale_fill_manual(values = c("Kept" = "steelblue", "Removed: QC thresholds" = "tomato",
                                 "Removed: doublets" = "goldenrod")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    scale_x_discrete(expand = expansion(add = 0.9)) +   # room for the outer bars' labels
    labs(title = sprintf("Cells per sample after filtering (total %s kept of %s)",
                         format(sum(removal_df$after), big.mark = ","), format(sum(removal_df$before), big.mark = ",")),
         y = "Cells", x = NULL, fill = NULL) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p, "03_filter_cells_removed.png", cfg, width = 9, height = 5.5)
}

log_msg("Step 03 complete")
