#!/usr/bin/env Rscript
# 04_filter.R
# Apply QC thresholds to the corrected counts of the 02_ambient run behind the
# 03_doublets run named by `input`, and remove the cells that run called doublets
# (filter.remove_doublets). Removed cells are kept, tagged with why they were removed.
# Writes (run folder): seurat_filtered.rds, seurat_filtered_out.rds, filtered_out_cells.csv,
#         filter_summary.json

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
cfg <- load_config(args$config_path, step = "04_filter", label = args$label)
set.seed(cfg$project$seed)

outdir <- run_dir(cfg)
qc <- cfg$qc$thresholds
remove_doublets <- isTRUE(cfg$filter$remove_doublets)

# --- Inputs: corrected object (02) and doublet calls (03) ---
amb_dir <- ancestor_dir(cfg, "02_ambient")
dbl_dir <- ancestor_dir(cfg, "03_doublets")
in_rds <- file.path(amb_dir, "seurat_decontaminated.rds")
obj <- read_ancestor_rds(cfg, "02_ambient", "seurat_decontaminated.rds")
calls <- read.csv(file.path(dbl_dir, "doublet_calls.csv"), stringsAsFactors = FALSE)
if (length(miss <- setdiff(colnames(obj), calls$barcode)))
  emit_error("doublet_calls.csv does not cover every cell of the ambient run", list(n_missing = length(miss)))
calls <- calls[match(colnames(obj), calls$barcode), ]
obj$scDblFinder_class <- calls$class
obj$scDblFinder_score <- calls$score
dbl_detection <- read_summary(dbl_dir, "03_doublets")$doublet_detection %||% "unknown"

total_before <- ncol(obj)
sample_ids <- unique(obj$sample_id)
before_counts <- table(obj$sample_id)
log_msg(sprintf("Cells before filtering: %d", total_before))

# --- Thresholds (a null bound means no filter on that side) ---
pass_nfeat_min <- obj$nFeature_RNA >= (qc$nFeature_min %||% -Inf)
pass_nfeat_max <- obj$nFeature_RNA <= (qc$nFeature_max %||% Inf)
pass_ncount_min <- obj$nCount_RNA >= (qc$nCount_min %||% -Inf)
pass_ncount_max <- obj$nCount_RNA <= (qc$nCount_max %||% Inf)
# percent.mt is NaN for 0-UMI cells; with a bound set they fail it.
pass_mt <- if (is.null(qc$percent_mt_max)) rep(TRUE, ncol(obj)) else
  !is.na(obj$percent.mt) & obj$percent.mt <= qc$percent_mt_max
pass_thr <- pass_nfeat_min & pass_nfeat_max & pass_ncount_min & pass_ncount_max & pass_mt
is_dbl <- remove_doublets & obj$scDblFinder_class == "doublet"
keep <- pass_thr & !is_dbl

# Cells failing each threshold (a cell can fail several; `any` is the union). Keys match qc.thresholds.
fails <- data.frame(nFeature_min = !pass_nfeat_min, nFeature_max = !pass_nfeat_max,
                    nCount_min = !pass_ncount_min, nCount_max = !pass_ncount_max,
                    percent_mt_max = !pass_mt, any = !pass_thr)
count_fails <- function(rows) as.list(vapply(fails[rows, , drop = FALSE], sum, integer(1)))

# Cells that fail on the corrected counts but would pass every threshold on the
# uncorrected totals (02 keeps them): ambient correction pushed them out.
has_uncorrected <- all(c("nCount_RNA_uncorrected", "nFeature_RNA_uncorrected") %in% colnames(obj@meta.data))
only_after_correction <- if (has_uncorrected)
  !pass_thr & qc_pass(obj$nFeature_RNA_uncorrected, obj$nCount_RNA_uncorrected, obj$percent.mt, qc) else NULL

reason <- ifelse(!pass_thr & is_dbl, "thresholds+doublet", ifelse(!pass_thr, "thresholds", ifelse(is_dbl, "doublet", NA)))
has_source <- "cell_source" %in% colnames(obj@meta.data)
source_before <- if (has_source) table(obj$sample_id, obj$cell_source) else NULL

# --- Pre-filter UMAP (qc.prefilter_umap): all input cells, coloured by this run's pass/fail ---
prefilter_emb <- NULL
if (isTRUE(cfg$qc$prefilter_umap) && args$make_plots) {
  log_msg("Computing pre-filter UMAP (qc.prefilter_umap = true)")
  prefilter_emb <- tryCatch({
    tmp <- obj
    tmp$qc_status <- ifelse(keep, "kept", ifelse(is.na(reason), "kept", reason))
    tmp <- NormalizeData(tmp, verbose = FALSE)
    tmp <- FindVariableFeatures(tmp, nfeatures = 2000, verbose = FALSE)
    tmp <- ScaleData(tmp, verbose = FALSE)
    npc <- min(20, ncol(tmp) - 1)
    tmp <- RunPCA(tmp, npcs = npc, verbose = FALSE, seed.use = cfg$project$seed)
    RunUMAP(tmp, dims = 1:npc, seed.use = cfg$project$seed, verbose = FALSE)
  }, error = function(e) { warn_msg(sprintf("pre-filter UMAP failed: %s", e$message)); NULL })
}

# --- Split kept / removed ---
n_filtered_out <- sum(!keep)
out_rds_path <- file.path(outdir, "seurat_filtered_out.rds")
out_csv_path <- file.path(outdir, "filtered_out_cells.csv")
if (n_filtered_out > 0) {
  ro <- obj[, !keep]
  for (k in setdiff(names(fails), "any")) ro[[paste0("failed_", k)]] <- fails[[k]][!keep]
  if (has_uncorrected) ro$passes_uncorrected <- only_after_correction[!keep]
  ro$removal_reason <- reason[!keep]
  save_run_rds(ro, cfg, "seurat_filtered_out.rds")
  keep_cols <- intersect(c("sample_id", "removal_reason", paste0("failed_", setdiff(names(fails), "any")),
                           "nCount_RNA", "nFeature_RNA", "percent.mt", "nCount_RNA_uncorrected",
                           "nFeature_RNA_uncorrected", "passes_uncorrected", "decontX_contamination",
                           "scDblFinder_class", "scDblFinder_score", "cell_source"), colnames(ro@meta.data))
  write.csv(data.frame(barcode = colnames(ro), ro@meta.data[, keep_cols, drop = FALSE]),
            out_csv_path, row.names = FALSE)
  rm(ro)
}
obj_all_meta <- obj@meta.data
obj <- obj[, keep]
total_after <- ncol(obj)
after_counts <- table(factor(obj$sample_id, levels = sample_ids))
log_msg(sprintf("Cells after filtering: %d (removed %d)", total_after, n_filtered_out))

pct_removed <- round(100 * (1 - total_after / total_before), 1)
if (pct_removed > 50) warn_msg(sprintf("%.1f%% of cells filtered; review thresholds", pct_removed))
for (sid in sample_ids) if (after_counts[[sid]] == 0)
  emit_error("All cells filtered for a sample", list(sample = sid, thresholds = qc))

rds_path <- save_run_rds(obj, cfg, "seurat_filtered.rds")

# --- JSON summary ---
source_after <- if (has_source) table(obj$sample_id, obj$cell_source) else NULL
by_source <- function(tab, sid) if (is.null(tab)) NULL else as.list(setNames(as.integer(tab[sid, ]), colnames(tab)))
per_sample_summary <- list()
for (sid in sample_ids) {
  r <- obj_all_meta$sample_id == sid
  b <- as.integer(before_counts[sid]); a <- as.integer(after_counts[sid])
  per_sample_summary[[sid]] <- c(
    list(n_cells_before = b,
         removed_by_threshold = count_fails(r)),
    if (has_uncorrected) list(n_removed_only_after_correction = sum(only_after_correction[r])),
    list(n_doublets_called = sum(obj_all_meta$scDblFinder_class[r] == "doublet"),
         n_removed_doublet_only = sum(is_dbl[r] & pass_thr[r]),
         n_cells_after = a,
         pct_removed = round(100 * (b - a) / b, 1)),
    if (has_source) list(n_cells_by_source = list(before = by_source(source_before, sid),
                                                  after = by_source(source_after, sid)))
  )
}

summary_data <- c(
  summary_header("04_filter", cfg),
  list(
    inputs = list(rds = rel_out(in_rds, cfg),
                  doublet_calls = rel_out(file.path(dbl_dir, "doublet_calls.csv"), cfg),
                  doublet_detection = dbl_detection),
    params = list(thresholds = qc, remove_doublets = remove_doublets),
    total = c(
      list(n_cells_before = total_before,
           removed_by_threshold = count_fails(rep(TRUE, nrow(fails)))),
      if (has_uncorrected) list(n_removed_only_after_correction = sum(only_after_correction)),
      list(n_doublets_called = sum(obj_all_meta$scDblFinder_class == "doublet"),
           n_removed_doublet_only = sum(is_dbl & pass_thr),
           n_cells_after = total_after,
           pct_removed = pct_removed),
      if (has_source) list(n_cells_by_source = list(before = as.list(colSums(source_before)),
                                                    after = as.list(colSums(source_after))))
    ),
    per_sample = per_sample_summary,
    outputs = list(rds = rel_out(rds_path, cfg),
                   filtered_out = if (n_filtered_out > 0)
                     list(n_cells = n_filtered_out, rds = rel_out(out_rds_path, cfg), csv = rel_out(out_csv_path, cfg))
                   else list(n_cells = 0L)),
    status = "success"
  )
)
write_json_summary(summary_data, file.path(outdir, STEPS[["04_filter"]]$summary))

# --- Plots ---
if (args$make_plots) {
  log_msg("Generating filter plots")
  md <- obj_all_meta
  md$stage <- "Before filter"
  kept_md <- md[keep, ]
  kept_md$stage <- "After filter"
  comb <- rbind(md, kept_md)
  comb$stage <- factor(comb$stage, levels = c("Before filter", "After filter"))
  violin <- function(y, title, ylab, lines, log = FALSE, file) {
    p <- ggplot(comb, aes(x = sample_id, y = .data[[y]], fill = stage)) +
      geom_violin(alpha = 0.6, scale = "width") +
      qc_threshold_lines("h", lines, color = "red") +
      labs(title = title, y = ylab, x = NULL, fill = NULL) +
      theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    if (log) p <- p + scale_y_log10()
    save_plot(p, file, cfg, width = 8, height = 5)
  }
  violin("nFeature_RNA", "Genes per cell: before/after filtering", "nFeature_RNA",
         c(qc$nFeature_min, qc$nFeature_max), FALSE, "04_filter_nFeature_before_after.png")
  violin("nCount_RNA", "UMI counts per cell: before/after filtering", "nCount_RNA (log10)",
         c(qc$nCount_min, qc$nCount_max), TRUE, "04_filter_nCount_before_after.png")
  violin("percent.mt", "Mitochondrial %: before/after filtering", "% MT",
         qc$percent_mt_max, FALSE, "04_filter_pct_mt_before_after.png")

  # Corrected counts vs genes with this run's thresholds
  p <- ggplot(md, aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) +
    geom_point(size = 0.3, alpha = 0.5) +
    scale_color_viridis_c() +
    facet_wrap(~sample_id) +
    qc_threshold_lines("h", c(qc$nFeature_min, qc$nFeature_max), color = "red", linewidth = 0.3) +
    qc_threshold_lines("v", c(qc$nCount_min, qc$nCount_max), color = "red", linewidth = 0.3) +
    scale_x_log10() + scale_y_log10() +
    labs(title = "UMI counts vs genes detected (corrected counts, this run's thresholds)", color = "% MT")
  save_plot(p, "04_scatter_nCount_nFeature.png", cfg, width = 10, height = 8)

  # Cells kept per sample, removed split into thresholds and doublets
  removal_df <- data.frame(
    sample = sample_ids,
    before = vapply(per_sample_summary, function(x) x$n_cells_before, integer(1)),
    thresholds = vapply(per_sample_summary, function(x) x$removed_by_threshold$any, integer(1)),
    doublets = vapply(per_sample_summary, function(x) x$n_removed_doublet_only, integer(1)),
    after = vapply(per_sample_summary, function(x) x$n_cells_after, integer(1)),
    stringsAsFactors = FALSE)
  removal_long <- rbind(
    data.frame(sample = removal_df$sample, count = removal_df$after, category = "Kept"),
    data.frame(sample = removal_df$sample, count = removal_df$thresholds, category = "Removed: QC thresholds"),
    data.frame(sample = removal_df$sample, count = removal_df$doublets, category = "Removed: doublets"))
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
    scale_x_discrete(expand = expansion(add = 0.9)) +
    labs(title = sprintf("Cells per sample after filtering (total %s kept of %s)",
                         format(sum(removal_df$after), big.mark = ","), format(sum(removal_df$before), big.mark = ",")),
         y = "Cells", x = NULL, fill = NULL) +
    theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p, "04_filter_cells_removed.png", cfg, width = 9, height = 5.5)

  if (!is.null(prefilter_emb)) {
    p <- DimPlot(prefilter_emb, reduction = "umap", group.by = "qc_status", order = TRUE) +
      ggtitle("Pre-filter UMAP — kept vs removal reason (all input cells)")
    save_plot(p, "04_umap_prefilter_qc.png", cfg, width = 9, height = 7)
    if (has_source) {
      p <- DimPlot(prefilter_emb, reduction = "umap", group.by = "cell_source", order = TRUE) +
        ggtitle("Pre-filter UMAP — cell source (all input cells)")
      save_plot(p, "04_umap_prefilter_cell_source.png", cfg, width = 9, height = 7)
    }
  }
}

log_msg("Step 04 complete")
