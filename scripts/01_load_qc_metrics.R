#!/usr/bin/env Rscript
# 01_load_qc_metrics.R
# Load Cell Ranger H5 files, merge, compute QC metrics.
# Writes: seurat_raw.rds, qc_metrics.csv, qc_summary.json

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})
source("scripts/utils.R")

args <- parse_args()
cfg <- load_config(args$config_path)
set.seed(cfg$project$seed)

outdir <- output_dir(cfg)

# --- Load samples ---
log_msg(sprintf("Loading %d samples", length(cfg$samples)))

obj_list <- list()
per_sample_stats <- list()

for (s in cfg$samples) {
  h5_path <- resolve_path(s$filtered_h5, cfg)

  if (!file.exists(h5_path)) {
    emit_error("Filtered H5 file not found", list(sample = s$id, path = h5_path))
  }

  log_msg(sprintf("  Loading %s from %s", s$id, h5_path))

  counts <- tryCatch(
    Read10X_h5(h5_path, use.names = TRUE, unique.features = TRUE),
    error = function(e) {
      emit_error("Failed to read H5 file", list(sample = s$id, path = h5_path, error = e$message))
    }
  )

  obj <- CreateSeuratObject(
    counts = counts,
    project = s$id,
    min.cells = 0,
    min.features = 0
  )

  obj <- attach_sample_metadata(obj, s)
  obj <- RenameCells(obj, add.cell.id = s$id)
  obj_list[[s$id]] <- obj

  per_sample_stats[[s$id]] <- list(
    cells_loaded = ncol(obj),
    genes_detected = sum(rowSums(LayerData(obj, layer = "counts")) > 0)
  )

  log_msg(sprintf("  %s: %d cells, %d genes detected",
                   s$id, ncol(obj), per_sample_stats[[s$id]]$genes_detected))
}

# --- Merge ---
log_msg("Merging samples")

if (length(obj_list) == 1) {
  obj <- obj_list[[1]]
} else {
  obj <- merge(obj_list[[1]], y = obj_list[2:length(obj_list)])
}

log_msg(sprintf("Merged object: %d cells x %d genes", ncol(obj), nrow(obj)))

# --- QC metrics ---
log_msg("Computing QC metrics")

mt_pattern <- cfg$qc$mt_pattern %||% "^MT-"
ribo_pattern <- cfg$qc$ribo_pattern %||% "^RP[SL]"

obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)
obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = ribo_pattern)
obj[["log10_nCount"]] <- log10(obj$nCount_RNA + 1)
obj[["complexity"]] <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)

# --- Export QC metrics CSV ---
qc_csv_path <- file.path(outdir, "qc_metrics.csv")
log_msg(sprintf("Writing QC metrics to %s", qc_csv_path))

qc_df <- data.frame(
  barcode = colnames(obj),
  sample_id = obj$sample_id,
  nCount_RNA = obj$nCount_RNA,
  nFeature_RNA = obj$nFeature_RNA,
  percent.mt = round(obj$percent.mt, 3),
  percent.ribo = round(obj$percent.ribo, 3),
  log10_nCount = round(obj$log10_nCount, 3),
  complexity = round(obj$complexity, 3),
  stringsAsFactors = FALSE
)
write.csv(qc_df, qc_csv_path, row.names = FALSE)

# --- Per-sample summary stats ---
for (sid in names(per_sample_stats)) {
  cells <- qc_df[qc_df$sample_id == sid, ]
  per_sample_stats[[sid]]$median_nFeature <- median(cells$nFeature_RNA)
  per_sample_stats[[sid]]$median_nCount <- median(cells$nCount_RNA)
  per_sample_stats[[sid]]$median_percent_mt <- round(median(cells$percent.mt), 3)
  per_sample_stats[[sid]]$median_percent_ribo <- round(median(cells$percent.ribo), 3)
  per_sample_stats[[sid]]$median_complexity <- round(median(cells$complexity, na.rm = TRUE), 3)
}

# --- Save RDS ---
rds_path <- file.path(outdir, "seurat_raw.rds")
log_msg(sprintf("Saving Seurat object to %s", rds_path))
saveRDS(obj, rds_path)

# --- JSON summary ---
mt_vals <- obj$percent.mt
summary_data <- c(
  summary_header("01_load_qc_metrics", cfg),
  list(
    samples_loaded = length(cfg$samples),
    per_sample = per_sample_stats,
    merged = list(
      total_cells = ncol(obj),
      total_genes = nrow(obj),
      median_nFeature = median(obj$nFeature_RNA),
      median_nCount = median(obj$nCount_RNA),
      percent_mt_quantiles = list(
        q25 = round(quantile(mt_vals, 0.25), 3),
        q50 = round(quantile(mt_vals, 0.50), 3),
        q75 = round(quantile(mt_vals, 0.75), 3),
        q95 = round(quantile(mt_vals, 0.95), 3),
        q99 = round(quantile(mt_vals, 0.99), 3)
      )
    ),
    outputs = list(
      rds = rds_path,
      csv = qc_csv_path
    ),
    status = "success"
  )
)

write_json_summary(summary_data, file.path(outdir, "qc_summary.json"))

# --- Plots ---
if (args$make_plots) {
  log_msg("Generating QC plots")

  theme_set(theme_bw(base_size = 12))
  mt_max <- cfg$qc$percent_mt_max

  # Violin: nFeature by sample
  p <- VlnPlot(obj, features = "nFeature_RNA", group.by = "sample_id", pt.size = 0) +
    geom_hline(yintercept = c(cfg$qc$nFeature_min, cfg$qc$nFeature_max),
               linetype = "dashed", color = "red") +
    ggtitle("Genes per cell") + NoLegend()
  save_plot(p, "01_violin_nFeature.png", cfg, width = 8, height = 5)

  # Violin: nCount by sample
  p <- VlnPlot(obj, features = "nCount_RNA", group.by = "sample_id", pt.size = 0) +
    geom_hline(yintercept = c(cfg$qc$nCount_min, cfg$qc$nCount_max),
               linetype = "dashed", color = "red") +
    ggtitle("UMI counts per cell") + NoLegend()
  save_plot(p, "01_violin_nCount.png", cfg, width = 8, height = 5)

  # Violin: percent.mt by sample
  p <- VlnPlot(obj, features = "percent.mt", group.by = "sample_id", pt.size = 0) +
    geom_hline(yintercept = mt_max, linetype = "dashed", color = "red") +
    ggtitle("Mitochondrial %") + NoLegend()
  save_plot(p, "01_violin_pct_mt.png", cfg, width = 8, height = 5)

  # Violin: percent.ribo by sample
  p <- VlnPlot(obj, features = "percent.ribo", group.by = "sample_id", pt.size = 0) +
    ggtitle("Ribosomal %") + NoLegend()
  save_plot(p, "01_violin_pct_ribo.png", cfg, width = 8, height = 5)

  # Scatter: nCount vs nFeature colored by percent.mt
  p <- ggplot(qc_df, aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) +
    geom_point(size = 0.3, alpha = 0.5) +
    scale_color_viridis_c() +
    facet_wrap(~sample_id) +
    geom_hline(yintercept = c(cfg$qc$nFeature_min, cfg$qc$nFeature_max),
               linetype = "dashed", color = "red", linewidth = 0.3) +
    geom_vline(xintercept = c(cfg$qc$nCount_min, cfg$qc$nCount_max),
               linetype = "dashed", color = "red", linewidth = 0.3) +
    labs(title = "UMI counts vs genes detected", color = "% MT")
  save_plot(p, "01_scatter_nCount_nFeature.png", cfg, width = 10, height = 8)

  # Scatter: complexity vs nCount
  p <- ggplot(qc_df, aes(x = nCount_RNA, y = complexity, color = sample_id)) +
    geom_point(size = 0.3, alpha = 0.5) +
    facet_wrap(~sample_id) +
    labs(title = "Complexity (log10 genes / log10 UMI) vs UMI counts") +
    NoLegend()
  save_plot(p, "01_scatter_complexity.png", cfg, width = 10, height = 8)

  # Density: percent.mt with threshold
  p <- ggplot(qc_df, aes(x = percent.mt, fill = sample_id)) +
    geom_density(alpha = 0.4) +
    geom_vline(xintercept = mt_max, linetype = "dashed", color = "red") +
    labs(title = "Mitochondrial % distribution", x = "% MT", fill = "Sample")
  save_plot(p, "01_density_pct_mt.png", cfg, width = 8, height = 5)
}

log_msg("Step 01 complete")
