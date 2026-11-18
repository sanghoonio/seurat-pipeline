#!/usr/bin/env Rscript
# 01_load_qc_metrics.R
# Load Cell Ranger H5 files, merge, compute QC metrics.
# Writes: seurat_raw.rds, qc_metrics.csv, qc_summary.json

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
cfg <- load_config(args$config_path, step = "01_load_qc_metrics", summary_file = "qc_summary.json")
set.seed(cfg$project$seed)

outdir <- output_dir(cfg)

# --- Load samples ---
log_msg(sprintf("Loading %d samples", length(cfg$samples)))

obj_list <- list()
per_sample_stats <- list()

for (s in cfg$samples) {
  filt_path <- resolve_path(sample_filtered_path(s), cfg)

  # If step 00 wrote a called-cell set (emptyDrops), use it: subset the RAW
  # matrix to the chosen barcodes (cell_calling.cells). Otherwise load the
  # Cell Ranger filtered matrix as before.
  bc <- read_called_barcodes(s$id, cfg)
  if (!is.null(bc)) {
    raw_path <- resolve_path(sample_raw_path(s), cfg)
    log_msg(sprintf("  Loading %s: %d called cells from raw %s", s$id, length(bc), raw_path))
    raw_counts <- read_10x_any(raw_path)
    raw_counts <- restrict_to_filtered_features(raw_counts, rownames(read_10x_any(filt_path)), s$id)
    keep <- intersect(bc, colnames(raw_counts))
    counts <- raw_counts[, keep, drop = FALSE]
    rm(raw_counts); gc()
  } else {
    log_msg(sprintf("  Loading %s (filtered) from %s", s$id, filt_path))
    counts <- read_10x_any(filt_path)
  }

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
    source_matrix = if (!is.null(bc)) "raw (called barcodes)" else "filtered",
    n_cells = ncol(obj),
    n_features = nrow(obj),
    n_genes_detected = sum(rowSums(LayerData(obj, layer = "counts")) > 0)
  )

  log_msg(sprintf("  %s: %d cells, %d genes detected",
                   s$id, ncol(obj), per_sample_stats[[s$id]]$n_genes_detected))
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

mt_pattern <- cfg$qc$metrics$mt_pattern
ribo_pattern <- cfg$qc$metrics$ribo_pattern

obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)
obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = ribo_pattern)
obj[["log10_nCount"]] <- log10(obj$nCount_RNA + 1)
obj[["complexity"]] <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)

# --- Cell-calling provenance (cell_source, emptydrops_FDR) from step 00 ---
obj <- attach_cell_calling_meta(obj, cfg)
has_source <- "cell_source" %in% colnames(obj@meta.data)
if (has_source) {
  log_msg(sprintf("Cell source: %s",
    paste(sprintf("%s=%d", names(table(obj$cell_source)), as.integer(table(obj$cell_source))),
          collapse = ", ")))
}

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
if (has_source) {
  qc_df$cell_source <- obj$cell_source
  qc_df$emptydrops_FDR <- signif(obj$emptydrops_FDR, 4)
}
write.csv(qc_df, qc_csv_path, row.names = FALSE)

# --- Per-sample summary stats ---
for (sid in names(per_sample_stats)) {
  cells <- qc_df[qc_df$sample_id == sid, ]
  if (has_source) per_sample_stats[[sid]]$n_cells_by_source <- as.list(table(cells$cell_source))
  per_sample_stats[[sid]]$quantiles <- qc_quantiles(cells)
  per_sample_stats[[sid]]$threshold_preview_pct <- threshold_preview(cells, cfg$qc$thresholds)
}

# --- Save RDS ---
rds_path <- file.path(outdir, "seurat_raw.rds")
log_msg(sprintf("Saving Seurat object to %s", rds_path))
saveRDS(obj, rds_path)

# --- JSON summary ---
summary_data <- c(
  summary_header("01_load_qc_metrics", cfg),
  list(
    inputs = list(cell_set = cfg$cell_calling$cells, matrices = sample_inputs(cfg)),
    params = list(metrics = cfg$qc$metrics, thresholds = cfg$qc$thresholds),
    n_samples_configured = length(cfg$samples),
    n_samples_loaded = length(obj_list),
    samples_loaded = I(names(obj_list)),
    per_sample = per_sample_stats,
    merged = c(
      list(n_cells = ncol(obj), n_features = nrow(obj)),
      if (has_source) list(n_cells_by_source = as.list(table(qc_df$cell_source))),
      list(quantiles = qc_quantiles(qc_df),
           threshold_preview_pct = threshold_preview(qc_df, cfg$qc$thresholds))
    ),
    outputs = list(
      rds = rel_out(rds_path, cfg),
      csv = rel_out(qc_csv_path, cfg)
    ),
    status = "success"
  )
)

write_json_summary(summary_data, file.path(outdir, "qc_summary.json"))

# --- Plots ---
if (args$make_plots) {
  log_msg("Generating QC plots")

  theme_set(theme_bw(base_size = 12))
  mt_max <- cfg$qc$thresholds$percent_mt_max

  # Bar: cells loaded per sample
  n_df <- as.data.frame(table(sample_id = qc_df$sample_id), responseName = "cells")
  p <- ggplot(n_df, aes(x = sample_id, y = cells, fill = sample_id)) +
    geom_col() +
    geom_text(aes(label = format(cells, big.mark = ",")), vjust = -0.3, size = 3.5) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(title = sprintf("Cells loaded per sample (total %s)", format(sum(n_df$cells), big.mark = ",")),
         x = NULL, y = "Cells") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + NoLegend()
  save_plot(p, "01_cells_per_sample.png", cfg, width = 8, height = 5)

  # Violin: nFeature by sample
  p <- VlnPlot(obj, features = "nFeature_RNA", group.by = "sample_id", pt.size = 0) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    qc_threshold_lines("h", c(cfg$qc$thresholds$nFeature_min, cfg$qc$thresholds$nFeature_max), color = "red") +
    ggtitle("Genes per cell") + NoLegend()
  save_plot(p, "01_violin_nFeature.png", cfg, width = 8, height = 5)

  # Violin: nCount by sample
  p <- VlnPlot(obj, features = "nCount_RNA", group.by = "sample_id", pt.size = 0) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    qc_threshold_lines("h", c(cfg$qc$thresholds$nCount_min, cfg$qc$thresholds$nCount_max), color = "red") +
    ggtitle("UMI counts per cell") + NoLegend()
  save_plot(p, "01_violin_nCount.png", cfg, width = 8, height = 5)

  # Violin: percent.mt by sample
  p <- VlnPlot(obj, features = "percent.mt", group.by = "sample_id", pt.size = 0) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    qc_threshold_lines("h", mt_max, color = "red") +
    ggtitle("Mitochondrial %") + NoLegend()
  save_plot(p, "01_violin_pct_mt.png", cfg, width = 8, height = 5)

  # Violin: percent.ribo by sample
  p <- VlnPlot(obj, features = "percent.ribo", group.by = "sample_id", pt.size = 0) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Ribosomal %") + NoLegend()
  save_plot(p, "01_violin_pct_ribo.png", cfg, width = 8, height = 5)

  # Scatter: nCount vs nFeature colored by percent.mt
  p <- ggplot(qc_df, aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) +
    geom_point(size = 0.3, alpha = 0.5) +
    scale_color_viridis_c() +
    facet_wrap(~sample_id) +
    qc_threshold_lines("h", c(cfg$qc$thresholds$nFeature_min, cfg$qc$thresholds$nFeature_max), color = "red", linewidth = 0.3) +
    qc_threshold_lines("v", c(cfg$qc$thresholds$nCount_min, cfg$qc$thresholds$nCount_max), color = "red", linewidth = 0.3) +
    scale_x_log10() + scale_y_log10() +
    labs(title = "UMI counts vs genes detected", color = "% MT")
  save_plot(p, "01_scatter_nCount_nFeature.png", cfg, width = 10, height = 8)

  # Scatter: complexity vs nCount
  p <- ggplot(qc_df, aes(x = nCount_RNA, y = complexity, color = sample_id)) +
    geom_point(size = 0.3, alpha = 0.5) +
    facet_wrap(~sample_id) +
    scale_x_log10() +
    labs(title = "Complexity (log10 genes / log10 UMI) vs UMI counts") +
    NoLegend()
  save_plot(p, "01_scatter_complexity.png", cfg, width = 10, height = 8)

  # Density: percent.mt with threshold
  p <- ggplot(qc_df, aes(x = percent.mt, fill = sample_id)) +
    geom_density(alpha = 0.4) +
    qc_threshold_lines("v", mt_max, color = "red") +
    labs(title = "Mitochondrial % distribution", x = "% MT", fill = "Sample")
  save_plot(p, "01_density_pct_mt.png", cfg, width = 8, height = 5)

  # Cell-calling provenance plots (only when step 00 tagged cells)
  if (has_source) {
    p <- ggplot(qc_df, aes(x = nCount_RNA, y = nFeature_RNA, color = cell_source)) +
      geom_point(size = 0.4, alpha = 0.6) +
      scale_x_log10() + scale_y_log10() +
      qc_threshold_lines("h", c(cfg$qc$thresholds$nFeature_min, cfg$qc$thresholds$nFeature_max), color = "grey30", linewidth = 0.3) +
      qc_threshold_lines("v", c(cfg$qc$thresholds$nCount_min, cfg$qc$thresholds$nCount_max), color = "grey30", linewidth = 0.3) +
      labs(title = "Counts vs genes by cell source (emptyDrops vs Cell Ranger filtered)",
           color = "Cell source")
    save_plot(p, "01_scatter_cell_source.png", cfg, width = 9, height = 7)

    p <- ggplot(qc_df, aes(x = percent.mt, fill = cell_source)) +
      geom_density(alpha = 0.4) +
      qc_threshold_lines("v", mt_max, color = "red") +
      labs(title = "Mitochondrial % by cell source", x = "% MT", fill = "Cell source")
    save_plot(p, "01_density_pct_mt_by_source.png", cfg, width = 8, height = 5)
  }

  # Optional PRE-FILTER UMAP (opt-in via qc.prefilter_umap): embeds ALL loaded
  # cells — including those step 03 will drop — so you can see where QC-fail and
  # emptyDrops-recovered cells sit in transcriptional space (the step-04 UMAP
  # only shows post-filter survivors).
  if (isTRUE(cfg$qc$prefilter_umap)) {
    log_msg("Computing pre-filter UMAP (qc.prefilter_umap = true)")
    qc <- cfg$qc$thresholds
    obj$qc_pass <- ifelse(
      qc_pass(obj$nFeature_RNA, obj$nCount_RNA, obj$percent.mt, qc), "pass", "fail")
    emb <- tryCatch({
      tmp <- obj
      if (packageVersion("Seurat") >= "5.0.0") tmp[["RNA"]] <- JoinLayers(tmp[["RNA"]])
      tmp <- NormalizeData(tmp, verbose = FALSE)
      tmp <- FindVariableFeatures(tmp, nfeatures = 2000, verbose = FALSE)
      tmp <- ScaleData(tmp, verbose = FALSE)
      npc <- min(20, ncol(tmp) - 1)
      tmp <- RunPCA(tmp, npcs = npc, verbose = FALSE, seed.use = cfg$project$seed)
      RunUMAP(tmp, dims = 1:npc, seed.use = cfg$project$seed, verbose = FALSE)
    }, error = function(e) { log_msg(sprintf("  pre-filter UMAP failed: %s", e$message)); NULL })
    if (!is.null(emb)) {
      p <- DimPlot(emb, reduction = "umap", group.by = "qc_pass", order = TRUE) +
        scale_color_manual(values = c(pass = "grey70", fail = "firebrick")) +
        ggtitle("Pre-filter UMAP — QC pass/fail (all loaded cells)")
      save_plot(p, "01_umap_prefilter_qc.png", cfg, width = 9, height = 7)
      if (has_source) {
        p <- DimPlot(emb, reduction = "umap", group.by = "cell_source", order = TRUE) +
          ggtitle("Pre-filter UMAP — cell source (all loaded cells)")
        save_plot(p, "01_umap_prefilter_cell_source.png", cfg, width = 9, height = 7)
      }
    }
  }
}

log_msg("Step 01 complete")
