#!/usr/bin/env Rscript
# 02_soupx.R
# Estimate and remove ambient RNA contamination using SoupX.
# Runs per-sample: needs both raw and filtered matrices to estimate the soup.
# Rebuilds the Seurat object from corrected counts rather than patching the merged object.
# Writes: seurat_decontaminated.rds, qc_metrics_decontaminated.csv, soupx_summary.json

suppressPackageStartupMessages({
  library(Seurat)
  library(SoupX)
  library(ggplot2)
})
source("scripts/utils.R")

args <- parse_args()
cfg <- load_config(args$config_path)
set.seed(cfg$project$seed)

outdir <- output_dir(cfg)

# --- Load raw merged object (for metadata) ---
raw_rds <- file.path(outdir, "seurat_raw.rds")
if (!file.exists(raw_rds)) {
  emit_error("Input RDS not found. Run step 01 first.", list(expected = raw_rds))
}
log_msg("Loading seurat_raw.rds")
obj <- readRDS(raw_rds)

soupx_cfg <- cfg$soupx %||% list(enabled = TRUE)
if (!isTRUE(soupx_cfg$enabled)) {
  log_msg("SoupX disabled in config, copying raw to decontaminated")
  if (packageVersion("Seurat") >= "5.0.0") {
    obj <- JoinLayers(obj)
  }
  saveRDS(obj, file.path(outdir, "seurat_decontaminated.rds"))
  write_json_summary(
    c(summary_header("02_soupx", cfg), list(skipped = TRUE, reason = "disabled in config", status = "success")),
    file.path(outdir, "soupx_summary.json")
  )
  log_msg("Step 02 complete (skipped)")
  quit(status = 0, save = "no")
}

contam_range <- soupx_cfg$contamination_range %||% c(0.01, 0.8)
round_to_int <- soupx_cfg$round_to_int %||% TRUE

# --- Run SoupX per sample, collect corrected counts ---
log_msg("Running SoupX per sample")

per_sample_results <- list()
corrected_obj_list <- list()

for (s in cfg$samples) {
  sid <- s$id
  log_msg(sprintf("  Processing %s", sid))

  filtered_h5 <- resolve_path(s$filtered_h5, cfg)
  raw_h5 <- if (!is.null(s$raw_h5)) resolve_path(s$raw_h5, cfg) else NA_character_

  if (is.na(raw_h5) || !file.exists(raw_h5)) {
    log_msg(sprintf("  WARNING: Raw H5 not found for %s, skipping SoupX", sid))
    per_sample_results[[sid]] <- list(status = "skipped", reason = "raw H5 not found")
    # Fall back to uncorrected counts so the sample is still included
    filt_counts <- Read10X_h5(filtered_h5, use.names = TRUE, unique.features = TRUE)
    sample_obj <- CreateSeuratObject(counts = filt_counts, project = sid,
                                      min.cells = 0, min.features = 0)
    sample_obj <- attach_sample_metadata(sample_obj, s)
    sample_obj <- RenameCells(sample_obj, add.cell.id = sid)
    corrected_obj_list[[sid]] <- sample_obj
    next
  }

  # Load matrices fresh from H5 (avoids v5 layer issues)
  filt_counts <- Read10X_h5(filtered_h5, use.names = TRUE, unique.features = TRUE)
  raw_counts <- Read10X_h5(raw_h5, use.names = TRUE, unique.features = TRUE)

  # Ensure gene names match
  common_genes <- intersect(rownames(filt_counts), rownames(raw_counts))
  filt_counts <- filt_counts[common_genes, ]
  raw_counts <- raw_counts[common_genes, ]

  # Build SoupChannel and estimate contamination
  sc <- tryCatch({
    sc <- SoupChannel(raw_counts, filt_counts)

    # Quick clustering for SoupX contamination estimation
    tmp <- CreateSeuratObject(filt_counts, min.cells = 0, min.features = 0)
    tmp <- NormalizeData(tmp, verbose = FALSE)
    tmp <- FindVariableFeatures(tmp, verbose = FALSE)
    tmp <- ScaleData(tmp, verbose = FALSE)
    n_pcs <- min(10, ncol(tmp) - 1)
    tmp <- RunPCA(tmp, npcs = n_pcs, verbose = FALSE)
    tmp <- FindNeighbors(tmp, dims = 1:n_pcs, verbose = FALSE)
    tmp <- FindClusters(tmp, resolution = 0.5, verbose = FALSE)

    sc <- setClusters(sc, setNames(as.character(Idents(tmp)), colnames(tmp)))
    sc <- autoEstCont(sc, forceAccept = TRUE)
    sc
  }, error = function(e) {
    log_msg(sprintf("  WARNING: SoupX estimation failed for %s: %s", sid, e$message))
    NULL
  })

  if (is.null(sc)) {
    per_sample_results[[sid]] <- list(status = "failed", reason = "estimation failed")
    corrected_counts <- filt_counts
  } else {
    rho <- sc$fit$rhoEst

    corrected_counts <- tryCatch({
      adjustCounts(sc, roundToInt = round_to_int)
    }, error = function(e) {
      log_msg(sprintf("  WARNING: adjustCounts failed for %s: %s", sid, e$message))
      NULL
    })

    if (is.null(corrected_counts)) {
      per_sample_results[[sid]] <- list(status = "failed", reason = "adjustCounts failed")
      corrected_counts <- filt_counts
    } else {
      n_changed <- sum(corrected_counts != filt_counts)
      n_total <- length(filt_counts@x)

      per_sample_results[[sid]] <- list(
        status = "success",
        estimated_contamination = round(rho, 4),
        cells_processed = ncol(corrected_counts),
        entries_changed_pct = round(100 * n_changed / max(n_total, 1), 2)
      )
      log_msg(sprintf("  %s: estimated contamination = %.3f, %d cells corrected",
                       sid, rho, ncol(corrected_counts)))
    }
  }

  # Build a per-sample Seurat object from corrected counts
  sample_obj <- CreateSeuratObject(
    counts = corrected_counts,
    project = sid,
    min.cells = 0,
    min.features = 0
  )
  sample_obj <- attach_sample_metadata(sample_obj, s)
  sample_obj <- RenameCells(sample_obj, add.cell.id = sid)

  corrected_obj_list[[sid]] <- sample_obj
}

# --- Merge corrected objects ---
log_msg("Merging corrected samples")

if (length(corrected_obj_list) == 1) {
  obj_new <- corrected_obj_list[[1]]
} else {
  obj_new <- merge(corrected_obj_list[[1]], y = corrected_obj_list[2:length(corrected_obj_list)])
}

if (packageVersion("Seurat") >= "5.0.0") {
  obj_new <- JoinLayers(obj_new)
}

log_msg(sprintf("Merged decontaminated object: %d cells x %d genes", ncol(obj_new), nrow(obj_new)))

# --- Recalculate QC metrics ---
mt_pattern <- cfg$qc$mt_pattern %||% "^MT-"
ribo_pattern <- cfg$qc$ribo_pattern %||% "^RP[SL]"

obj_new[["percent.mt"]] <- PercentageFeatureSet(obj_new, pattern = mt_pattern)
obj_new[["percent.ribo"]] <- PercentageFeatureSet(obj_new, pattern = ribo_pattern)
obj_new[["log10_nCount"]] <- log10(obj_new$nCount_RNA + 1)
obj_new[["complexity"]] <- log10(obj_new$nFeature_RNA) / log10(obj_new$nCount_RNA)

# --- Save ---
rds_path <- file.path(outdir, "seurat_decontaminated.rds")
log_msg(sprintf("Saving decontaminated object to %s", rds_path))
saveRDS(obj_new, rds_path)

# --- Updated QC metrics CSV ---
qc_csv_path <- file.path(outdir, "qc_metrics_decontaminated.csv")
qc_df <- data.frame(
  barcode = colnames(obj_new),
  sample_id = obj_new$sample_id,
  nCount_RNA = obj_new$nCount_RNA,
  nFeature_RNA = obj_new$nFeature_RNA,
  percent.mt = round(obj_new$percent.mt, 3),
  percent.ribo = round(obj_new$percent.ribo, 3),
  log10_nCount = round(obj_new$log10_nCount, 3),
  complexity = round(obj_new$complexity, 3),
  stringsAsFactors = FALSE
)
write.csv(qc_df, qc_csv_path, row.names = FALSE)

# --- JSON summary ---
successful_samples <- names(per_sample_results)[
  sapply(per_sample_results, function(x) identical(x$status, "success"))
]

summary_data <- c(
  summary_header("02_soupx", cfg),
  list(
    samples_processed = length(cfg$samples),
    samples_successful = length(successful_samples),
    per_sample = per_sample_results,
    merged_cells = ncol(obj_new),
    merged_genes = nrow(obj_new),
    outputs = list(
      rds = rds_path,
      csv = qc_csv_path
    ),
    status = if (length(successful_samples) > 0) "success" else "warning"
  )
)

write_json_summary(summary_data, file.path(outdir, "soupx_summary.json"))

# --- Plots ---
if (args$make_plots) {
  log_msg("Generating SoupX plots")

  pre_csv <- file.path(outdir, "qc_metrics.csv")
  if (file.exists(pre_csv)) {
    pre_df <- read.csv(pre_csv, stringsAsFactors = FALSE)
    pre_df$stage <- "Before SoupX"
    post_df <- qc_df
    post_df$stage <- "After SoupX"
    comb <- rbind(pre_df[, c("sample_id", "nCount_RNA", "nFeature_RNA", "stage")],
                  post_df[, c("sample_id", "nCount_RNA", "nFeature_RNA", "stage")])
    comb$stage <- factor(comb$stage, levels = c("Before SoupX", "After SoupX"))

    p <- ggplot(comb, aes(x = sample_id, y = nCount_RNA, fill = stage)) +
      geom_violin(alpha = 0.6, scale = "width") +
      labs(title = "UMI counts before/after SoupX", y = "nCount_RNA", x = NULL, fill = NULL) +
      theme_bw(base_size = 12)
    save_plot(p, "02_soupx_nCount_before_after.png", cfg, width = 8, height = 5)

    p <- ggplot(comb, aes(x = sample_id, y = nFeature_RNA, fill = stage)) +
      geom_violin(alpha = 0.6, scale = "width") +
      labs(title = "Genes detected before/after SoupX", y = "nFeature_RNA", x = NULL, fill = NULL) +
      theme_bw(base_size = 12)
    save_plot(p, "02_soupx_nFeature_before_after.png", cfg, width = 8, height = 5)
  }

  contam_df <- data.frame(
    sample = names(per_sample_results),
    contamination = sapply(per_sample_results, function(x) x$estimated_contamination %||% NA),
    stringsAsFactors = FALSE
  )
  contam_df <- contam_df[!is.na(contam_df$contamination), ]

  if (nrow(contam_df) > 0) {
    p <- ggplot(contam_df, aes(x = sample, y = contamination)) +
      geom_col(fill = "steelblue") +
      geom_text(aes(label = sprintf("%.1f%%", contamination * 100)), vjust = -0.3) +
      labs(title = "Estimated ambient RNA contamination (SoupX)", y = "Contamination fraction", x = NULL) +
      theme_bw(base_size = 12) +
      ylim(0, max(contam_df$contamination) * 1.3)
    save_plot(p, "02_soupx_contamination.png", cfg, width = 6, height = 5)
  }
}

log_msg("Step 02 complete")
