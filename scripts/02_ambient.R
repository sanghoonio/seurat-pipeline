#!/usr/bin/env Rscript
# 02_ambient.R
# Estimate and remove ambient RNA contamination, per sample, with SoupX (one
# contamination fraction per sample) or DecontX (one per cell); ambient.method.
# Both use the raw matrix: SoupX for the soup profile, DecontX for background droplets.
# Rebuilds the Seurat object from corrected counts rather than patching the merged object.
# Writes: seurat_decontaminated.rds, qc_metrics_decontaminated.csv, ambient_summary.json

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
cfg <- load_config(args$config_path, step = "02_ambient", summary_file = "ambient_summary.json")
set.seed(cfg$project$seed)

outdir <- output_dir(cfg)

# --- Load raw merged object (for metadata) ---
raw_rds <- file.path(outdir, "seurat_raw.rds")
if (!file.exists(raw_rds)) {
  emit_error("Input RDS not found. Run step 01 first.", list(expected = raw_rds))
}
log_msg("Loading seurat_raw.rds")
obj <- readRDS(raw_rds)

amb_cfg <- cfg$ambient   # defaults filled by load_config()
if (!isTRUE(amb_cfg$enabled)) {
  log_msg("Ambient correction disabled in config, copying raw to decontaminated")
  if (packageVersion("Seurat") >= "5.0.0") {
    obj <- JoinLayers(obj)
  }
  saveRDS(obj, file.path(outdir, "seurat_decontaminated.rds"))
  write_json_summary(
    c(summary_header("02_ambient", cfg), list(status = "skipped", reason = "ambient.enabled is false")),
    file.path(outdir, "ambient_summary.json")
  )
  log_msg("Step 02 complete (skipped)")
  quit(status = 0, save = "no")
}

method <- amb_cfg$method   # validated by load_config()
pkg <- c(soupx = "SoupX", decontx = "celda")[[method]]
if (!requireNamespace(pkg, quietly = TRUE))
  emit_error(sprintf("ambient.method is '%s' but the '%s' package is not installed", method, pkg),
             list(fix = if (method == "decontx") "BiocManager::install('celda')" else "install.packages('SoupX')"))
log_msg(sprintf("Ambient RNA method: %s", method))
dx_cfg <- amb_cfg$decontx %||% list()

contam_range <- unlist(amb_cfg$soupx$contamination_range %||% c(0.01, 0.8))
round_to_int <- amb_cfg$round_to_int %||% TRUE
# Samples to leave uncorrected (by sample id), e.g. Flex libraries where autoEstCont
# has no soup markers to work with and returns nothing or a pinned-at-ceiling rho.
skip_samples <- as.character(amb_cfg$skip_samples %||% character(0))
if (length(skip_samples)) log_msg(sprintf("Ambient correction skipped for samples: %s", paste(skip_samples, collapse = ", ")))

# --- Method helpers ---
quick_clusters <- function(cell_counts) {
  tmp <- CreateSeuratObject(cell_counts, min.cells = 0, min.features = 0)
  tmp <- NormalizeData(tmp, verbose = FALSE)
  tmp <- FindVariableFeatures(tmp, verbose = FALSE)
  tmp <- ScaleData(tmp, verbose = FALSE)
  n_pcs <- min(10, ncol(tmp) - 1)
  tmp <- RunPCA(tmp, npcs = n_pcs, verbose = FALSE)
  tmp <- FindNeighbors(tmp, dims = 1:n_pcs, verbose = FALSE)
  tmp <- FindClusters(tmp, resolution = 0.5, verbose = FALSE)
  setNames(as.character(Idents(tmp)), colnames(tmp))
}

# Each returns list(corrected, contamination, cell_contamination, extra, reason);
# corrected = NULL on failure.
run_soupx <- function(raw_counts, cell_counts, clusters, contam_range, round_to_int, sid) {
  sc <- tryCatch({
    sc <- SoupX::SoupChannel(raw_counts, cell_counts)
    sc <- SoupX::setClusters(sc, clusters)
    # contaminationRange caps forceAccept, which otherwise accepts absurd estimates
    # (e.g. 0.8) on low-RNA cell sets. doPlot = FALSE avoids a stray Rplots.pdf.
    SoupX::autoEstCont(sc, forceAccept = TRUE, contaminationRange = contam_range, doPlot = FALSE)
  }, error = function(e) {
    warn_msg(sprintf("%s: SoupX estimation failed: %s", sid, e$message)); NULL
  })
  if (is.null(sc)) return(list(corrected = NULL, reason = "estimation failed"))
  rho <- sc$fit$rhoEst
  corrected <- tryCatch(SoupX::adjustCounts(sc, roundToInt = round_to_int), error = function(e) {
    warn_msg(sprintf("%s: adjustCounts failed: %s", sid, e$message)); NULL
  })
  if (is.null(corrected)) return(list(corrected = NULL, reason = "adjustCounts failed"))
  # Flag estimates on a configured bound: at the ceiling the fit is usually
  # unreliable; at the floor the correction is close to a no-op.
  at_bound <- "none"
  if (rho >= contam_range[2] - 1e-3) at_bound <- "upper"
  if (rho <= contam_range[1] + 1e-3) at_bound <- "lower"
  if (at_bound != "none")
    warn_msg(sprintf("%s: contamination %.3f is at the %s bound of contamination_range [%s]",
                    sid, rho, at_bound, paste(contam_range, collapse = ", ")))
  list(corrected = corrected, contamination = rho, cell_contamination = NULL,
       extra = list(at_contamination_bound = at_bound))
}

run_decontx <- function(raw_counts, cell_counts, clusters, dx_cfg, round_to_int, seed, sid) {
  # Background = raw-matrix droplets that are not called cells; they anchor the
  # contamination profile. Without it, contamination is modelled only from the
  # other cell populations.
  bg <- NULL
  if (isTRUE(dx_cfg$use_background %||% TRUE)) {
    bg_bc <- setdiff(colnames(raw_counts), colnames(cell_counts))
    bg <- raw_counts[, bg_bc, drop = FALSE]
    bg <- bg[, Matrix::colSums(bg) > 0, drop = FALSE]
    log_msg(sprintf("  %s: DecontX background = %d non-cell droplets", sid, ncol(bg)))
  }
  fit <- tryCatch(
    celda::decontX(cell_counts, z = clusters[colnames(cell_counts)], background = bg,
                   maxIter = dx_cfg$max_iter %||% 500,
                   delta = unlist(dx_cfg$delta %||% c(10, 10)),
                   estimateDelta = dx_cfg$estimate_delta %||% TRUE,
                   seed = seed, verbose = FALSE),
    error = function(e) {
      warn_msg(sprintf("%s: DecontX failed: %s", sid, e$message)); NULL
    })
  if (is.null(fit)) return(list(corrected = NULL, reason = "decontX failed"))
  corrected <- fit$decontXcounts
  if (round_to_int) corrected <- round(corrected)
  corrected <- Matrix::drop0(as(corrected, "CsparseMatrix"))   # rounding leaves explicit zeros
  contam <- setNames(fit$contamination, colnames(cell_counts))
  list(corrected = corrected, contamination = median(contam), cell_contamination = contam,
       extra = list(contamination_mean = round(mean(contam), 4),
                    contamination_q90 = round(unname(quantile(contam, 0.9)), 4),
                    background_droplets = if (is.null(bg)) 0L else ncol(bg)))
}

# --- Correct each sample, collect corrected counts ---
log_msg(sprintf("Running %s per sample", method))

per_sample_results <- list()
corrected_obj_list <- list()

for (s in cfg$samples) {
  sid <- s$id
  log_msg(sprintf("  Processing %s", sid))

  filtered_path <- resolve_path(sample_filtered_path(s), cfg)
  raw_field <- sample_raw_path(s)
  raw_path <- if (!is.null(raw_field)) resolve_path(raw_field, cfg) else NA_character_
  raw_exists <- !is.na(raw_path) && (file.exists(raw_path) || dir.exists(raw_path))

  # Cell set to decontaminate: emptyDrops called set (step 00) if present,
  # else Cell Ranger filtered.
  called_bc <- read_called_barcodes(sid, cfg)

  if (!raw_exists) {
    warn_msg(sprintf("%s: raw matrix not found, counts left uncorrected", sid))
    per_sample_results[[sid]] <- list(status = "skipped", method = method, reason = "raw matrix not found")
    # Fall back to uncorrected counts so the sample is still included
    cell_counts <- read_10x_any(filtered_path)
    sample_obj <- CreateSeuratObject(counts = cell_counts, project = sid,
                                      min.cells = 0, min.features = 0)
    sample_obj <- attach_sample_metadata(sample_obj, s)
    sample_obj <- RenameCells(sample_obj, add.cell.id = sid)
    corrected_obj_list[[sid]] <- sample_obj
    next
  }

  # Load raw once; derive the cell matrix from it (called set) or from filtered.
  raw_counts <- read_10x_any(raw_path)
  filt_counts <- read_10x_any(filtered_path)
  # Same feature restriction as steps 00/01: raw contributes barcodes only.
  raw_counts <- restrict_to_filtered_features(raw_counts, rownames(filt_counts), sid)
  if (!is.null(called_bc)) {
    cell_counts <- raw_counts[, intersect(called_bc, colnames(raw_counts)), drop = FALSE]
  } else {
    cell_counts <- filt_counts
  }
  rm(filt_counts)

  # Ensure gene names match
  common_genes <- intersect(rownames(cell_counts), rownames(raw_counts))
  cell_counts <- cell_counts[common_genes, ]
  raw_counts <- raw_counts[common_genes, ]

  # Per-sample skip: keep the called-cell counts uncorrected.
  if (sid %in% skip_samples) {
    log_msg(sprintf("  %s: listed in ambient.skip_samples, counts left uncorrected", sid))
    per_sample_results[[sid]] <- list(status = "skipped", method = method, reason = "listed in ambient.skip_samples",
                                      n_cells = ncol(cell_counts))
    sample_obj <- CreateSeuratObject(counts = cell_counts, project = sid, min.cells = 0, min.features = 0)
    sample_obj <- attach_sample_metadata(sample_obj, s)
    sample_obj <- RenameCells(sample_obj, add.cell.id = sid)
    corrected_obj_list[[sid]] <- sample_obj
    rm(raw_counts); gc()
    next
  }

  # Quick per-sample clustering: SoupX uses it to find soup-marker genes, DecontX
  # to define the cell populations contamination is modelled as coming from.
  clusters <- tryCatch(quick_clusters(cell_counts), error = function(e) {
    warn_msg(sprintf("%s: quick clustering failed: %s", sid, e$message)); NULL
  })

  res <- if (is.null(clusters)) {
    list(corrected = NULL, reason = "quick clustering failed")
  } else if (method == "decontx") {
    run_decontx(raw_counts, cell_counts, clusters, dx_cfg, round_to_int, cfg$project$seed, sid)
  } else {
    run_soupx(raw_counts, cell_counts, clusters, contam_range, round_to_int, sid)
  }

  cell_contam <- NULL
  if (is.null(res$corrected)) {
    per_sample_results[[sid]] <- list(status = "failed", method = method, reason = res$reason)
    corrected_counts <- cell_counts
  } else {
    corrected_counts <- res$corrected
    cell_contam <- res$cell_contamination
    n_changed <- sum(corrected_counts != cell_counts)
    n_total <- length(cell_counts@x)

    # How much was removed, and from which genes. The top genes are the main
    # sanity check: they should be plausible ambient transcripts (highly
    # expressed, cell-type-restricted), not housekeeping genes across the board.
    umi_before <- sum(cell_counts)
    umi_after <- sum(corrected_counts)
    gene_before <- Matrix::rowSums(cell_counts)
    removed_by_gene <- gene_before[rownames(corrected_counts)] - Matrix::rowSums(corrected_counts)
    top <- head(sort(removed_by_gene[removed_by_gene > 0], decreasing = TRUE), 20)
    top_removed <- lapply(names(top), function(g) list(
      gene = g,
      umis_removed = round(unname(top[[g]])),
      pct_of_gene = round(100 * unname(top[[g]]) / max(gene_before[[g]], 1), 1)
    ))

    per_sample_results[[sid]] <- c(
      list(status = "success", method = method,
           estimated_contamination = round(res$contamination, 4),
           contamination_basis = if (method == "decontx") "median of per-cell estimates" else "one estimate per sample"),
      res$extra,
      list(n_cells = ncol(corrected_counts),
           entries_changed_pct = round(100 * n_changed / max(n_total, 1), 2),
           umis_removed_pct = round(100 * (1 - umi_after / max(umi_before, 1)), 2),
           median_genes_lost_per_cell = median(Matrix::colSums(cell_counts > 0)[colnames(corrected_counts)] -
                                               Matrix::colSums(corrected_counts > 0)),
           top_removed_genes = top_removed)
    )
    log_msg(sprintf("  %s: %s contamination = %.3f, %.1f%% of UMIs removed, %d cells corrected",
                    sid, method, res$contamination, per_sample_results[[sid]]$umis_removed_pct,
                    ncol(corrected_counts)))
    log_msg(sprintf("    top removed genes: %s", paste(head(names(top), 10), collapse = ", ")))
  }

  # Build a per-sample Seurat object from corrected counts
  sample_obj <- CreateSeuratObject(
    counts = corrected_counts,
    project = sid,
    min.cells = 0,
    min.features = 0
  )
  if (!is.null(cell_contam)) sample_obj$decontX_contamination <- unname(cell_contam[colnames(sample_obj)])
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
mt_pattern <- cfg$qc$metrics$mt_pattern
ribo_pattern <- cfg$qc$metrics$ribo_pattern

obj_new[["percent.mt"]] <- PercentageFeatureSet(obj_new, pattern = mt_pattern)
obj_new[["percent.ribo"]] <- PercentageFeatureSet(obj_new, pattern = ribo_pattern)
obj_new[["log10_nCount"]] <- log10(obj_new$nCount_RNA + 1)
obj_new[["complexity"]] <- log10(obj_new$nFeature_RNA) / log10(obj_new$nCount_RNA)

# Re-attach cell-calling provenance (the rebuild from corrected counts drops
# metadata) so cell_source / emptydrops_FDR flow to steps 03-05.
obj_new <- attach_cell_calling_meta(obj_new, cfg)

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
if ("decontX_contamination" %in% colnames(obj_new@meta.data))
  qc_df$decontX_contamination <- round(obj_new$decontX_contamination, 4)
write.csv(qc_df, qc_csv_path, row.names = FALSE)

# QC metric distributions after correction (compare with step 01's quantiles)
for (sid in names(per_sample_results)) {
  per_sample_results[[sid]]$quantiles <- qc_quantiles(qc_df[qc_df$sample_id == sid, , drop = FALSE])
}

# --- JSON summary ---
successful_samples <- names(per_sample_results)[
  sapply(per_sample_results, function(x) identical(x$status, "success"))
]

summary_data <- c(
  summary_header("02_ambient", cfg, packages = if (method == "decontx") "celda" else "SoupX"),
  list(
    inputs = list(rds = rel_out(raw_rds, cfg), cell_set = cfg$cell_calling$cells,
                  matrices = sample_inputs(cfg)),
    params = c(list(method = method, round_to_int = round_to_int, skip_samples = I(skip_samples)),
               setNames(list(amb_cfg[[method]]), method)),
    n_samples = length(per_sample_results),
    n_samples_corrected = length(successful_samples),
    per_sample = per_sample_results,
    merged = list(
      n_cells = ncol(obj_new),
      n_features = nrow(obj_new),
      quantiles = qc_quantiles(qc_df),
      # Step 03 filters these corrected counts, so this preview is the one to
      # tune against (step 01's is before correction).
      threshold_preview_pct = threshold_preview(qc_df, cfg$qc$thresholds)
    ),
    outputs = list(
      rds = rel_out(rds_path, cfg),
      csv = rel_out(qc_csv_path, cfg)
    ),
    status = "success"   # per-sample failures are recorded as warnings
  )
)

write_json_summary(summary_data, file.path(outdir, "ambient_summary.json"))

# --- Plots ---
if (args$make_plots) {
  log_msg("Generating ambient correction plots")

  pre_csv <- file.path(outdir, "qc_metrics.csv")
  if (file.exists(pre_csv)) {
    pre_df <- read.csv(pre_csv, stringsAsFactors = FALSE)
    pre_df$stage <- "Before correction"
    post_df <- qc_df
    post_df$stage <- "After correction"
    comb <- rbind(pre_df[, c("sample_id", "nCount_RNA", "nFeature_RNA", "stage")],
                  post_df[, c("sample_id", "nCount_RNA", "nFeature_RNA", "stage")])
    comb$stage <- factor(comb$stage, levels = c("Before correction", "After correction"))

    p <- ggplot(comb, aes(x = sample_id, y = nCount_RNA, fill = stage)) +
      geom_violin(alpha = 0.6, scale = "width") +
      scale_y_log10() +
      labs(title = sprintf("UMI counts before/after ambient correction (%s)", method), y = "nCount_RNA", x = NULL, fill = NULL) +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    save_plot(p, "02_ambient_nCount_before_after.png", cfg, width = 8, height = 5)

    p <- ggplot(comb, aes(x = sample_id, y = nFeature_RNA, fill = stage)) +
      geom_violin(alpha = 0.6, scale = "width") +
      scale_y_log10() +
      labs(title = sprintf("Genes detected before/after ambient correction (%s)", method), y = "nFeature_RNA", x = NULL, fill = NULL) +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    save_plot(p, "02_ambient_nFeature_before_after.png", cfg, width = 8, height = 5)
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
      labs(title = if (method == "decontx") "Ambient RNA contamination (DecontX, median per cell)"
                   else "Estimated ambient RNA contamination (SoupX)",
           y = "Contamination fraction", x = NULL) +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ylim(0, max(contam_df$contamination) * 1.3)
    save_plot(p, "02_ambient_contamination.png", cfg, width = 6, height = 5)
  }

  # DecontX: per-cell contamination, and how it tracks genes detected per cell
  # (low-RNA cells are expected to carry a larger ambient share).
  if ("decontX_contamination" %in% colnames(qc_df)) {
    p <- ggplot(qc_df, aes(x = nFeature_RNA, y = decontX_contamination)) +
      geom_point(size = 0.2, alpha = 0.3) +
      scale_x_log10() +
      facet_wrap(~sample_id) +
      labs(title = "DecontX per-cell contamination vs genes detected",
           x = "Genes detected (after correction)", y = "Contamination fraction") +
      theme_bw(base_size = 11)
    save_plot(p, "02_decontx_contamination_vs_nFeature.png", cfg, width = 10, height = 8)
  }
}

log_msg("Step 02 complete")
