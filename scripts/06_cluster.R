#!/usr/bin/env Rscript
# 06_cluster.R
# Integrate (optional), embed, cluster and find cluster markers for the cells of
# the 05_normalize run named by `input`: Harmony / RPCA / CCA on its PCA, or scVI
# on raw counts; neighbors, UMAP, clustering at every resolution; FindAllMarkers
# at the default resolution (07_annotate runs read them).
# RPCA/CCA integrate per-batch SCT layers: the 05 run needs normalization.split_by
# set to integration.batch_var.
# Writes (run folder): seurat_clustered.rds, umap_coordinates.csv, cluster_markers.csv,
#         cluster_summary.json

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})
.utils_dir <- local({
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else "scripts"
})
source(file.path(.utils_dir, "utils.R"))
options(seurat_pipeline_dir = .utils_dir)  # lets run_scvi_integration find scripts/scvi/ by default

args <- parse_args()
cfg <- load_config(args$config_path, step = "06_cluster", label = args$label)
set.seed(cfg$project$seed)

options(future.globals.maxSize = 4 * 1024^3)  # 4 GiB

outdir <- run_dir(cfg)
norm_dir <- ancestor_dir(cfg, "05_normalize")
input_rds <- file.path(norm_dir, "seurat_normalized.rds")
obj <- read_ancestor_rds(cfg, "05_normalize", "seurat_normalized.rds")
log_msg(sprintf("Input: %d cells x %d genes", ncol(obj), nrow(obj)))

norm_cfg <- cfg$normalization   # the 05 run's settings
red_cfg <- cfg$reduction
clust_cfg <- cfg$clustering
int_cfg <- cfg$integration
batch_var <- int_cfg$batch_var

# --- Integration ---
integration_info <- list(enabled = FALSE)
reduction_to_use <- "pca"
if (isTRUE(int_cfg$enabled) && !(batch_var %in% colnames(obj@meta.data)))
  emit_error("integration.batch_var is not a metadata column", list(batch_var = batch_var))
# Integration needs more than one level of the batch variable.
n_batches <- if (batch_var %in% colnames(obj@meta.data)) length(unique(obj@meta.data[[batch_var]])) else 1
do_integrate <- isTRUE(int_cfg$enabled) && n_batches > 1
if (isTRUE(int_cfg$enabled) && !do_integrate)
  warn_msg(sprintf("integration.enabled is true but %s has a single level; integration skipped", batch_var))
int_method <- if (do_integrate) int_cfg$method else "none"

if (int_method %in% c("rpca", "cca")) {
  if (!identical(norm_cfg$split_by, batch_var))
    emit_error(sprintf("integration.method %s needs a 05_normalize run with normalization.split_by: %s", int_method, batch_var),
               list(normalize_split_by = norm_cfg$split_by %||% "null"))
  log_msg(sprintf("Integrating layers with %s", int_method))
  int_result <- tryCatch({
    list(obj = IntegrateLayers(obj, method = switch(int_method, rpca = RPCAIntegration, cca = CCAIntegration),
                               normalization.method = "SCT", orig.reduction = "pca",
                               new.reduction = "integrated.dr", verbose = FALSE), success = TRUE)
  }, error = function(e) list(success = FALSE, error = e$message))
  # A requested integration that fails at runtime -> HARD STOP.
  if (!int_result$success) emit_error(sprintf("Integration method '%s' failed", int_method), list(error = int_result$error))
  obj <- int_result$obj
  reduction_to_use <- "integrated.dr"
  integration_info <- list(enabled = TRUE, method = int_method, status = "success")
} else if (int_method == "harmony") {
  # Missing dependency -> HARD STOP (outside the tryCatch).
  if (!requireNamespace("harmony", quietly = TRUE))
    emit_error("Integration method 'harmony' requires the 'harmony' package", list(fix = "install.packages('harmony')"))
  log_msg("Integrating with Harmony")
  int_result <- tryCatch({
    library(harmony)
    obj <- RunHarmony(obj, group.by.vars = batch_var, reduction = "pca", reduction.save = "harmony", verbose = FALSE)
    list(success = TRUE)
  }, error = function(e) list(success = FALSE, error = e$message))
  if (!int_result$success) emit_error("Integration method 'harmony' failed", list(error = int_result$error))
  reduction_to_use <- "harmony"
  integration_info <- list(enabled = TRUE, method = "harmony", status = "success")
} else if (int_method == "scvi") {
  log_msg("Integrating with scVI (external count-based model)")
  scvi_res <- run_scvi_integration(obj, batch_var, int_cfg$scvi %||% list(), cfg)
  if (!isTRUE(scvi_res$success)) emit_error("Integration method 'scvi' failed", list(error = scvi_res$error))
  obj[["scvi"]] <- CreateDimReducObject(embeddings = scvi_res$embedding, key = "SCVI_", assay = DefaultAssay(obj))
  reduction_to_use <- "scvi"
  integration_info <- list(enabled = TRUE, method = "scvi", status = "success", scvi = scvi_res$summary,
                           cache_hit = isTRUE(scvi_res$cache_hit))
  log_msg(sprintf("scVI integration successful%s", if (isTRUE(scvi_res$cache_hit)) " (cached)" else ""))
}
log_msg(sprintf("Using reduction: %s", reduction_to_use))

pca_df <- read.csv(file.path(norm_dir, "pca_variance.csv"))
pca_cumvar <- pca_df$cumulative_variance_pct

# --- UMAP ---
# Cap by the dims actually in the reduction used (scvi has n_latent dims, which may
# be fewer than npcs), not by the PCA length.
dims_use <- min(red_cfg$dims_use, ncol(Embeddings(obj, reduction_to_use)))
log_msg(sprintf("Running UMAP with %d dims", dims_use))

obj <- RunUMAP(obj, dims = 1:dims_use, reduction = reduction_to_use,
               n.neighbors = red_cfg$umap$n_neighbors,
               min.dist = red_cfg$umap$min_dist,
               seed.use = cfg$project$seed, verbose = FALSE)

# --- Clustering ---
log_msg("Finding neighbors and clustering")
obj <- FindNeighbors(obj, dims = 1:dims_use, reduction = reduction_to_use, verbose = FALSE)

resolutions <- as.numeric(unlist(clust_cfg$resolutions))
default_res <- clust_cfg$default_resolution   # validated to be in resolutions

clustering_results <- list()
for (r in resolutions) {
  obj <- FindClusters(obj, resolution = r,
                      algorithm = c(louvain = 1, leiden = 4)[[clust_cfg$algorithm]],
                      random.seed = cfg$project$seed, verbose = FALSE)

  n_clust <- length(levels(Idents(obj)))
  # Key like the config writes it (1.0, not "1")
  res_key <- if (r == round(r)) sprintf("%.1f", r) else as.character(r)
  clustering_results[[res_key]] <- list(n_clusters = n_clust)
  log_msg(sprintf("  Resolution %.1f: %d clusters", r, n_clust))
}

# Set default resolution — find the actual column name
possible_cols <- grep(paste0("_snn_res\\.", default_res, "$"), colnames(obj@meta.data), value = TRUE)
if (length(possible_cols) > 0) {
  Idents(obj) <- possible_cols[1]
  log_msg(sprintf("Set default identity to %s", possible_cols[1]))
} else {
  emit_error("Clustering column for default_resolution not found",
             list(default_resolution = default_res, columns = grep("_snn_res\\.", colnames(obj@meta.data), value = TRUE)))
}

# --- Cells per cluster at default resolution ---
cluster_table <- table(Idents(obj))
cells_per_cluster <- as.list(setNames(as.integer(cluster_table), names(cluster_table)))

# --- Per-sample cluster composition ---
composition <- table(obj$sample_id, Idents(obj))
composition_list <- list()
for (sid in rownames(composition)) {
  composition_list[[sid]] <- as.list(setNames(as.integer(composition[sid, ]), colnames(composition)))
}

# --- Join RNA layers (07_annotate and marker finding need one layer) ---
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

# --- Cluster markers at the default resolution (RNA, log-normalized) ---
has_rna_data <- tryCatch(ncol(LayerData(obj, assay = "RNA", layer = "data")) > 0, error = function(e) FALSE)
if (!has_rna_data) obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)
mk <- cfg$markers
log_msg("Finding cluster markers")
markers <- tryCatch(
  FindAllMarkers(obj, assay = "RNA", only.pos = isTRUE(mk$only_pos), min.pct = mk$min_pct,
                 logfc.threshold = mk$logfc_threshold, test.use = mk$test, verbose = FALSE),
  error = function(e) { warn_msg(sprintf("FindAllMarkers failed: %s", e$message)); NULL })
markers_csv_path <- file.path(outdir, "cluster_markers.csv")
if (!is.null(markers) && nrow(markers) > 0) {
  write.csv(markers, markers_csv_path, row.names = FALSE)
  log_msg(sprintf("Found %d markers across %d clusters", nrow(markers), length(unique(markers$cluster))))
} else {
  warn_msg("No cluster markers found")
}

# --- Save ---
rds_path <- save_run_rds(obj, cfg, "seurat_clustered.rds")

# --- UMAP coordinates CSV (for AI/plotting without R) ---
umap_coords <- Embeddings(obj, "umap")
umap_df <- data.frame(
  barcode = rownames(umap_coords),
  UMAP_1 = round(umap_coords[, 1], 4),
  UMAP_2 = round(umap_coords[, 2], 4),
  sample_id = obj$sample_id,
  cluster = as.character(Idents(obj)),
  stringsAsFactors = FALSE
)
umap_csv_path <- file.path(outdir, "umap_coordinates.csv")
write.csv(umap_df, umap_csv_path, row.names = FALSE)

# --- JSON summary ---
integration_summary <- c(
  list(requested = isTRUE(int_cfg$enabled), method = int_cfg$method,
       batch_var = batch_var, n_batches = n_batches),
  if (do_integrate) list(status = "applied")
  else if (isTRUE(int_cfg$enabled)) list(status = "skipped", reason = sprintf("%s has a single level", batch_var))
  else list(status = "disabled", reason = "integration.enabled is false"),
  integration_info[setdiff(names(integration_info), c("enabled", "method", "status"))]
)

# Per cluster: size and the sample contributing most, against that sample's
# share of all cells (a cluster far above expectation is donor-dominated).
sample_share <- prop.table(table(obj$sample_id))
per_cluster <- lapply(setNames(nm = names(cells_per_cluster)), function(cl) {
  counts <- setNames(as.integer(composition[, cl]), rownames(composition))
  top <- names(which.max(counts))
  list(n_cells = cells_per_cluster[[cl]], top_sample = top,
       top_sample_pct = round(100 * max(counts) / sum(counts), 1),
       top_sample_expected_pct = round(100 * unname(sample_share[top]), 1))
})

summary_data <- c(
  summary_header("06_cluster", cfg, packages = if (int_method == "harmony") "harmony"),
  list(
    inputs = list(rds = rel_out(input_rds, cfg)),
    params = list(integration = int_cfg, reduction = red_cfg[c("dims_use", "umap")], clustering = clust_cfg,
                  markers = cfg$markers),
    n_cells = ncol(obj),
    n_features = nrow(obj),               # default assay (SCT or RNA)
    n_features_rna = nrow(obj[["RNA"]]),
    integration = integration_summary,
    reduction_used = reduction_to_use,    # embedding for neighbors, UMAP, clustering
    pca = list(
      npcs_computed = length(pca_cumvar),
      cumulative_variance_at_dims_used = round(pca_cumvar[min(dims_use, length(pca_cumvar))], 2)
    ),
    umap = list(
      dims_used = dims_use,
      n_neighbors = red_cfg$umap$n_neighbors,
      min_dist = red_cfg$umap$min_dist
    ),
    clustering = list(
      algorithm = clust_cfg$algorithm,
      resolutions = clustering_results,
      default_resolution = default_res,
      default_n_clusters = length(cluster_table)
    ),
    per_cluster = per_cluster,
    cluster_counts_by_sample = composition_list,
    cluster_markers = list(n_markers = if (is.null(markers)) 0L else nrow(markers)),
    outputs = list(
      rds = rel_out(rds_path, cfg),
      umap_csv = rel_out(umap_csv_path, cfg),
      markers_csv = if (file.exists(markers_csv_path)) rel_out(markers_csv_path, cfg) else NULL
    ),
    status = "success"
  )
)

write_json_summary(summary_data, file.path(outdir, STEPS[["06_cluster"]]$summary))

# --- Plots ---
if (args$make_plots) {
  log_msg("Generating clustering plots")

  # Elbow plot (the 05 run's PCA) with this run's dims_use
  p <- ggplot(pca_df, aes(x = PC, y = variance_pct)) +
    geom_point() + geom_line() +
    geom_vline(xintercept = dims_use, linetype = "dashed", color = "red") +
    labs(title = "PCA elbow plot (dashed: dims used)", x = "PC", y = "% of total variance (variable features)") +
    theme_bw(base_size = 12)
  save_plot(p, "06_elbow_plot.png", cfg, width = 7, height = 5)

  # UMAP by cluster
  p <- DimPlot(obj, reduction = "umap", label = TRUE, repel = TRUE) +
    ggtitle("UMAP — clusters") + NoLegend()
  save_plot(p, "06_umap_clusters.png", cfg, width = 8, height = 7)

  # UMAP by sample
  p <- DimPlot(obj, reduction = "umap", group.by = "sample_id") +
    ggtitle("UMAP — samples")
  save_plot(p, "06_umap_samples.png", cfg, width = 9, height = 7)

  # UMAP by integration batch variable (e.g. chemistry) — the effect integration
  # is correcting. Skipped when batch_var is sample_id (already plotted above).
  if (do_integrate && batch_var %in% colnames(obj@meta.data) && batch_var != "sample_id") {
    p <- DimPlot(obj, reduction = "umap", group.by = batch_var) +
      ggtitle(sprintf("UMAP — %s (integration batch)", batch_var))
    save_plot(p, "06_umap_batch.png", cfg, width = 9, height = 7)
  }

  # UMAP split by sample: grid and figure size scale with the number of samples so
  # panels stay roughly square (~3.3 in) instead of being squeezed into a fixed height.
  n_samples <- length(unique(obj$sample_id))
  split_ncol <- min(4L, n_samples)
  split_nrow <- ceiling(n_samples / split_ncol)
  p <- DimPlot(obj, reduction = "umap", split.by = "sample_id", ncol = split_ncol) +
    ggtitle("UMAP split by sample") + NoLegend()
  save_plot(p, "06_umap_split_sample.png", cfg,
            width = 3.3 * split_ncol + 0.5, height = 3.3 * split_nrow + 0.6)

  # QC metrics on the UMAP and per cluster: is any cluster driven by quality
  # (depth, complexity, MT) rather than biology?
  qc_feats <- intersect(c("log10_nCount", "nFeature_RNA", "percent.mt"), colnames(obj@meta.data))
  if (length(qc_feats)) {
    p <- FeaturePlot(obj, features = qc_feats, reduction = "umap", ncol = length(qc_feats),
                     order = TRUE, raster = FALSE)
    save_plot(p, "06_umap_qc_metrics.png", cfg, width = 5 * length(qc_feats), height = 4.8)

    clust_levels <- levels(Idents(obj))
    qc_long <- do.call(rbind, lapply(qc_feats, function(f) data.frame(
      cluster = factor(as.character(Idents(obj)), levels = clust_levels),
      metric = f, value = obj@meta.data[[f]])))
    p <- ggplot(qc_long, aes(x = cluster, y = value, fill = cluster)) +
      geom_violin(scale = "width", linewidth = 0.2) +
      facet_wrap(~metric, ncol = 1, scales = "free_y") +
      labs(title = "QC metrics per cluster", x = "Cluster", y = NULL) +
      theme_bw(base_size = 11) + NoLegend()
    save_plot(p, "06_violin_qc_per_cluster.png", cfg,
              width = max(8, 0.35 * length(clust_levels) + 2), height = 8)
  }

  # Donor x cluster heatmap: % of each sample's own cells in each cluster, so
  # differences in sample size don't masquerade as donor-specific clusters.
  if (n_samples > 1) {
    donor_pct <- 100 * prop.table(table(obj$sample_id, Idents(obj)), 1)
    dp <- as.data.frame(donor_pct)
    colnames(dp) <- c("sample", "cluster", "pct")
    dp$cluster <- factor(dp$cluster, levels = levels(Idents(obj)))
    p <- ggplot(dp, aes(x = cluster, y = sample, fill = pct)) +
      geom_tile(color = "grey90") +
      geom_text(aes(label = ifelse(pct >= 1, sprintf("%.0f", pct), "")), size = 2.6) +
      scale_fill_gradient(low = "white", high = "firebrick", name = "% of\nsample") +
      labs(title = "% of each sample's cells in each cluster", x = "Cluster", y = NULL) +
      theme_minimal(base_size = 11)
    save_plot(p, "06_heatmap_sample_cluster.png", cfg,
              width = max(8, 0.35 * nlevels(dp$cluster) + 3), height = 1 + 0.4 * n_samples)
  }

  # Cell-calling provenance: where do emptyDrops-recovered cells land?
  if ("cell_source" %in% colnames(obj@meta.data)) {
    p <- DimPlot(obj, reduction = "umap", group.by = "cell_source", order = TRUE) +
      ggtitle("UMAP — cell source (emptyDrops-recovered vs filtered)")
    save_plot(p, "06_umap_cell_source.png", cfg, width = 9, height = 7)

    if ("emptydrops_FDR" %in% colnames(obj@meta.data)) {
      obj$neg_log10_FDR <- -log10(pmax(obj$emptydrops_FDR, 1e-6))
      p <- FeaturePlot(obj, features = "neg_log10_FDR", reduction = "umap") +
        ggtitle("UMAP — emptyDrops -log10(FDR)  (grey = not tested, total <= lower)")
      save_plot(p, "06_umap_emptydrops_fdr.png", cfg, width = 9, height = 7)
    }
  }

  # Cluster sizes bar chart
  cluster_df <- data.frame(
    cluster = names(cells_per_cluster),
    cells = unlist(cells_per_cluster),
    stringsAsFactors = FALSE
  )
  cluster_df$cluster <- factor(cluster_df$cluster, levels = cluster_df$cluster[order(cluster_df$cells, decreasing = TRUE)])

  p <- ggplot(cluster_df, aes(x = cluster, y = cells)) +
    geom_col(fill = "steelblue") +
    geom_text(aes(label = cells), vjust = -0.3, size = 3) +
    labs(title = "Cells per cluster", x = "Cluster", y = "Cells") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p, "06_cluster_sizes.png", cfg, width = 8, height = 5)
}

log_msg("Step 06 complete")
