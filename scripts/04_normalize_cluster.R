#!/usr/bin/env Rscript
# 04_normalize_cluster.R
# Normalize, integrate (optional), PCA, UMAP, cluster at multiple resolutions.
# Writes: seurat_clustered.rds, pca_variance.csv, umap_coordinates.csv, cluster_summary.json

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})
source("scripts/utils.R")

args <- parse_args()
cfg <- load_config(args$config_path)
set.seed(cfg$project$seed)

options(future.globals.maxSize = 4 * 1024^3)  # 4 GiB

outdir <- output_dir(cfg)

# --- Load input ---
input_rds <- file.path(outdir, "seurat_filtered.rds")
if (!file.exists(input_rds)) {
  emit_error("Input RDS not found. Run step 03 first.", list(expected = input_rds))
}

log_msg("Loading seurat_filtered.rds")
obj <- readRDS(input_rds)
log_msg(sprintf("Input: %d cells x %d genes", ncol(obj), nrow(obj)))

norm_cfg <- cfg$normalization
red_cfg <- cfg$reduction
clust_cfg <- cfg$clustering
int_cfg <- cfg$integration %||% list(enabled = FALSE)

seurat_v5 <- packageVersion("Seurat") >= "5.0.0"
log_msg(sprintf("Seurat version: %s (v5 API: %s)", packageVersion("Seurat"), seurat_v5))

# --- Normalization + Integration ---
integration_info <- list(enabled = FALSE)
norm_method <- norm_cfg$method %||% "SCTransform"
reduction_to_use <- "pca"

n_samples <- length(unique(obj$sample_id))
do_integrate <- isTRUE(int_cfg$enabled) && n_samples > 1

if (do_integrate) {
  log_msg(sprintf("Integration enabled (%s), %d samples", int_cfg$method, n_samples))
}

if (norm_method == "SCTransform") {
  log_msg("Normalizing with SCTransform")

  vst_flavor <- norm_cfg$sct_vst_flavor %||% "v2"
  n_features <- norm_cfg$sct_variable_features %||% 3000
  vars_regress <- norm_cfg$sct_vars_to_regress

  obj <- SCTransform(obj, vars.to.regress = vars_regress,
                     variable.features.n = n_features,
                     vst.flavor = vst_flavor, verbose = FALSE,
                     seed.use = cfg$project$seed)

  obj <- RunPCA(obj, npcs = red_cfg$npcs, verbose = FALSE, seed.use = cfg$project$seed)

  if (do_integrate) {
    int_method <- int_cfg$method %||% "harmony"

    if (int_method == "harmony") {
      log_msg("Integrating with Harmony")

      int_result <- tryCatch({
        if (!requireNamespace("harmony", quietly = TRUE)) {
          stop("harmony package not installed")
        }
        library(harmony)
        obj <- RunHarmony(obj, group.by.vars = "sample_id", reduction = "pca",
                          reduction.save = "harmony", verbose = FALSE)
        list(success = TRUE)
      }, error = function(e) {
        log_msg(sprintf("WARNING: Harmony failed: %s. Proceeding without integration.", e$message))
        list(success = FALSE, error = e$message)
      })

      if (int_result$success) {
        reduction_to_use <- "harmony"
        integration_info <- list(enabled = TRUE, method = "harmony", status = "success")
        log_msg("Harmony integration successful")
      } else {
        reduction_to_use <- "pca"
        integration_info <- list(enabled = TRUE, method = "harmony", status = "failed",
                                 error = int_result$error)
      }

    } else if (seurat_v5) {
      log_msg(sprintf("Integrating layers with %s", int_method))

      obj[["RNA"]] <- split(obj[["RNA"]], f = obj$sample_id)
      obj[["SCT"]] <- split(obj[["SCT"]], f = obj$sample_id)
      obj <- RunPCA(obj, npcs = red_cfg$npcs, verbose = FALSE, seed.use = cfg$project$seed)

      int_result <- tryCatch({
        integrated_obj <- IntegrateLayers(
          obj,
          method = switch(int_method,
            "rpca" = RPCAIntegration,
            "cca" = CCAIntegration
          ),
          orig.reduction = "pca",
          new.reduction = "integrated.pca",
          verbose = FALSE
        )
        list(obj = integrated_obj, success = TRUE)
      }, error = function(e) {
        log_msg(sprintf("WARNING: Integration failed: %s. Proceeding without integration.", e$message))
        list(obj = NULL, success = FALSE, error = e$message)
      })

      if (int_result$success) {
        obj <- int_result$obj
        reduction_to_use <- "integrated.pca"
        integration_info <- list(enabled = TRUE, method = int_method, status = "success")
        log_msg("Integration successful")
      } else {
        reduction_to_use <- "pca"
        integration_info <- list(enabled = TRUE, method = int_method, status = "failed",
                                 error = int_result$error)
      }

    } else {
      # Seurat v4: split, SCT per sample, find anchors, integrate
      log_msg("Using Seurat v4 integration workflow")
      obj_list <- SplitObject(obj, split.by = "sample_id")

      obj_list <- lapply(obj_list, function(x) {
        SCTransform(x, vars.to.regress = vars_regress,
                    variable.features.n = n_features,
                    vst.flavor = vst_flavor, verbose = FALSE,
                    seed.use = cfg$project$seed)
      })

      features <- SelectIntegrationFeatures(obj_list, nfeatures = n_features)
      obj_list <- PrepSCTIntegration(obj_list, anchor.features = features)
      obj_list <- lapply(obj_list, RunPCA, features = features, verbose = FALSE, seed.use = cfg$project$seed)

      int_method <- int_cfg$method %||% "rpca"
      anchors <- tryCatch(
        FindIntegrationAnchors(obj_list, normalization.method = "SCT",
                               anchor.features = features,
                               reduction = int_method, dims = 1:red_cfg$npcs),
        error = function(e) {
          log_msg(sprintf("WARNING: FindIntegrationAnchors failed: %s", e$message))
          NULL
        }
      )

      if (!is.null(anchors)) {
        obj <- IntegrateData(anchors, normalization.method = "SCT", dims = 1:red_cfg$npcs)
        DefaultAssay(obj) <- "integrated"
        obj <- RunPCA(obj, npcs = red_cfg$npcs, verbose = FALSE, seed.use = cfg$project$seed)
        integration_info <- list(enabled = TRUE, method = int_method, status = "success")
      } else {
        integration_info <- list(enabled = TRUE, method = int_method, status = "failed",
                                 error = "anchor finding failed")
      }
    }
  }  # end do_integrate

} else {
  # LogNormalize
  log_msg("Normalizing with LogNormalize")
  obj <- NormalizeData(obj, scale.factor = norm_cfg$scale_factor %||% 10000, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = norm_cfg$log_variable_features %||% 2000, verbose = FALSE)
  obj <- ScaleData(obj, vars.to.regress = norm_cfg$sct_vars_to_regress, verbose = FALSE)
  obj <- RunPCA(obj, npcs = red_cfg$npcs, verbose = FALSE, seed.use = cfg$project$seed)

  if (do_integrate) {
    log_msg("WARNING: Integration with LogNormalize not fully implemented, proceeding without")
    integration_info <- list(enabled = TRUE, status = "skipped", reason = "LogNormalize integration not implemented")
  }
}

log_msg(sprintf("Using reduction: %s", reduction_to_use))

# --- PCA variance ---
pca_stdev <- obj[[reduction_to_use]]@stdev
pca_var <- pca_stdev^2
pca_var_pct <- round(100 * pca_var / sum(pca_var), 3)
pca_cumvar <- round(cumsum(pca_var_pct), 3)

pca_df <- data.frame(
  PC = seq_along(pca_stdev),
  stdev = round(pca_stdev, 4),
  variance_pct = pca_var_pct,
  cumulative_variance_pct = pca_cumvar
)
pca_csv_path <- file.path(outdir, "pca_variance.csv")
write.csv(pca_df, pca_csv_path, row.names = FALSE)
log_msg(sprintf("PCA: cumulative variance at PC%d = %.1f%%", red_cfg$dims_use, pca_cumvar[red_cfg$dims_use]))

# --- Elbow heuristic ---
suggested_dims <- which(pca_var_pct < 0.5)[1]
if (is.na(suggested_dims)) suggested_dims <- red_cfg$npcs
suggested_dims <- max(suggested_dims - 1, 5)

# --- UMAP ---
dims_use <- min(red_cfg$dims_use, length(pca_stdev))
log_msg(sprintf("Running UMAP with %d dims", dims_use))

obj <- RunUMAP(obj, dims = 1:dims_use, reduction = reduction_to_use,
               n.neighbors = red_cfg$umap_n_neighbors %||% 30,
               min.dist = red_cfg$umap_min_dist %||% 0.3,
               seed.use = cfg$project$seed, verbose = FALSE)

# --- Clustering ---
log_msg("Finding neighbors and clustering")
obj <- FindNeighbors(obj, dims = 1:dims_use, reduction = reduction_to_use, verbose = FALSE)

resolutions <- clust_cfg$resolutions %||% c(0.5)
default_res <- clust_cfg$default_resolution %||% 0.5

clustering_results <- list()
for (r in resolutions) {
  obj <- FindClusters(obj, resolution = r,
                      algorithm = switch(clust_cfg$algorithm %||% "louvain",
                                         "louvain" = 1, "leiden" = 4, 1),
                      random.seed = cfg$project$seed, verbose = FALSE)

  n_clust <- length(levels(Idents(obj)))
  clustering_results[[as.character(r)]] <- list(n_clusters = n_clust)
  log_msg(sprintf("  Resolution %.1f: %d clusters", r, n_clust))
}

# Set default resolution — find the actual column name
possible_cols <- grep(paste0("_snn_res\\.", default_res, "$"), colnames(obj@meta.data), value = TRUE)
if (length(possible_cols) > 0) {
  Idents(obj) <- possible_cols[1]
  log_msg(sprintf("Set default identity to %s", possible_cols[1]))
} else {
  log_msg("Using last clustering resolution as default identity")
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

# --- Join RNA layers for downstream use (step 05 needs this) ---
if (seurat_v5) {
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
  log_msg("Joined RNA assay layers for downstream compatibility")
}

# --- Save ---
rds_path <- file.path(outdir, "seurat_clustered.rds")
log_msg(sprintf("Saving clustered object to %s", rds_path))
saveRDS(obj, rds_path)

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
summary_data <- c(
  summary_header("04_normalize_cluster", cfg),
  list(
    cells = ncol(obj),
    genes = nrow(obj),
    normalization = list(
      method = norm_method,
      variable_features = norm_cfg$sct_variable_features %||% norm_cfg$log_variable_features,
      vars_regressed = norm_cfg$sct_vars_to_regress
    ),
    integration = integration_info,
    pca = list(
      npcs_computed = length(pca_stdev),
      dims_used = dims_use,
      cumulative_variance_at_dims_used = pca_cumvar[dims_use],
      suggested_dims = suggested_dims
    ),
    umap = list(
      dims = dims_use,
      n_neighbors = red_cfg$umap_n_neighbors %||% 30,
      min_dist = red_cfg$umap_min_dist %||% 0.3
    ),
    clustering = list(
      algorithm = clust_cfg$algorithm %||% "louvain",
      resolutions = clustering_results,
      default_resolution = default_res,
      default_n_clusters = length(cluster_table)
    ),
    cells_per_cluster = cells_per_cluster,
    per_sample_composition = composition_list,
    outputs = list(
      rds = rds_path,
      pca_csv = pca_csv_path,
      umap_csv = umap_csv_path
    ),
    status = "success"
  )
)

write_json_summary(summary_data, file.path(outdir, "cluster_summary.json"))

# --- Plots ---
if (args$make_plots) {
  log_msg("Generating clustering plots")

  # Elbow plot
  p <- ggplot(pca_df, aes(x = PC, y = variance_pct)) +
    geom_point() + geom_line() +
    geom_vline(xintercept = dims_use, linetype = "dashed", color = "red") +
    labs(title = "PCA elbow plot", x = "PC", y = "% Variance explained") +
    theme_bw(base_size = 12)
  save_plot(p, "04_elbow_plot.png", cfg, width = 7, height = 5)

  # UMAP by cluster
  p <- DimPlot(obj, reduction = "umap", label = TRUE, repel = TRUE) +
    ggtitle("UMAP — clusters") + NoLegend()
  save_plot(p, "04_umap_clusters.png", cfg, width = 8, height = 7)

  # UMAP by sample
  p <- DimPlot(obj, reduction = "umap", group.by = "sample_id") +
    ggtitle("UMAP — samples")
  save_plot(p, "04_umap_samples.png", cfg, width = 9, height = 7)

  # UMAP split by sample
  p <- DimPlot(obj, reduction = "umap", split.by = "sample_id", ncol = 2) +
    ggtitle("UMAP split by sample") + NoLegend()
  save_plot(p, "04_umap_split_sample.png", cfg, width = 10, height = 10)

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
  save_plot(p, "04_cluster_sizes.png", cfg, width = 8, height = 5)
}

log_msg("Step 04 complete")
