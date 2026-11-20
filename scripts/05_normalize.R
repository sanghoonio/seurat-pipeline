#!/usr/bin/env Rscript
# 05_normalize.R
# Normalize the cells of the 04_filter run named by `input` (SCTransform or
# LogNormalize), select variable features, and run PCA. Integration, UMAP and
# clustering are 06_cluster's, so their variants reuse one normalization.
#
# normalization.split_by: SCTransform each level of a metadata column separately
# (per-layer models). Needed by 06_cluster runs with integration.method rpca/cca,
# which integrate those layers; harmony, scVI and no integration use one joint model.
#
# Writes (run folder): seurat_normalized.rds, pca_variance.csv, normalize_summary.json

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
cfg <- load_config(args$config_path, step = "05_normalize", label = args$label)
set.seed(cfg$project$seed)
options(future.globals.maxSize = 4 * 1024^3)  # 4 GiB

outdir <- run_dir(cfg)
norm_cfg <- cfg$normalization
npcs <- cfg$reduction$npcs
in_rds <- file.path(ancestor_dir(cfg, "04_filter"), "seurat_filtered.rds")
obj <- read_ancestor_rds(cfg, "04_filter", "seurat_filtered.rds")
log_msg(sprintf("Input: %d cells x %d genes", ncol(obj), nrow(obj)))

if (packageVersion("Seurat") < "5.0.0") emit_error("05_normalize needs Seurat >= 5")
split_by <- norm_cfg$split_by
if (!is.null(split_by)) {
  if (norm_cfg$method != "SCTransform")
    emit_error("normalization.split_by needs normalization.method: SCTransform")
  if (!split_by %in% colnames(obj@meta.data))
    emit_error("normalization.split_by is not a metadata column", list(split_by = split_by))
}

if (norm_cfg$method == "SCTransform") {
  vars_regress <- unlist(norm_cfg$vars_to_regress)   # [] -> NULL
  if (!is.null(split_by)) {
    log_msg(sprintf("SCTransform per level of %s (%d levels)", split_by, length(unique(obj@meta.data[[split_by]]))))
    obj[["RNA"]] <- split(obj[["RNA"]], f = obj@meta.data[[split_by]])
  } else {
    log_msg("SCTransform (one joint model)")
  }
  obj <- SCTransform(obj, vars.to.regress = vars_regress,
                     variable.features.n = norm_cfg$sctransform$variable_features,
                     vst.flavor = norm_cfg$sctransform$vst_flavor, verbose = FALSE,
                     seed.use = cfg$project$seed)
} else {
  log_msg("LogNormalize")
  obj <- NormalizeData(obj, scale.factor = norm_cfg$lognormalize$scale_factor, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = norm_cfg$lognormalize$variable_features, verbose = FALSE)
  obj <- ScaleData(obj, vars.to.regress = unlist(norm_cfg$vars_to_regress), verbose = FALSE)
}
obj <- RunPCA(obj, npcs = npcs, verbose = FALSE, seed.use = cfg$project$seed)

# --- PCA variance ---
# Shares are of the TOTAL variance of the matrix PCA ran on (scaled data for the
# variable features), so they don't depend on how many PCs were computed.
pca_stdev <- obj[["pca"]]@stdev
total_var <- Misc(obj[["pca"]], "total.variance")
if (is.null(total_var)) {
  sd_mat <- GetAssayData(obj, assay = DefaultAssay(obj[["pca"]]), layer = "scale.data")
  feats <- intersect(rownames(Loadings(obj[["pca"]])), rownames(sd_mat))
  total_var <- sum(apply(sd_mat[feats, , drop = FALSE], 1, var))
  rm(sd_mat)
}
pca_var_pct <- 100 * pca_stdev^2 / total_var
pca_df <- data.frame(PC = seq_along(pca_stdev), stdev = round(pca_stdev, 4),
                     variance_pct = round(pca_var_pct, 3),
                     cumulative_variance_pct = round(cumsum(pca_var_pct), 3))
pca_csv_path <- file.path(outdir, "pca_variance.csv")
write.csv(pca_df, pca_csv_path, row.names = FALSE)

rds_path <- save_run_rds(obj, cfg, "seurat_normalized.rds")

summary_data <- c(
  summary_header("05_normalize", cfg, packages = if (norm_cfg$method == "SCTransform") "glmGamPoi"),
  list(
    inputs = list(rds = rel_out(in_rds, cfg)),
    params = list(normalization = norm_cfg, npcs = npcs),
    n_cells = ncol(obj),
    n_features = nrow(obj),
    normalization = list(
      method = norm_cfg$method,
      split_by = split_by,
      variable_features = length(VariableFeatures(obj)),
      vars_regressed = I(unlist(norm_cfg$vars_to_regress) %||% character(0))
    ),
    pca = list(npcs_computed = length(pca_stdev),
               total_variance_basis = "scaled data, variable features",
               cumulative_variance_pct = I(round(cumsum(pca_var_pct), 2))),
    outputs = list(rds = rel_out(rds_path, cfg), pca_csv = rel_out(pca_csv_path, cfg)),
    status = "success"
  )
)
write_json_summary(summary_data, file.path(outdir, STEPS[["05_normalize"]]$summary))

if (args$make_plots) {
  p <- ggplot(pca_df, aes(x = PC, y = variance_pct)) +
    geom_point() + geom_line() +
    labs(title = "PCA elbow plot", x = "PC", y = "% of total variance (variable features)") +
    theme_bw(base_size = 12)
  save_plot(p, "05_elbow_plot.png", cfg, width = 7, height = 5)
}

log_msg("Step 05 complete")
