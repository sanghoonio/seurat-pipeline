#!/usr/bin/env Rscript
# 07_annotate.R
# Annotate the clusters of the 06_cluster run named by `input`: Azimuth, SingleR
# (one or more celldex references), CellTypist (one or more models, via the
# bundled Python env), a custom Seurat reference (annotation.custom_reference),
# and marker-panel module scores. Each method runs only when configured and
# writes its own metadata columns and per-cell CSV under annotation/ in the run
# folder. Builds a per-cluster table with a suggested label, using the 06 run's
# cluster markers. Optional: label the cells the 04_filter run removed
# (annotation.removed_cells), split diagnostics by sample metadata
# (annotation.compare_by). With a marker panel, also reports how much of each
# marker set the 02_ambient run's correction removed.
# Writes (run folder): seurat_final.rds, cluster_annotations.csv, annotation/*.csv,
#         filtered_out_annotation.csv (removed_cells), annotation_summary.json

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})
.utils_dir <- local({
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else "scripts"
})
source(file.path(.utils_dir, "utils.R"))
options(seurat_pipeline_dir = .utils_dir)  # lets run_celltypist find the bundled env

args <- parse_args()
cfg <- load_config(args$config_path, step = "07_annotate", label = args$label)
set.seed(cfg$project$seed)

options(future.globals.maxSize = 4 * 1024^3)  # 4 GiB

outdir <- run_dir(cfg)
ann_dir <- file.path(outdir, "annotation")   # per-method, per-cell label CSVs
dir.create(ann_dir, showWarnings = FALSE)
safe_name <- function(x) gsub("[^[:alnum:]_]+", "_", x)
method_label_cols <- character(0)            # label columns added by SingleR/CellTypist/custom reference (UMAPs)

# --- Load input ---
clust_dir <- ancestor_dir(cfg, "06_cluster")
input_rds <- file.path(clust_dir, "seurat_clustered.rds")
obj <- read_ancestor_rds(cfg, "06_cluster", "seurat_clustered.rds")
log_msg(sprintf("Input: %d cells, %d clusters", ncol(obj), length(levels(Idents(obj)))))

ann_cfg <- cfg$annotation %||% list()
annotation_methods <- list()

# Determine which assay to use for marker finding and scoring.
# Prefer RNA with joined layers; fall back to SCT.
rna_assay <- "RNA"
if (rna_assay %in% Assays(obj)) {
  if (packageVersion("Seurat") >= "5.0.0") {
    tryCatch({
      obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
    }, error = function(e) {
      log_msg(sprintf("JoinLayers RNA note: %s", e$message))
    })
  }
  has_rna_data <- tryCatch({
    d <- LayerData(obj, assay = "RNA", layer = "data")
    ncol(d) > 0
  }, error = function(e) FALSE)

  if (!has_rna_data) {
    log_msg("RNA assay not normalized, running NormalizeData")
    DefaultAssay(obj) <- "RNA"
    obj <- NormalizeData(obj, verbose = FALSE)
  }
} else {
  rna_assay <- "SCT"
}

log_msg(sprintf("Using assay '%s' for marker analysis", rna_assay))
DefaultAssay(obj) <- rna_assay

# Resolved per-cell label at the primary level: references such as lungref use
# the literal "None" where a cell type has no subdivision at this level (e.g.
# smooth muscle at ann_level_3). Fill those from the nearest coarser level
# (label columns run coarse -> fine) instead of treating "None" as a cell type.
resolve_azimuth_label <- function(md, label_cols, primary_col) {
  lab <- as.character(md[[primary_col]])
  n_fallback <- 0L
  for (col in rev(head(label_cols, match(primary_col, label_cols) - 1))) {
    miss <- is.na(lab) | lab %in% c("None", "none", "")
    if (!any(miss)) break
    lab[miss] <- as.character(md[[col]][miss])
    n_fallback <- n_fallback + sum(miss & !(lab %in% c("None", "none", "")))
  }
  list(label = lab, n_fallback = n_fallback)
}

# --- Azimuth ---
azimuth_info <- list(attempted = FALSE)

if (!is.null(ann_cfg$azimuth$reference)) {
  # Reference requested but dep missing -> HARD STOP (don't silently drop annotation).
  if (!requireNamespace("Azimuth", quietly = TRUE))
    emit_error("annotation.azimuth.reference is set but the 'Azimuth' package is not installed",
               list(fix = "install Azimuth (or set annotation.azimuth.reference: null to skip)"))
  {
    log_msg(sprintf("Running Azimuth with reference: %s", ann_cfg$azimuth$reference))
    azimuth_info <- tryCatch({
      suppressPackageStartupMessages(library(Azimuth))
      obj <- RunAzimuth(obj, reference = ann_cfg$azimuth$reference)

      pred_cols <- grep("^predicted\\.", colnames(obj@meta.data), value = TRUE)
      score_cols <- grep("\\.score$", pred_cols, value = TRUE)
      label_cols <- setdiff(pred_cols, score_cols)
      label_cols <- setdiff(label_cols, "predicted.id")

      info <- list(attempted = TRUE, status = "success", reference = ann_cfg$azimuth$reference,
                   label_columns = I(label_cols), score_columns = I(score_cols))

      if (length(label_cols) > 0) {
        primary_lvl <- ann_cfg$azimuth$primary_level
        if (!is.null(primary_lvl)) {
          primary_col <- paste0("predicted.", primary_lvl)
          # Hard stop (quits even inside this tryCatch): a misnamed level would
          # otherwise silently fall back to another level's labels.
          if (!primary_col %in% label_cols)
            emit_error("annotation.azimuth.primary_level is not a level of this reference",
                       list(primary_level = primary_lvl,
                            available = sub("^predicted\\.", "", label_cols)))
        } else {
          primary_col <- label_cols[1]
          log_msg(sprintf("annotation.azimuth.primary_level not set; using %s", primary_col))
        }
        info$primary_level <- sub("^predicted\\.", "", primary_col)
        info$primary_label_column <- primary_col
        rl <- resolve_azimuth_label(obj@meta.data, label_cols, primary_col)
        lab <- rl$label
        obj$azimuth_label <- lab
        info$n_cells_label_from_coarser_level <- rl$n_fallback
        info$cell_types_found <- I(sort(unique(lab)))
        score_col <- paste0(primary_col, ".score")
        if (score_col %in% colnames(obj@meta.data)) {
          info$median_prediction_score <- round(median(obj@meta.data[[score_col]], na.rm = TRUE), 3)
          info$low_confidence_pct <- round(100 * mean(obj@meta.data[[score_col]] < 0.5, na.rm = TRUE), 1)
        }
      }
      info
    }, error = function(e) {
      warn_msg(sprintf("Azimuth failed: %s", e$message))
      list(attempted = TRUE, status = "failed", error = e$message)
    })
  }
}
annotation_methods$azimuth <- azimuth_info

# --- SingleR ---
# One run per celldex reference; each writes singler_<ref> / singler_<ref>_pruned
# and outputs/annotation/singler_<ref>.csv. The first successful reference also
# fills singler_label, the SingleR fallback in suggested_label.
singler_refs <- as.character(unlist(ann_cfg$singler$reference))
singler_level <- ann_cfg$singler$labels %||% "main"
singler_info <- list(attempted = FALSE)

if (length(singler_refs)) {
  # Reference requested but deps missing -> HARD STOP (don't silently drop annotation).
  .singler_missing <- setdiff(c("SingleR", "celldex", "SingleCellExperiment"),
                              rownames(installed.packages()))
  if (length(.singler_missing) > 0)
    emit_error("annotation.singler.reference is set but required packages are not installed",
               list(missing = .singler_missing,
                    fix = "BiocManager::install(the missing packages)  (or set annotation.singler.reference: null)"))
  suppressPackageStartupMessages({
    library(SingleR)
    library(celldex)
    library(SingleCellExperiment)
  })
  sce <- as.SingleCellExperiment(obj)
  singler_info <- list(attempted = TRUE, labels = singler_level, references = list())
  for (r in singler_refs) {
    log_msg(sprintf("Running SingleR with reference: %s (%s labels)", r, singler_level))
    singler_info$references[[r]] <- tryCatch({
      ref <- do.call(r, list())
      lab_col <- if (singler_level == "fine") "label.fine" else "label.main"
      present <- intersect(rownames(ref), rownames(sce))
      res <- SingleR(test = sce, ref = ref, labels = ref[[lab_col]])
      col <- paste0("singler_", safe_name(r))
      obj[[col]] <- res$labels
      obj[[paste0(col, "_pruned")]] <- res$pruned.labels
      method_label_cols <- c(method_label_cols, col)
      if (!"singler_label" %in% colnames(obj@meta.data)) {
        obj$singler_label <- res$labels
        obj$singler_pruned_label <- res$pruned.labels
      }
      csv <- file.path(ann_dir, paste0(col, ".csv"))
      write.csv(data.frame(barcode = colnames(obj), label = res$labels, pruned_label = res$pruned.labels,
                           delta_next = round(res$delta.next, 4)), csv, row.names = FALSE)
      list(status = "success", label_column = col,
           n_reference_genes = nrow(ref), n_reference_genes_present = length(present),
           gene_coverage_pct = round(100 * length(present) / nrow(ref), 1),
           n_pruned = sum(is.na(res$pruned.labels)),
           cell_types_found = I(sort(unique(res$labels))), csv = rel_out(csv, cfg))
    }, error = function(e) {
      warn_msg(sprintf("SingleR (%s) failed: %s", r, e$message))
      list(status = "failed", error = e$message)
    })
  }
  st <- vapply(singler_info$references, function(x) x$status, "")
  singler_info$status <- if (all(st == "success")) "success" else if (any(st == "success")) "partial" else "failed"
  rm(sce)
}
annotation_methods$singler <- singler_info

# --- CellTypist ---
# Each model writes celltypist_<model> (majority-voting label when enabled, else
# the per-cell label), celltypist_<model>_conf, and outputs/annotation/celltypist_<model>.csv.
ct_models <- as.character(unlist(ann_cfg$celltypist$models))
celltypist_info <- list(attempted = FALSE)
if (length(ct_models)) {
  ct <- run_celltypist(obj, ct_models, ann_cfg$celltypist, ann_dir)
  if (!isTRUE(ct$success)) {
    warn_msg(sprintf("CellTypist failed: %s", ct$error))
    celltypist_info <- list(attempted = TRUE, status = "failed", error = ct$error)
  } else {
    celltypist_info <- list(attempted = TRUE, status = "success",
                            majority_voting = isTRUE(ann_cfg$celltypist$majority_voting),
                            celltypist_version = ct$summary$celltypist_version, models = list())
    for (k in names(ct$per_model)) {
      d <- ct$per_model[[k]]
      col <- paste0("celltypist_", safe_name(k))
      obj[[col]] <- if ("majority" %in% colnames(d)) d$majority else d$label
      obj[[paste0(col, "_conf")]] <- d$conf
      method_label_cols <- c(method_label_cols, col)
      celltypist_info$models[[k]] <- c(list(label_column = col, csv = rel_out(ct$csv[[k]], cfg),
                                            median_conf = round(median(d$conf, na.rm = TRUE), 3)),
                                       ct$summary$models[[k]])
    }
  }
}
annotation_methods$celltypist <- celltypist_info

# --- Custom Seurat reference (e.g. LungMAP CellRef; built by scripts/refs/) ---
# Transfers each of label_columns: customref_<name>_<col> and _score, and
# outputs/annotation/customref_<name>.csv (<name> = reference file name).
cr_cfg <- ann_cfg$custom_reference
customref_info <- list(attempted = FALSE)
if (!is.null(cr_cfg$path)) {
  cr_path <- resolve_path(cr_cfg$path, cfg)
  if (!file.exists(cr_path)) emit_error("annotation.custom_reference.path not found", list(path = cr_path))
  cr_name <- safe_name(sub("\\.rds$", "", basename(cr_path), ignore.case = TRUE))
  cr_cols <- as.character(unlist(cr_cfg$label_columns))
  cr_dims <- seq_len(as.integer(cr_cfg$dims %||% 50))
  log_msg(sprintf("Mapping to custom reference %s (%s)", cr_name, paste(cr_cols, collapse = ", ")))
  customref_info <- tryCatch({
    ref <- readRDS(cr_path)
    if (length(miss <- setdiff(cr_cols, colnames(ref@meta.data))))
      emit_error("annotation.custom_reference.label_columns not in the reference", list(missing = miss,
                 available = colnames(ref@meta.data)))
    # Light query from raw counts; the reference's SCT model recomputes residuals.
    q <- CreateSeuratObject(counts = LayerData(obj, assay = "RNA", layer = "counts"), min.cells = 0, min.features = 0)
    ref_feats <- VariableFeatures(ref, assay = DefaultAssay(ref))
    present <- intersect(ref_feats, rownames(q))
    anchors <- FindTransferAnchors(reference = ref, query = q, normalization.method = "SCT",
                                   reference.reduction = "pca", dims = cr_dims, recompute.residuals = TRUE,
                                   verbose = FALSE)
    q <- TransferData(anchorset = anchors, reference = ref, query = q,
                      refdata = setNames(as.list(cr_cols), cr_cols), dims = cr_dims, verbose = FALSE)
    out <- data.frame(barcode = colnames(obj))
    per_col <- list()
    for (cc in cr_cols) {
      col <- paste0("customref_", cr_name, "_", safe_name(cc))
      lab <- q@meta.data[colnames(obj), paste0("predicted.", cc)]
      sc <- q@meta.data[colnames(obj), paste0("predicted.", cc, ".score")]
      obj[[col]] <- lab
      obj[[paste0(col, "_score")]] <- sc
      out[[cc]] <- lab
      out[[paste0(cc, "_score")]] <- round(sc, 4)
      method_label_cols <- c(method_label_cols, col)
      per_col[[cc]] <- list(label_column = col, median_score = round(median(sc, na.rm = TRUE), 3),
                            low_confidence_pct = round(100 * mean(sc < 0.5, na.rm = TRUE), 1),
                            cell_types_found = I(sort(unique(lab))))
    }
    csv <- file.path(ann_dir, paste0("customref_", cr_name, ".csv"))
    write.csv(out, csv, row.names = FALSE)
    rm(ref, q, anchors)
    list(attempted = TRUE, status = "success", reference = cr_path, name = cr_name, dims = length(cr_dims),
         n_reference_features = length(ref_feats), n_reference_features_present = length(present),
         gene_coverage_pct = round(100 * length(present) / max(length(ref_feats), 1), 1),
         label_columns = per_col, csv = rel_out(csv, cfg))
  }, error = function(e) {
    warn_msg(sprintf("Custom reference mapping failed: %s", e$message))
    list(attempted = TRUE, status = "failed", reference = cr_path, error = e$message)
  })
}
annotation_methods$custom_reference <- customref_info

# --- Cluster markers (from the 06_cluster run) ---
markers_csv_path <- file.path(clust_dir, "cluster_markers.csv")
markers <- if (file.exists(markers_csv_path)) read.csv(markers_csv_path, stringsAsFactors = FALSE) else NULL
markers_info <- list()
top_markers_per_cluster <- list()
if (!is.null(markers) && nrow(markers) > 0) {
  for (cl in unique(markers$cluster)) {
    cl_markers <- markers[markers$cluster == cl, ]
    cl_markers <- cl_markers[order(-cl_markers$avg_log2FC), ]
    top_markers_per_cluster[[as.character(cl)]] <- I(head(cl_markers$gene, 5))
  }
  markers_info <- list(n_markers = nrow(markers), csv = rel_out(markers_csv_path, cfg),
                       top_markers_per_cluster = top_markers_per_cluster)
} else {
  warn_msg("The 06_cluster run has no cluster markers")
  markers_info <- list(n_markers = 0L)
}

# --- Marker gene scoring from config ---
log_msg("Scoring config marker genes")

marker_score_info <- list(attempted = FALSE)
flat_markers <- list()
if (!is.null(ann_cfg$markers$genes)) {
  # Cell type -> genes present in the data
  flat_markers <- list()
  for (ct in names(panel <- flatten_markers(ann_cfg$markers$genes))) {
    present <- intersect(panel[[ct]], rownames(obj))
    if (length(present) == 0) {
      warn_msg(sprintf("marker cell type '%s': none of its genes are in the data; not scored", ct))
      next
    }
    flat_markers[[ct]] <- present
  }

  if (length(flat_markers) > 0) {
    for (ct in names(flat_markers)) {
      score_name <- paste0("score_", gsub("[^[:alnum:]]", "_", ct))
      obj <- tryCatch({
        AddModuleScore(obj, features = list(flat_markers[[ct]]),
                       name = score_name, seed = cfg$project$seed)
      }, error = function(e) {
        warn_msg(sprintf("module score failed for %s: %s", ct, e$message))
        obj
      })
    }

    # Map score columns back to the cell-type names as written in the config
    # (score column names are sanitized), so override_azimuth matches them.
    score_map <- setNames(names(flat_markers), paste0("score_", gsub("[^[:alnum:]]", "_", names(flat_markers)), "1"))
    score_cols <- intersect(names(score_map), colnames(obj@meta.data))
    if (length(score_cols) > 0) {
      cluster_scores <- aggregate(obj@meta.data[, score_cols, drop = FALSE],
                                  by = list(cluster = as.character(Idents(obj))), FUN = mean)

      min_marker_score <- ann_cfg$markers$min_score %||% 0
      marker_labels <- character(nrow(cluster_scores))
      marker_best <- marker_margin <- numeric(nrow(cluster_scores))
      for (i in seq_len(nrow(cluster_scores))) {
        scores <- sort(setNames(as.numeric(cluster_scores[i, score_cols]), score_cols), decreasing = TRUE)
        # No label unless the best module score clears annotation.markers.min_score
        marker_labels[i] <- if (scores[1] > min_marker_score) score_map[[names(scores)[1]]] else NA_character_
        marker_best[i] <- scores[1]
        # Margin over the runner-up: a near-zero margin means an arbitrary pick.
        marker_margin[i] <- if (length(scores) > 1) scores[1] - scores[2] else NA_real_
      }
      names(marker_labels) <- names(marker_best) <- names(marker_margin) <- cluster_scores$cluster

      marker_score_info <- list(
        attempted = TRUE,
        cell_types_scored = I(names(flat_markers)),
        genes_per_type = as.list(vapply(flat_markers, length, integer(1))),
        marker_based_labels = as.list(marker_labels)
      )
    }
  }
}

# --- Build cluster annotation table ---
log_msg("Building cluster annotation table")

clusters <- levels(Idents(obj))
if (length(clusters) == 0) clusters <- sort(unique(as.character(Idents(obj))))

ann_table <- data.frame(cluster = clusters, stringsAsFactors = FALSE)

if (isTRUE(azimuth_info$status == "success") && "azimuth_label" %in% colnames(obj@meta.data)) {
  cl_ids <- as.character(Idents(obj))
  top_lab <- tapply(obj$azimuth_label, cl_ids, function(x) {
    tt <- sort(table(x), decreasing = TRUE)
    c(label = names(tt)[1], frac = round(unname(tt[1]) / length(x), 3))
  })
  ann_table$azimuth_label <- vapply(top_lab[ann_table$cluster], function(v) v[["label"]], "")
  # Share of the cluster carrying that label: low = a mixed cluster.
  ann_table$azimuth_label_fraction <- as.numeric(vapply(top_lab[ann_table$cluster], function(v) v[["frac"]], ""))
  sc <- obj@meta.data[[paste0(azimuth_info$primary_label_column, ".score")]]
  if (!is.null(sc)) {
    ann_table$azimuth_median_score <- round(as.numeric(tapply(sc, cl_ids, median)[ann_table$cluster]), 3)
    ann_table$azimuth_low_conf_pct <- round(100 * as.numeric(tapply(sc < 0.5, cl_ids, mean)[ann_table$cluster]), 1)
  }
} else {
  ann_table$azimuth_label <- NA
}

if ("singler_label" %in% colnames(obj@meta.data)) {   # first successful SingleR reference
  sr_labels <- tapply(obj$singler_label, as.character(Idents(obj)), function(x) {
    tt <- sort(table(x), decreasing = TRUE)
    names(tt)[1]
  })
  ann_table$singler_label <- sr_labels[ann_table$cluster]
} else {
  ann_table$singler_label <- NA
}

if (isTRUE(marker_score_info$attempted) && !is.null(marker_score_info$marker_based_labels)) {
  ann_table$marker_based_label <- unlist(marker_score_info$marker_based_labels[ann_table$cluster])
  ann_table$marker_score <- round(unname(marker_best[ann_table$cluster]), 3)
  ann_table$marker_margin <- round(unname(marker_margin[ann_table$cluster]), 3)
} else {
  ann_table$marker_based_label <- NA
}

ann_table$top_markers <- sapply(ann_table$cluster, function(cl) {
  if (cl %in% names(top_markers_per_cluster)) {
    paste(top_markers_per_cluster[[cl]], collapse = ", ")
  } else {
    NA_character_
  }
})

# Suggested label: prefer Azimuth > marker-based > SingleR > "Unknown", EXCEPT
# for cell types listed in annotation.markers.override_azimuth — those are cell types the
# Azimuth reference cannot represent (it has no class for them), so a
# marker-based call for them wins over Azimuth.
override_types <- unlist(ann_cfg$markers$override_azimuth)
is_label <- function(x) !is.na(x) && !(x %in% c("None", "none", ""))
pick_label <- function(az, mk, sr) {
  if (is_label(mk) && mk %in% override_types) return(c(mk, "markers (override)"))
  if (is_label(az)) return(c(az, "azimuth"))
  if (is_label(mk)) return(c(mk, "markers"))
  if (is_label(sr)) return(c(sr, "singler"))
  c("Unknown", "none")
}
picked <- mapply(pick_label, ann_table$azimuth_label, ann_table$marker_based_label,
                 ann_table$singler_label, USE.NAMES = FALSE)
ann_table$suggested_label <- picked[1, ]
ann_table$label_source <- picked[2, ]

cluster_counts <- table(as.character(Idents(obj)))
ann_table$n_cells <- as.integer(cluster_counts[ann_table$cluster])

ann_csv_path <- file.path(outdir, "cluster_annotations.csv")
write.csv(ann_table, ann_csv_path, row.names = FALSE)
log_msg(sprintf("Wrote cluster annotations to %s", ann_csv_path))

# --- Add suggested labels to object metadata ---
label_map <- setNames(ann_table$suggested_label, ann_table$cluster)
obj$cell_type <- unname(label_map[as.character(Idents(obj))])
# Cluster ids at the default resolution, kept for plots after Idents changes.
cluster_ids <- as.character(Idents(obj))
cluster_levels <- unique(cluster_ids)
cluster_levels <- if (all(grepl("^[0-9]+$", cluster_levels)))
  as.character(sort(as.integer(cluster_levels))) else sort(cluster_levels)

# --- Save final object ---
rds_path <- save_run_rds(obj, cfg, "seurat_final.rds")

# --- Cell type composition ---
composition_overall <- as.list(sort(table(obj$cell_type), decreasing = TRUE))
composition_overall <- lapply(composition_overall, as.integer)

cell_type_counts_by_sample <- list()
for (sid in unique(obj$sample_id)) {
  ct <- table(obj$cell_type[obj$sample_id == sid])
  cell_type_counts_by_sample[[sid]] <- as.list(setNames(as.integer(ct), names(ct)))
}

# --- Removed cells (annotation.removed_cells) ---
# Same Azimuth reference and marker panel, on the cells the 04_filter run removed, to see
# whether they are debris or real cell types the thresholds cut.
low_conf <- 0.5
removed_min_umis <- 10   # cells with fewer UMIs are listed but not annotated
removed_info <- list(enabled = isTRUE(ann_cfg$removed_cells))
removed_df <- NULL
if (isTRUE(ann_cfg$removed_cells)) {
  fo_path <- file.path(ancestor_dir(cfg, "04_filter"), "seurat_filtered_out.rds")
  if (!file.exists(fo_path)) {
    log_msg("The 04_filter run removed no cells; nothing to annotate")
    removed_info$status <- "skipped"
    removed_info$reason <- "no removed cells"
  } else {
    log_msg("Annotating removed cells")
    ro <- read_ancestor_rds(cfg, "04_filter", "seurat_filtered_out.rds")
    removed_df <- data.frame(barcode = colnames(ro), sample_id = ro$sample_id,
                             removal_reason = ro$removal_reason, nCount_RNA = ro$nCount_RNA,
                             stringsAsFactors = FALSE)
    # Near-empty barcodes carry no type information and can break Azimuth's SCTransform.
    ro <- ro[, ro$nCount_RNA >= removed_min_umis]
    DefaultAssay(ro) <- "RNA"
    ro <- NormalizeData(ro, verbose = FALSE)
    idx <- match(colnames(ro), removed_df$barcode)
    if (isTRUE(azimuth_info$status == "success") && ncol(ro) > 0) {
      az <- tryCatch(RunAzimuth(ro, reference = ann_cfg$azimuth$reference)@meta.data, error = function(e) {
        warn_msg(sprintf("Azimuth on removed cells failed: %s", e$message)); NULL })
      if (!is.null(az) && azimuth_info$primary_label_column %in% colnames(az)) {
        az <- az[colnames(ro), , drop = FALSE]
        removed_df$azimuth_label <- NA_character_
        removed_df$azimuth_score <- NA_real_
        removed_df$azimuth_label[idx] <- resolve_azimuth_label(az, unlist(azimuth_info$label_columns),
                                                               azimuth_info$primary_label_column)$label
        removed_df$azimuth_score[idx] <- az[[paste0(azimuth_info$primary_label_column, ".score")]]
      }
      rm(az)
    }
    if (length(flat_markers) && ncol(ro) > 0) {
      sc_mat <- sapply(names(flat_markers), function(ct) {
        g <- intersect(flat_markers[[ct]], rownames(ro))
        r <- tryCatch(AddModuleScore(ro, features = list(g), name = "mscore", seed = cfg$project$seed)$mscore1,
                      error = function(e) rep(NA_real_, ncol(ro)))
        r
      })
      sc_mat <- matrix(sc_mat, nrow = ncol(ro), dimnames = list(NULL, names(flat_markers)))
      best <- apply(sc_mat, 1, function(v) if (all(is.na(v))) NA_integer_ else which.max(v))
      best_score <- sc_mat[cbind(seq_len(nrow(sc_mat)), ifelse(is.na(best), 1L, best))]
      min_marker_score <- ann_cfg$markers$min_score %||% 0
      removed_df$marker_best <- NA_character_
      removed_df$marker_best_score <- NA_real_
      removed_df$marker_best[idx] <- ifelse(!is.na(best) & best_score > min_marker_score, colnames(sc_mat)[best], NA)
      removed_df$marker_best_score[idx] <- round(ifelse(is.na(best), NA, best_score), 3)
    }
    removed_csv_path <- file.path(outdir, "filtered_out_annotation.csv")
    write.csv(removed_df, removed_csv_path, row.names = FALSE)
    top_counts <- function(x, n = 20) {
      tt <- sort(table(x, useNA = "no"), decreasing = TRUE)
      as.list(setNames(as.integer(head(tt, n)), names(head(tt, n))))
    }
    by_reason <- lapply(split(removed_df, removed_df$removal_reason), function(d) c(
      list(n_cells = nrow(d), n_not_annotated_low_umi = sum(d$nCount_RNA < removed_min_umis)),
      if (!is.null(d$azimuth_label)) list(
        azimuth_median_score = round(median(d$azimuth_score, na.rm = TRUE), 3),
        azimuth_low_confidence_pct = round(100 * mean(d$azimuth_score < low_conf, na.rm = TRUE), 1),
        azimuth_label_counts = top_counts(d$azimuth_label)),
      if (!is.null(d$marker_best)) list(
        marker_best_counts = top_counts(d$marker_best),
        marker_no_label = sum(is.na(d$marker_best) & d$nCount_RNA >= removed_min_umis))
    ))
    removed_info <- c(removed_info, list(status = "success", n_cells = nrow(removed_df), by_reason = by_reason,
                                         csv = rel_out(removed_csv_path, cfg)))
    rm(ro)
  }
}

# --- Annotation by sample metadata (annotation.compare_by) ---
# A reference built from one condition can label another condition's cells
# differently or with lower confidence; compare score and composition by level.
compare_min_cells <- 20
compare_score_gap <- 0.2
compare_info <- list()
cmp_label_col <- if ("azimuth_label" %in% colnames(obj@meta.data)) "azimuth_label" else "cell_type"
cmp_score <- if (isTRUE(azimuth_info$status == "success"))
  obj@meta.data[[paste0(azimuth_info$primary_label_column, ".score")]] else NULL
for (col in unlist(ann_cfg$compare_by)) {
  if (!col %in% colnames(obj@meta.data)) {
    warn_msg(sprintf("annotation.compare_by: '%s' is not a metadata column; skipped", col))
    next
  }
  g <- as.character(obj@meta.data[[col]])
  lab <- as.character(obj@meta.data[[cmp_label_col]])
  lvls <- sort(unique(g[!is.na(g)]))
  per_level <- lapply(setNames(nm = lvls), function(l) {
    i <- !is.na(g) & g == l
    c(list(n_cells = sum(i), n_samples = length(unique(obj$sample_id[i]))),
      if (!is.null(cmp_score)) list(
        azimuth_median_score = round(median(cmp_score[i], na.rm = TRUE), 3),
        azimuth_low_confidence_pct = round(100 * mean(cmp_score[i] < low_conf, na.rm = TRUE), 1)),
      list(label_pct = as.list(round(100 * sort(prop.table(table(lab[i])), decreasing = TRUE), 2))))
  })
  # Per-label median score by level, for labels with enough cells in every level;
  # a large gap suggests the label means something different across levels.
  label_scores <- list(); flagged <- character(0)
  if (!is.null(cmp_score) && length(lvls) > 1) {
    for (lb in sort(unique(lab))) {
      n_by <- vapply(lvls, function(l) sum(lab == lb & !is.na(g) & g == l), integer(1))
      if (any(n_by < compare_min_cells)) next
      med <- vapply(lvls, function(l) median(cmp_score[lab == lb & !is.na(g) & g == l], na.rm = TRUE), numeric(1))
      label_scores[[lb]] <- as.list(round(setNames(med, lvls), 3))
      if (diff(range(med)) >= compare_score_gap) flagged <- c(flagged, lb)
    }
  }
  if (length(flagged))
    warn_msg(sprintf("compare_by %s: median Azimuth score differs by >= %.1f between levels for: %s",
                     col, compare_score_gap, paste(flagged, collapse = ", ")))
  compare_info[[col]] <- list(label_column = cmp_label_col, levels = per_level,
                              label_median_score_by_level = label_scores,
                              labels_with_score_gap = I(flagged))
}

# --- Marker genes removed by ambient correction (02_ambient run on the lineage) ---
# Share of each marker cell type's UMIs (over its genes) that correction removed,
# per sample. A high value means correction is taking that cell type's own signal.
marker_removal <- NULL
amb_csv <- file.path(ancestor_dir(cfg, "02_ambient"), "ambient_removed_by_gene.csv")
if (length(flat_markers) && file.exists(amb_csv)) {
  g <- read.csv(amb_csv, stringsAsFactors = FALSE)
  marker_removal <- lapply(split(g, g$sample_id), function(d) lapply(flat_markers, function(genes) {
    d <- d[d$gene %in% genes, ]
    if (!nrow(d) || sum(d$umis_before) == 0) NA_real_ else round(100 * sum(d$umis_removed) / sum(d$umis_before), 1)
  }))
}

# --- JSON summary ---
summary_data <- c(
  summary_header("07_annotate", cfg,
                 packages = c(if (!is.null(ann_cfg$azimuth$reference)) c("Azimuth", "SeuratData"),
                              if (length(singler_refs)) c("SingleR", "celldex"))),
  list(
    inputs = list(rds = rel_out(input_rds, cfg)),
    params = list(azimuth = ann_cfg$azimuth,
                  singler = list(reference = I(singler_refs), labels = singler_level),
                  celltypist = list(models = I(ct_models), majority_voting = isTRUE(ann_cfg$celltypist$majority_voting)),
                  custom_reference = list(path = cr_cfg$path, label_columns = I(as.character(unlist(cr_cfg$label_columns))),
                                          dims = cr_cfg$dims),
                  markers = list(override_azimuth = I(override_types %||% character(0)),
                                 min_score = ann_cfg$markers$min_score,
                                 n_cell_types = length(unlist(lapply(ann_cfg$markers$genes, names)))),
                  low_confidence_score = low_conf,
                  removed_cells = isTRUE(ann_cfg$removed_cells), removed_min_umis = removed_min_umis,
                  compare_by = I(as.character(unlist(ann_cfg$compare_by))),
                  compare_min_cells = compare_min_cells, compare_score_gap = compare_score_gap),
    n_cells = ncol(obj),
    n_clusters = length(clusters),
    methods = c(annotation_methods, list(markers = marker_score_info)),
    cluster_markers = markers_info,
    cluster_annotations = lapply(seq_len(nrow(ann_table)), function(i) {
      as.list(ann_table[i, ])
    }),
    cell_type_composition = composition_overall,
    cell_type_counts_by_sample = cell_type_counts_by_sample,
    removed_cells = removed_info,
    compare_by = compare_info,
    marker_genes_removed_by_ambient_pct = marker_removal,
    outputs = list(
      rds = rel_out(rds_path, cfg),
      cluster_annotations_csv = rel_out(ann_csv_path, cfg),
      annotation_dir = rel_out(ann_dir, cfg)
    ),
    status = "success"
  )
)

write_json_summary(summary_data, file.path(outdir, STEPS[["07_annotate"]]$summary))

# --- Plots ---
if (args$make_plots) {
  log_msg("Generating annotation plots")

  # Bottom legend that wraps into rows (a single side column overflows with many
  # cell types). ncol caps columns per row; height grows with the row count.
  legend_cols <- 5
  bottom_legend <- function(n) {
    list(theme(legend.position = "bottom", legend.text = element_text(size = 8)),
         guides(color = guide_legend(ncol = legend_cols, byrow = TRUE,
                                     override.aes = list(size = 3))),
         height = 7 + 0.45 * ceiling(n / legend_cols))
  }

  # UMAP colored by cell type (suggested label) — keep on-plot labels + add legend
  Idents(obj) <- "cell_type"
  ct_n <- length(unique(obj$cell_type))
  bl <- bottom_legend(ct_n)
  p <- DimPlot(obj, reduction = "umap", label = TRUE, repel = TRUE, label.size = 3) +
    ggtitle("UMAP — cell type annotations") + bl[[1]] + bl[[2]]
  save_plot(p, "07_umap_cell_types.png", cfg, width = 9, height = bl$height)

  # UMAP for each config-specified Azimuth level
  cfg_levels <- unlist(ann_cfg$azimuth$plot_levels)
  if (length(cfg_levels)) {
    for (lvl in cfg_levels) {
      col_name <- paste0("predicted.", lvl)
      if (col_name %in% colnames(obj@meta.data)) {
        n_types <- length(unique(obj@meta.data[[col_name]]))
        bl <- bottom_legend(n_types)
        p <- DimPlot(obj, reduction = "umap", group.by = col_name,
                     label = TRUE, repel = TRUE, label.size = 2.5) +
          ggtitle(paste0("UMAP — Azimuth ", lvl)) + bl[[1]] + bl[[2]]
        save_plot(p, paste0("07_umap_azimuth_", lvl, ".png"), cfg,
                  width = 10, height = bl$height)
      }
    }
  }

  # UMAP for each label column from SingleR, CellTypist and the custom reference
  for (col in unique(method_label_cols)) {
    n_types <- length(unique(obj@meta.data[[col]]))
    bl <- bottom_legend(n_types)
    p <- DimPlot(obj, reduction = "umap", group.by = col, label = n_types <= 40, repel = TRUE, label.size = 2.5) +
      ggtitle(paste0("UMAP — ", col)) + bl[[1]] + bl[[2]]
    save_plot(p, paste0("07_umap_", col, ".png"), cfg, width = 10, height = bl$height)
  }

  # UMAP split by sample, colored by the per-cell primary-level label (Azimuth,
  # "None" filled from the coarser level), else the per-cluster suggested label.
  # Mirrors 04_umap_split_sample.png; a legend replaces on-plot labels.
  split_col <- if ("azimuth_label" %in% colnames(obj@meta.data)) "azimuth_label" else "cell_type"
  n_samples <- length(unique(obj$sample_id))
  if (n_samples > 1) {
    split_ncol <- min(4L, n_samples)
    split_nrow <- ceiling(n_samples / split_ncol)
    n_types <- length(unique(obj@meta.data[[split_col]]))
    bl <- bottom_legend(n_types)
    p <- DimPlot(obj, reduction = "umap", group.by = split_col, split.by = "sample_id",
                 ncol = split_ncol) +
      ggtitle(sprintf("UMAP split by sample: %s",
                      if (split_col == "azimuth_label") sub("^predicted\\.", "", azimuth_info$primary_label_column)
                      else "suggested cell type")) +
      bl[[1]] + bl[[2]]
    save_plot(p, "07_umap_split_sample_cell_type.png", cfg,
              width = 3.3 * split_ncol + 0.5, height = 3.3 * split_nrow + 0.6 + (bl$height - 7))
  }

  # Stacked bar: cell type composition per sample
  comp_df <- as.data.frame(table(obj$sample_id, obj$cell_type))
  colnames(comp_df) <- c("sample", "cell_type", "count")
  totals <- tapply(comp_df$count, comp_df$sample, sum)
  comp_df$proportion <- comp_df$count / totals[comp_df$sample]

  p <- ggplot(comp_df, aes(x = sample, y = proportion, fill = cell_type)) +
    geom_col() +
    labs(title = "Cell type composition per sample", x = NULL, y = "Proportion", fill = "Cell type") +
    theme_bw(base_size = 12) +
    theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p, "07_composition_per_sample.png", cfg, width = 10, height = 6)

  if (isTRUE(azimuth_info$status == "success") && !is.null(azimuth_info$primary_label_column)) {
    md <- obj@meta.data

    # Azimuth confidence: primary-level prediction score on the UMAP and per cluster
    score_col <- paste0(azimuth_info$primary_label_column, ".score")
    if (score_col %in% colnames(md)) {
      p <- FeaturePlot(obj, features = score_col, reduction = "umap", raster = FALSE) +
        scale_colour_viridis_c(limits = c(0, 1), name = "score") +
        ggtitle(sprintf("Azimuth prediction score (%s)", sub("^predicted\\.", "", azimuth_info$primary_label_column)))
      save_plot(p, "07_umap_azimuth_score.png", cfg, width = 8, height = 7)

      sc_df <- data.frame(cluster = factor(cluster_ids, levels = cluster_levels), score = md[[score_col]])
      p <- ggplot(sc_df, aes(x = cluster, y = score, fill = cluster)) +
        geom_violin(scale = "width", linewidth = 0.2) +
        geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
        labs(title = "Azimuth prediction score per cluster", x = "Cluster", y = "Score") +
        theme_bw(base_size = 11) + NoLegend()
      save_plot(p, "07_violin_azimuth_score_per_cluster.png", cfg,
                width = max(8, 0.35 * length(cluster_levels) + 2), height = 5)
    }

    # Cluster x label heatmap: share of each cluster's cells per Azimuth label at
    # the last of azimuth.plot_levels (else the primary level), so mixed clusters show.
    lvl_cols <- intersect(paste0("predicted.", unlist(ann_cfg$azimuth$plot_levels)), colnames(md))
    hm_col <- if (length(lvl_cols)) lvl_cols[length(lvl_cols)] else azimuth_info$primary_label_column
    frac <- 100 * prop.table(table(factor(cluster_ids, levels = cluster_levels), md[[hm_col]]), 1)
    frac <- frac[, apply(frac, 2, max) >= 3, drop = FALSE]   # drop labels never >= 3% of a cluster
    if (ncol(frac) > 0) {
      main_lab <- apply(frac, 1, which.max)
      row_ord <- order(main_lab, -apply(frac, 1, max))
      col_ord <- order(apply(frac[row_ord, , drop = FALSE], 2, which.max))
      n_by_cluster <- table(factor(cluster_ids, levels = cluster_levels))
      row_lab <- sprintf("%s (n=%d)", rownames(frac), as.integer(n_by_cluster[rownames(frac)]))
      hm <- as.data.frame(as.table(frac))
      colnames(hm) <- c("cluster", "label", "pct")
      hm$cluster <- factor(row_lab[match(hm$cluster, rownames(frac))], levels = rev(row_lab[row_ord]))
      hm$label <- factor(hm$label, levels = colnames(frac)[col_ord])
      p <- ggplot(hm, aes(x = label, y = cluster, fill = pct)) +
        geom_tile(color = "grey90") +
        geom_text(aes(label = ifelse(pct >= 10, sprintf("%.0f", pct), "")), size = 2.6) +
        scale_fill_gradient(low = "white", high = "navy", limits = c(0, 100), name = "% of\ncluster") +
        labs(title = sprintf("Cluster x Azimuth label (%s)", sub("^predicted\\.", "", hm_col)),
             x = NULL, y = NULL) +
        theme_minimal(base_size = 10) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      save_plot(p, "07_heatmap_cluster_label.png", cfg,
                width = max(8, 0.3 * ncol(frac) + 4), height = max(6, 0.28 * nrow(frac) + 2.5))
    }
  }

  # Marker-gene UMIs removed by ambient correction, per sample
  if (!is.null(marker_removal)) {
    mk <- do.call(rbind, lapply(names(marker_removal), function(sid) data.frame(
      sample = sid, cell_type = names(marker_removal[[sid]]),
      pct = vapply(marker_removal[[sid]], function(x) x %||% NA_real_, numeric(1)), stringsAsFactors = FALSE)))
    mk$cell_type <- factor(mk$cell_type, levels = rev(names(flat_markers)))
    p <- ggplot(mk, aes(x = sample, y = cell_type, fill = pct)) +
      geom_tile(color = "grey90") +
      geom_text(aes(label = ifelse(is.na(pct), "", sprintf("%.0f", pct))), size = 2.6) +
      scale_fill_gradient(low = "white", high = "firebrick", limits = c(0, 100), na.value = "grey95",
                          name = "% of marker\nUMIs removed") +
      labs(title = "Marker-gene UMIs removed by ambient correction", x = NULL, y = NULL) +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    save_plot(p, "07_ambient_marker_genes_removed.png", cfg,
              width = max(6, 0.5 * length(unique(mk$sample)) + 4), height = max(5, 0.25 * nlevels(mk$cell_type) + 2))
  }

  # Removed cells: label composition by removal reason (Azimuth and markers)
  if (!is.null(removed_df) && nrow(removed_df)) {
    lump <- function(x, n = 15) {
      x[is.na(x)] <- "(no label)"
      top <- names(sort(table(x), decreasing = TRUE))[seq_len(min(n, length(unique(x))))]
      ifelse(x %in% top, x, "Other")
    }
    parts <- list()
    if (!is.null(removed_df$azimuth_label))
      parts$az <- data.frame(reason = removed_df$removal_reason, label = lump(removed_df$azimuth_label),
                             method = "Azimuth")
    if (!is.null(removed_df$marker_best))
      parts$mk <- data.frame(reason = removed_df$removal_reason, label = lump(removed_df$marker_best),
                             method = "Marker panel")
    if (length(parts)) {
      rd <- do.call(rbind, parts)
      p <- ggplot(rd, aes(x = reason, fill = label)) +
        geom_bar(position = "fill") +
        facet_wrap(~method) +
        labs(title = sprintf("Labels of removed cells (n = %s)", format(nrow(removed_df), big.mark = ",")),
             x = "Removal reason", y = "Proportion", fill = NULL) +
        theme_bw(base_size = 11) +
        theme(legend.text = element_text(size = 8))
      save_plot(p, "07_removed_cells_labels.png", cfg, width = 11, height = 6)
    }
  }

  # Annotation by sample metadata: score and composition per level
  for (col in names(compare_info)) {
    g <- as.character(obj@meta.data[[col]])
    keep <- !is.na(g)
    lab <- as.character(obj@meta.data[[cmp_label_col]])
    cd <- as.data.frame(table(level = g[keep], label = lab[keep]))
    cd$prop <- cd$Freq / ave(cd$Freq, cd$level, FUN = sum)
    p <- ggplot(cd, aes(x = level, y = prop, fill = label)) +
      geom_col() +
      labs(title = sprintf("Composition by %s (%s)", col, cmp_label_col), x = col, y = "Proportion", fill = NULL) +
      theme_bw(base_size = 11) +
      theme(legend.text = element_text(size = 8))
    save_plot(p, sprintf("07_compare_%s_composition.png", col), cfg, width = 10, height = 6)
    if (!is.null(cmp_score)) {
      sd <- data.frame(level = g[keep], label = lab[keep], score = cmp_score[keep])
      ok_labels <- names(compare_info[[col]]$label_median_score_by_level)
      sd <- sd[sd$label %in% ok_labels, ]
      if (nrow(sd)) {
        p <- ggplot(sd, aes(x = label, y = score, fill = level)) +
          geom_boxplot(outlier.size = 0.2, linewidth = 0.2) +
          geom_hline(yintercept = low_conf, linetype = "dashed", color = "red") +
          labs(title = sprintf("Azimuth score by label and %s (labels with >= %d cells per level)",
                               col, compare_min_cells), x = NULL, y = "Score", fill = col) +
          theme_bw(base_size = 10) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        save_plot(p, sprintf("07_compare_%s_score.png", col), cfg,
                  width = max(8, 0.35 * length(ok_labels) + 3), height = 5.5)
      }
    }
  }

  # Dot plot of marker genes (from config) x clusters
  if (!is.null(ann_cfg$markers$genes)) {
    dot_genes <- c()
    for (compartment in names(ann_cfg$markers$genes)) {
      for (cell_type in names(ann_cfg$markers$genes[[compartment]])) {
        genes <- ann_cfg$markers$genes[[compartment]][[cell_type]]
        present <- intersect(genes, rownames(obj))
        dot_genes <- c(dot_genes, head(present, 2))
      }
    }
    dot_genes <- unique(dot_genes)

    if (length(dot_genes) > 0) {
      # Clusters at the default resolution (the same ones the labels and markers use)
      Idents(obj) <- factor(cluster_ids, levels = cluster_levels)

      p <- DotPlot(obj, features = dot_genes) +
        coord_flip() +
        theme_bw(base_size = 10) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle("Marker gene expression by cluster")
      save_plot(p, "07_dotplot_markers.png", cfg, width = 12, height = max(6, length(dot_genes) * 0.25))
    }
  }

  # Heatmap of top markers per cluster
  tryCatch({
    if (!is.null(markers) && nrow(markers) > 0) {
      top_genes <- c()
      for (cl in unique(as.character(markers$cluster))) {
        cl_markers <- markers[as.character(markers$cluster) == cl, ]
        cl_markers <- cl_markers[order(-cl_markers$avg_log2FC), ]
        top_genes <- c(top_genes, head(cl_markers$gene, 3))
      }
      top_genes <- unique(top_genes)

      if (length(top_genes) > 0) {
        heatmap_assay <- if ("SCT" %in% Assays(obj)) "SCT" else DefaultAssay(obj)
        DefaultAssay(obj) <- heatmap_assay
        Idents(obj) <- factor(cluster_ids, levels = cluster_levels)   # default resolution
        available_genes <- rownames(GetAssayData(obj, layer = "scale.data"))
        heatmap_genes <- intersect(top_genes, available_genes)
        if (length(heatmap_genes) > 0) {
          p <- DoHeatmap(obj, features = heatmap_genes, size = 3) +
            ggtitle("Top 3 markers per cluster") +
            theme(axis.text.y = element_text(size = 6))
          save_plot(p, "07_heatmap_top_markers.png", cfg, width = 12, height = 10)
        } else {
          log_msg("Heatmap skipped: no marker genes found in scale.data")
        }
      }
    }
  }, error = function(e) {
    log_msg(sprintf("Heatmap skipped: %s", e$message))
  })
}

log_msg("Step 07 complete")
