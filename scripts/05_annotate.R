#!/usr/bin/env Rscript
# 05_annotate.R
# Annotate clusters: optional Azimuth, optional SingleR, marker gene scoring,
# FindAllMarkers. Builds a per-cluster annotation table with a suggested label.
# Writes: seurat_final.rds, cluster_markers.csv, cluster_annotations.csv, annotation_summary.json

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
cfg <- load_config(args$config_path, step = "05_annotate", summary_file = "annotation_summary.json")
set.seed(cfg$project$seed)

options(future.globals.maxSize = 4 * 1024^3)  # 4 GiB

outdir <- output_dir(cfg)

# --- Load input ---
input_rds <- file.path(outdir, "seurat_clustered.rds")
if (!file.exists(input_rds)) {
  emit_error("Input RDS not found. Run step 04 first.", list(expected = input_rds))
}

log_msg("Loading seurat_clustered.rds")
obj <- readRDS(input_rds)
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
        # Resolved per-cell label: references such as lungref use the literal
        # "None" where a cell type has no subdivision at this level (e.g.
        # smooth muscle at ann_level_3). Fill those from the nearest coarser
        # level (label columns run coarse -> fine) instead of treating "None"
        # as a cell type.
        lab <- as.character(obj@meta.data[[primary_col]])
        n_fallback <- 0L
        for (col in rev(head(label_cols, match(primary_col, label_cols) - 1))) {
          miss <- is.na(lab) | lab %in% c("None", "none", "")
          if (!any(miss)) break
          lab[miss] <- as.character(obj@meta.data[[col]][miss])
          n_fallback <- n_fallback + sum(miss & !(lab %in% c("None", "none", "")))
        }
        obj$azimuth_label <- lab
        info$n_cells_label_from_coarser_level <- n_fallback
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
singler_info <- list(attempted = FALSE)

if (!is.null(ann_cfg$singler$reference)) {
  # Reference requested but deps missing -> HARD STOP (don't silently drop annotation).
  .singler_missing <- setdiff(c("SingleR", "celldex", "SingleCellExperiment"),
                              rownames(installed.packages()))
  if (length(.singler_missing) > 0)
    emit_error("annotation.singler.reference is set but required packages are not installed",
               list(missing = .singler_missing,
                    fix = "BiocManager::install(the missing packages)  (or set annotation.singler.reference: null)"))
  {
    log_msg(sprintf("Running SingleR with reference: %s", ann_cfg$singler$reference))
    singler_info <- tryCatch({
      suppressPackageStartupMessages({
        library(SingleR)
        library(celldex)
        library(SingleCellExperiment)
      })
      sce <- as.SingleCellExperiment(obj)
      ref <- do.call(ann_cfg$singler$reference, list())
      results <- SingleR(test = sce, ref = ref, labels = ref$label.main)
      obj$singler_label <- results$labels
      obj$singler_pruned_label <- results$pruned.labels
      list(attempted = TRUE, status = "success", reference = ann_cfg$singler$reference,
           cell_types_found = I(sort(unique(results$labels))))
    }, error = function(e) {
      warn_msg(sprintf("SingleR failed: %s", e$message))
      list(attempted = TRUE, status = "failed", error = e$message)
    })
  }
}
annotation_methods$singler <- singler_info

# --- FindAllMarkers ---
log_msg("Finding cluster markers")

marker_params <- list(only_pos = TRUE, min_pct = 0.25, logfc_threshold = 0.25, test = "wilcox")
markers <- tryCatch({
  FindAllMarkers(obj, only.pos = marker_params$only_pos, min.pct = marker_params$min_pct,
                 logfc.threshold = marker_params$logfc_threshold, test.use = marker_params$test,
                 verbose = FALSE)
}, error = function(e) {
  warn_msg(sprintf("FindAllMarkers failed: %s", e$message))
  NULL
})

markers_info <- list()
top_markers_per_cluster <- list()

if (!is.null(markers) && nrow(markers) > 0) {
  markers_csv_path <- file.path(outdir, "cluster_markers.csv")
  write.csv(markers, markers_csv_path, row.names = FALSE)
  log_msg(sprintf("Found %d markers across %d clusters", nrow(markers), length(unique(markers$cluster))))

  for (cl in unique(markers$cluster)) {
    cl_markers <- markers[markers$cluster == cl, ]
    cl_markers <- cl_markers[order(-cl_markers$avg_log2FC), ]
    top_markers_per_cluster[[as.character(cl)]] <- I(head(cl_markers$gene, 5))
  }

  markers_info <- list(
    n_markers = nrow(markers),
    top_markers_per_cluster = top_markers_per_cluster
  )
} else {
  warn_msg("No cluster markers found")
  markers_info <- list(n_markers = 0L)
}

# --- Marker gene scoring from config ---
log_msg("Scoring config marker genes")

marker_score_info <- list(attempted = FALSE)
if (!is.null(ann_cfg$markers$genes)) {
  # Flatten marker dictionary: cell_type -> genes
  # Accept compartment -> cell type -> genes, or a flat cell type -> genes map.
  entries <- list()
  for (key in names(ann_cfg$markers$genes)) {
    val <- ann_cfg$markers$genes[[key]]
    if (is.list(val)) {
      for (ct in names(val)) entries[[length(entries) + 1]] <- list(ct = ct, genes = unlist(val[[ct]]))
    } else {
      entries[[length(entries) + 1]] <- list(ct = key, genes = unlist(val))
    }
  }
  flat_markers <- list()
  for (e in entries) {
    if (e$ct %in% names(flat_markers)) {
      warn_msg(sprintf("marker cell type '%s' is defined more than once; keeping the first", e$ct))
      next
    }
    present <- intersect(e$genes, rownames(obj))
    if (length(present) == 0) {
      warn_msg(sprintf("marker cell type '%s': none of its genes are in the data; not scored", e$ct))
      next
    }
    flat_markers[[e$ct]] <- present
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

if (isTRUE(singler_info$status == "success") && "singler_label" %in% colnames(obj@meta.data)) {
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
rds_path <- file.path(outdir, "seurat_final.rds")
log_msg(sprintf("Saving final annotated object to %s", rds_path))
saveRDS(obj, rds_path)

# --- Cell type composition ---
composition_overall <- as.list(sort(table(obj$cell_type), decreasing = TRUE))
composition_overall <- lapply(composition_overall, as.integer)

cell_type_counts_by_sample <- list()
for (sid in unique(obj$sample_id)) {
  ct <- table(obj$cell_type[obj$sample_id == sid])
  cell_type_counts_by_sample[[sid]] <- as.list(setNames(as.integer(ct), names(ct)))
}

# --- JSON summary ---
summary_data <- c(
  summary_header("05_annotate", cfg,
                 packages = c(if (!is.null(ann_cfg$azimuth$reference)) c("Azimuth", "SeuratData"),
                              if (!is.null(ann_cfg$singler$reference)) c("SingleR", "celldex"))),
  list(
    inputs = list(rds = rel_out(input_rds, cfg)),
    params = list(azimuth = ann_cfg$azimuth, singler = ann_cfg$singler,
                  markers = list(override_azimuth = I(override_types %||% character(0)),
                                 min_score = ann_cfg$markers$min_score,
                                 n_cell_types = length(unlist(lapply(ann_cfg$markers$genes, names)))),
                  find_all_markers = marker_params, low_confidence_score = 0.5),
    n_cells = ncol(obj),
    n_clusters = length(clusters),
    methods = c(annotation_methods, list(markers = marker_score_info)),
    cluster_markers = markers_info,
    cluster_annotations = lapply(seq_len(nrow(ann_table)), function(i) {
      as.list(ann_table[i, ])
    }),
    cell_type_composition = composition_overall,
    cell_type_counts_by_sample = cell_type_counts_by_sample,
    outputs = list(
      rds = rel_out(rds_path, cfg),
      cluster_annotations_csv = rel_out(ann_csv_path, cfg),
      cluster_markers_csv = if (exists("markers_csv_path")) rel_out(markers_csv_path, cfg) else NULL
    ),
    status = "success"
  )
)

write_json_summary(summary_data, file.path(outdir, "annotation_summary.json"))

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
  save_plot(p, "05_umap_cell_types.png", cfg, width = 9, height = bl$height)

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
        save_plot(p, paste0("05_umap_azimuth_", lvl, ".png"), cfg,
                  width = 10, height = bl$height)
      }
    }
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
    save_plot(p, "05_umap_split_sample_cell_type.png", cfg,
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
  save_plot(p, "05_composition_per_sample.png", cfg, width = 10, height = 6)

  if (isTRUE(azimuth_info$status == "success") && !is.null(azimuth_info$primary_label_column)) {
    md <- obj@meta.data

    # Azimuth confidence: primary-level prediction score on the UMAP and per cluster
    score_col <- paste0(azimuth_info$primary_label_column, ".score")
    if (score_col %in% colnames(md)) {
      p <- FeaturePlot(obj, features = score_col, reduction = "umap", raster = FALSE) +
        scale_colour_viridis_c(limits = c(0, 1), name = "score") +
        ggtitle(sprintf("Azimuth prediction score (%s)", sub("^predicted\\.", "", azimuth_info$primary_label_column)))
      save_plot(p, "05_umap_azimuth_score.png", cfg, width = 8, height = 7)

      sc_df <- data.frame(cluster = factor(cluster_ids, levels = cluster_levels), score = md[[score_col]])
      p <- ggplot(sc_df, aes(x = cluster, y = score, fill = cluster)) +
        geom_violin(scale = "width", linewidth = 0.2) +
        geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
        labs(title = "Azimuth prediction score per cluster", x = "Cluster", y = "Score") +
        theme_bw(base_size = 11) + NoLegend()
      save_plot(p, "05_violin_azimuth_score_per_cluster.png", cfg,
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
      save_plot(p, "05_heatmap_cluster_label.png", cfg,
                width = max(8, 0.3 * ncol(frac) + 4), height = max(6, 0.28 * nrow(frac) + 2.5))
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
      save_plot(p, "05_dotplot_markers.png", cfg, width = 12, height = max(6, length(dot_genes) * 0.25))
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
          save_plot(p, "05_heatmap_top_markers.png", cfg, width = 12, height = 10)
        } else {
          log_msg("Heatmap skipped: no marker genes found in scale.data")
        }
      }
    }
  }, error = function(e) {
    log_msg(sprintf("Heatmap skipped: %s", e$message))
  })
}

log_msg("Step 05 complete")
