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
cfg <- load_config(args$config_path)
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

if (!is.null(ann_cfg$azimuth_reference)) {
  # Reference requested but dep missing -> HARD STOP (don't silently drop annotation).
  if (!requireNamespace("Azimuth", quietly = TRUE))
    emit_error("annotation.azimuth_reference is set but the 'Azimuth' package is not installed",
               list(fix = "install Azimuth (or set annotation.azimuth_reference: null to skip)"))
  {
    log_msg(sprintf("Running Azimuth with reference: %s", ann_cfg$azimuth_reference))
    azimuth_info <- tryCatch({
      suppressPackageStartupMessages(library(Azimuth))
      obj <- RunAzimuth(obj, reference = ann_cfg$azimuth_reference)

      pred_cols <- grep("^predicted\\.", colnames(obj@meta.data), value = TRUE)
      score_cols <- grep("\\.score$", pred_cols, value = TRUE)
      label_cols <- setdiff(pred_cols, score_cols)
      label_cols <- setdiff(label_cols, "predicted.id")

      info <- list(attempted = TRUE, status = "success", reference = ann_cfg$azimuth_reference,
                   label_columns = label_cols, score_columns = score_cols)

      if (length(label_cols) > 0) {
        cfg_levels <- ann_cfg$azimuth_levels
        if (!is.null(cfg_levels)) {
          cfg_level_cols <- paste0("predicted.", cfg_levels)
          primary_col <- NULL
          for (lp in cfg_level_cols) {
            if (lp %in% label_cols) { primary_col <- lp; break }
          }
          if (is.null(primary_col)) primary_col <- label_cols[1]
        } else {
          primary_col <- label_cols[1]
        }
        info$primary_label_column <- primary_col
        info$cell_types_found <- sort(unique(obj@meta.data[[primary_col]]))
        score_col <- paste0(primary_col, ".score")
        if (score_col %in% colnames(obj@meta.data)) {
          info$median_prediction_score <- round(median(obj@meta.data[[score_col]], na.rm = TRUE), 3)
          info$low_confidence_pct <- round(100 * mean(obj@meta.data[[score_col]] < 0.5, na.rm = TRUE), 1)
        }
      }
      info
    }, error = function(e) {
      log_msg(sprintf("WARNING: Azimuth failed: %s", e$message))
      list(attempted = TRUE, status = "failed", error = e$message)
    })
  }
}
annotation_methods$azimuth <- azimuth_info

# --- SingleR ---
singler_info <- list(attempted = FALSE)

if (!is.null(ann_cfg$singler_reference)) {
  # Reference requested but deps missing -> HARD STOP (don't silently drop annotation).
  .singler_missing <- setdiff(c("SingleR", "celldex", "SingleCellExperiment"),
                              rownames(installed.packages()))
  if (length(.singler_missing) > 0)
    emit_error("annotation.singler_reference is set but required packages are not installed",
               list(missing = .singler_missing,
                    fix = "BiocManager::install(the missing packages)  (or set annotation.singler_reference: null)"))
  {
    log_msg(sprintf("Running SingleR with reference: %s", ann_cfg$singler_reference))
    singler_info <- tryCatch({
      suppressPackageStartupMessages({
        library(SingleR)
        library(celldex)
        library(SingleCellExperiment)
      })
      sce <- as.SingleCellExperiment(obj)
      ref <- do.call(ann_cfg$singler_reference, list())
      results <- SingleR(test = sce, ref = ref, labels = ref$label.main)
      obj$singler_label <- results$labels
      obj$singler_pruned_label <- results$pruned.labels
      list(attempted = TRUE, status = "success", reference = ann_cfg$singler_reference,
           cell_types_found = sort(unique(results$labels)))
    }, error = function(e) {
      log_msg(sprintf("WARNING: SingleR failed: %s", e$message))
      list(attempted = TRUE, status = "failed", error = e$message)
    })
  }
}
annotation_methods$singler <- singler_info

# --- FindAllMarkers ---
log_msg("Finding cluster markers")

markers <- tryCatch({
  FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25,
                 logfc.threshold = 0.25, test.use = "wilcox",
                 verbose = FALSE)
}, error = function(e) {
  log_msg(sprintf("WARNING: FindAllMarkers failed: %s", e$message))
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
    top_markers_per_cluster[[as.character(cl)]] <- head(cl_markers$gene, 5)
  }

  markers_info <- list(
    total_markers = nrow(markers),
    markers_csv = markers_csv_path,
    top_markers_per_cluster = top_markers_per_cluster
  )
} else {
  log_msg("WARNING: No markers found")
  markers_info <- list(total_markers = 0, status = "no markers found")
}

# --- Marker gene scoring from config ---
log_msg("Scoring config marker genes")

marker_score_info <- list()
if (!is.null(ann_cfg$marker_genes)) {
  # Flatten marker dictionary: cell_type -> genes
  flat_markers <- list()
  for (compartment in names(ann_cfg$marker_genes)) {
    for (cell_type in names(ann_cfg$marker_genes[[compartment]])) {
      genes <- ann_cfg$marker_genes[[compartment]][[cell_type]]
      present <- intersect(genes, rownames(obj))
      if (length(present) > 0) {
        flat_markers[[cell_type]] <- present
      }
    }
  }

  if (length(flat_markers) > 0) {
    for (ct in names(flat_markers)) {
      score_name <- paste0("score_", gsub("[^[:alnum:]]", "_", ct))
      obj <- tryCatch({
        AddModuleScore(obj, features = list(flat_markers[[ct]]),
                       name = score_name, seed = cfg$project$seed)
      }, error = function(e) {
        log_msg(sprintf("  Module score failed for %s: %s", ct, e$message))
        obj
      })
    }

    score_cols <- grep("^score_.*1$", colnames(obj@meta.data), value = TRUE)
    if (length(score_cols) > 0) {
      cluster_scores <- aggregate(obj@meta.data[, score_cols, drop = FALSE],
                                  by = list(cluster = as.character(Idents(obj))), FUN = mean)

      marker_labels <- character(nrow(cluster_scores))
      for (i in seq_len(nrow(cluster_scores))) {
        scores <- as.numeric(cluster_scores[i, -1])
        names(scores) <- score_cols
        best <- names(which.max(scores))
        marker_labels[i] <- gsub("^score_|1$", "", best)
      }
      names(marker_labels) <- cluster_scores$cluster

      marker_score_info <- list(
        cell_types_scored = names(flat_markers),
        genes_per_type = sapply(flat_markers, length),
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

if (isTRUE(azimuth_info$status == "success") && !is.null(azimuth_info$primary_label_column)) {
  az_col <- azimuth_info$primary_label_column
  az_labels <- tapply(obj@meta.data[[az_col]], as.character(Idents(obj)), function(x) {
    tt <- sort(table(x), decreasing = TRUE)
    names(tt)[1]
  })
  ann_table$azimuth_label <- az_labels[ann_table$cluster]
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

if (length(marker_score_info) > 0 && !is.null(marker_score_info$marker_based_labels)) {
  ann_table$marker_based_label <- unlist(marker_score_info$marker_based_labels[ann_table$cluster])
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
# for cell types listed in annotation.marker_override — those are cell types the
# Azimuth reference cannot represent (it has no class for them), so a
# marker-based call for them wins over Azimuth.
override_types <- cfg$annotation$marker_override %||% character(0)
pick_label <- function(az, mk, sr) {
  if (!is.na(mk) && mk %in% override_types) return(mk)
  if (!is.na(az)) return(az)
  if (!is.na(mk)) return(mk)
  if (!is.na(sr)) return(sr)
  "Unknown"
}
ann_table$suggested_label <- mapply(pick_label, ann_table$azimuth_label,
                                    ann_table$marker_based_label, ann_table$singler_label,
                                    USE.NAMES = FALSE)

cluster_counts <- table(as.character(Idents(obj)))
ann_table$n_cells <- as.integer(cluster_counts[ann_table$cluster])

ann_csv_path <- file.path(outdir, "cluster_annotations.csv")
write.csv(ann_table, ann_csv_path, row.names = FALSE)
log_msg(sprintf("Wrote cluster annotations to %s", ann_csv_path))

# --- Add suggested labels to object metadata ---
label_map <- setNames(ann_table$suggested_label, ann_table$cluster)
obj$cell_type <- unname(label_map[as.character(Idents(obj))])

# --- Save final object ---
rds_path <- file.path(outdir, "seurat_final.rds")
log_msg(sprintf("Saving final annotated object to %s", rds_path))
saveRDS(obj, rds_path)

# --- Cell type composition ---
composition_overall <- as.list(sort(table(obj$cell_type), decreasing = TRUE))
composition_overall <- lapply(composition_overall, as.integer)

composition_per_sample <- list()
for (sid in unique(obj$sample_id)) {
  ct <- table(obj$cell_type[obj$sample_id == sid])
  composition_per_sample[[sid]] <- as.list(setNames(as.integer(ct), names(ct)))
}

# --- JSON summary ---
summary_data <- c(
  summary_header("05_annotate", cfg),
  list(
    cells = ncol(obj),
    n_clusters = length(clusters),
    methods = annotation_methods,
    markers = markers_info,
    marker_scoring = marker_score_info,
    cluster_annotations = lapply(seq_len(nrow(ann_table)), function(i) {
      as.list(ann_table[i, ])
    }),
    cell_type_composition = composition_overall,
    per_sample_composition = composition_per_sample,
    outputs = list(
      rds = rds_path,
      cluster_annotations_csv = ann_csv_path,
      cluster_markers_csv = if (exists("markers_csv_path")) markers_csv_path else NULL
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
  cfg_levels <- ann_cfg$azimuth_levels
  if (!is.null(cfg_levels)) {
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

  # Dot plot of marker genes (from config) x clusters
  if (!is.null(ann_cfg$marker_genes)) {
    dot_genes <- c()
    for (compartment in names(ann_cfg$marker_genes)) {
      for (cell_type in names(ann_cfg$marker_genes[[compartment]])) {
        genes <- ann_cfg$marker_genes[[compartment]][[cell_type]]
        present <- intersect(genes, rownames(obj))
        dot_genes <- c(dot_genes, head(present, 2))
      }
    }
    dot_genes <- unique(dot_genes)

    if (length(dot_genes) > 0) {
      possible_cols <- grep("_snn_res\\.", colnames(obj@meta.data), value = TRUE)
      if (length(possible_cols) > 0) Idents(obj) <- obj@meta.data[[possible_cols[1]]]

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
        possible_cols <- grep("_snn_res\\.", colnames(obj@meta.data), value = TRUE)
        if (length(possible_cols) > 0) {
          Idents(obj) <- factor(obj@meta.data[[possible_cols[1]]])
        }
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
