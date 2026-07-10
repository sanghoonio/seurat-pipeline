#!/usr/bin/env Rscript
# 00_cell_calling.R
# Call cells from the RAW matrix with DropletUtils::emptyDrops to recover
# low-RNA cells (e.g. neutrophils) that Cell Ranger's filtered set discards.
# Runs once per sample, writes the called-cell barcodes (tagged by source) plus
# a rich diagnostic summary; steps 01/02 read the chosen set (cell_calling.use).
#
# No-op when cell_calling.method != "emptydrops" (writes a "skipped" summary),
# so default pipeline behavior is unchanged.
#
# Writes: outputs/cell_calling/<id>_barcodes.csv,
#         outputs/cell_calling/<id>_barcode_ranks.csv,
#         cell_calling_summary.json

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
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

outdir <- output_dir(cfg)

cc <- cfg$cell_calling %||% list(method = "filtered")
method <- cc$method %||% "filtered"

# --- No-op path: cell-calling disabled ---
if (method != "emptydrops") {
  log_msg("cell_calling.method != 'emptydrops' — skipping (downstream uses Cell Ranger filtered set)")
  write_json_summary(
    c(summary_header("00_cell_calling", cfg),
      list(skipped = TRUE, reason = "cell_calling disabled", method = method, status = "success")),
    file.path(outdir, "cell_calling_summary.json")
  )
  log_msg("Step 00 complete (skipped)")
  quit(status = 0, save = "no")
}

# --- emptyDrops parameters ---
use_set      <- cc$use %||% "union"
fdr          <- cc$fdr %||% 0.01
lower        <- cc$lower %||% 100
niters       <- cc$niters %||% 10000
retain       <- cc$retain          # NULL -> emptyDrops auto (knee)
alpha        <- cc$alpha           # NULL -> emptyDrops auto (MLE)
test_ambient <- isTRUE(cc$test_ambient)
mt_pattern   <- cfg$qc$mt_pattern %||% "^MT-"

cc_dir <- file.path(outdir, "cell_calling")
dir.create(cc_dir, recursive = TRUE, showWarnings = FALSE)

log_msg(sprintf("emptyDrops cell-calling: fdr=%.3g lower=%g niters=%g use=%s",
                fdr, lower, niters, use_set))

per_sample <- list()

for (s in cfg$samples) {
  sid <- s$id
  raw_path  <- resolve_path(sample_raw_path(s), cfg)
  filt_path <- resolve_path(sample_filtered_path(s), cfg)
  log_msg(sprintf("  [%s] raw: %s", sid, raw_path))

  raw_counts <- read_10x_any(raw_path)
  filt_bc    <- colnames(read_10x_any(filt_path))

  # --- emptyDrops ---
  cc_res <- call_cells_emptydrops(raw_counts, fdr = fdr, lower = lower,
                                  niters = niters, retain = retain, alpha = alpha,
                                  test_ambient = test_ambient, seed = cfg$project$seed)
  ed     <- cc_res$ed
  ed_bc  <- cc_res$emptydrops_barcodes

  # --- Set operations vs Cell Ranger filtered ---
  filt_in_raw     <- intersect(filt_bc, colnames(raw_counts))
  union_bc        <- union(ed_bc, filt_in_raw)
  intersection_bc <- intersect(ed_bc, filt_in_raw)
  emptydrops_only <- setdiff(ed_bc, filt_in_raw)
  filtered_only   <- setdiff(filt_in_raw, ed_bc)

  # --- Write the called-barcode table (downstream reads the `use` subset) ---
  bc_df <- data.frame(
    barcode       = union_bc,
    in_filtered   = union_bc %in% filt_in_raw,
    in_emptydrops = union_bc %in% ed_bc,
    stringsAsFactors = FALSE
  )
  write.csv(bc_df, file.path(cc_dir, paste0(sid, "_barcodes.csv")), row.names = FALSE)

  # --- Barcode-rank-plot data (curve + knee/inflection) ---
  br <- DropletUtils::barcodeRanks(raw_counts, lower = lower)
  knee       <- S4Vectors::metadata(br)$knee
  inflection <- S4Vectors::metadata(br)$inflection
  # Record per-barcode rank/total/FDR for the informative (tested) region only,
  # i.e. barcodes with total >= lower (the millions of sub-lower empties are a
  # flat tail and would bloat the file). Knee/inflection are reported in JSON.
  keep_rank <- which(br$total >= lower)
  rk_bc <- rownames(br)[keep_rank]
  ranks_df <- data.frame(
    barcode       = rk_bc,
    total         = br$total[keep_rank],
    rank          = br$rank[keep_rank],
    FDR           = ed[rk_bc, "FDR"],
    in_filtered   = rk_bc %in% filt_in_raw,
    in_emptydrops = rk_bc %in% ed_bc,
    stringsAsFactors = FALSE
  )
  ranks_df <- ranks_df[order(ranks_df$rank), ]
  write.csv(ranks_df, file.path(cc_dir, paste0(sid, "_barcode_ranks.csv")), row.names = FALSE)

  # --- Stats on the recovered (emptydrops_only) fraction: are they real cells? ---
  ed_only_stats <- list(n = length(emptydrops_only))
  if (length(emptydrops_only) > 0) {
    m  <- raw_counts[, emptydrops_only, drop = FALSE]
    nc <- Matrix::colSums(m)
    nf <- Matrix::colSums(m > 0)
    mt_rows <- grepl(mt_pattern, rownames(m))
    mt_pct  <- if (any(mt_rows)) 100 * Matrix::colSums(m[mt_rows, , drop = FALSE]) / pmax(nc, 1) else rep(0, ncol(m))
    ed_only_stats$median_nCount   <- as.numeric(median(nc))
    ed_only_stats$median_nFeature <- as.numeric(median(nf))
    ed_only_stats$median_pct_mt   <- round(as.numeric(median(mt_pct)), 3)
  }

  # --- Diagnostic: does `lower` itself swallow cell-like barcodes? ---
  sub_lower_celllike <- NA_integer_
  if (test_ambient) {
    sub_lower <- !is.na(ed$FDR) & ed$Total <= lower & ed$FDR <= fdr
    sub_lower_celllike <- sum(sub_lower)
  }

  per_sample[[sid]] <- list(
    n_filtered        = length(filt_bc),
    n_filtered_in_raw = length(filt_in_raw),
    n_emptydrops      = length(ed_bc),
    n_union           = length(union_bc),
    n_intersection    = length(intersection_bc),
    n_emptydrops_only = length(emptydrops_only),
    n_filtered_only   = length(filtered_only),
    n_used            = switch(use_set,
                               union = length(union_bc),
                               emptydrops = length(ed_bc),
                               filtered = length(filt_in_raw)),
    emptydrops_only_stats = ed_only_stats,
    barcode_rank = list(knee = as.numeric(knee), inflection = as.numeric(inflection)),
    sub_lower_celllike = sub_lower_celllike
  )

  log_msg(sprintf("  [%s] filtered=%d emptyDrops=%d union=%d (+%d recovered, -%d filtered-only)",
                  sid, length(filt_bc), length(ed_bc), length(union_bc),
                  length(emptydrops_only), length(filtered_only)))

  # --- Optional barcode-rank plot ---
  if (args$make_plots) {
    p <- ggplot(ranks_df, aes(x = rank, y = total)) +
      geom_line(color = "grey50") +
      geom_point(aes(color = in_emptydrops), size = 0.4, alpha = 0.5) +
      geom_hline(yintercept = c(knee, inflection), linetype = "dashed",
                 color = c("dodgerblue", "forestgreen")) +
      scale_x_log10() + scale_y_log10() +
      scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick"),
                         name = "emptyDrops cell") +
      labs(title = sprintf("%s — barcode rank (knee=blue, inflection=green)", sid),
           x = "Barcode rank", y = "Total UMI") +
      theme_bw(base_size = 12)
    save_plot(p, paste0("00_barcode_rank_", sid, ".png"), cfg, width = 7, height = 5)
  }

  rm(raw_counts, ed, br); gc()
}

# --- JSON summary ---
summary_data <- c(
  summary_header("00_cell_calling", cfg),
  list(
    method = method,
    use = use_set,
    params = list(fdr = fdr, lower = lower, niters = niters,
                  retain = retain %||% "auto", alpha = alpha %||% "auto",
                  test_ambient = test_ambient),
    per_sample = per_sample,
    outputs = list(barcodes_dir = cc_dir),
    status = "success"
  )
)
write_json_summary(summary_data, file.path(outdir, "cell_calling_summary.json"))

log_msg("Step 00 complete")
