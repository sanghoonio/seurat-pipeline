#!/usr/bin/env Rscript
# Build a Seurat reference for annotation.custom_reference (step 05) from the
# output of export_h5ad.py: SCTransform (v2) + PCA, saved as one RDS. Step 05
# maps query cells onto it with FindTransferAnchors / TransferData.
#
# Usage: Rscript scripts/refs/build_seurat_reference.R --dir <export dir> --out <ref.rds> [--dims 50] [--nfeatures 3000]
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

a <- commandArgs(trailingOnly = TRUE)
arg <- function(flag, default = NULL) { i <- match(flag, a); if (is.na(i)) default else a[i + 1] }
dir <- arg("--dir"); out <- arg("--out")
dims <- as.integer(arg("--dims", 50)); nfeat <- as.integer(arg("--nfeatures", 3000))
if (is.null(dir) || is.null(out)) stop("usage: build_seurat_reference.R --dir <export dir> --out <ref.rds> [--dims 50] [--nfeatures 3000]")

counts <- as(readMM(file.path(dir, "counts.mtx")), "CsparseMatrix")
rownames(counts) <- readLines(file.path(dir, "genes.txt"))
colnames(counts) <- readLines(file.path(dir, "barcodes.txt"))
obs <- read.csv(file.path(dir, "obs.csv"), check.names = FALSE, stringsAsFactors = FALSE)
rownames(obs) <- obs$barcode
obs$barcode <- NULL
message(sprintf("reference: %d genes x %d cells", nrow(counts), ncol(counts)))

ref <- CreateSeuratObject(counts = counts, meta.data = obs[colnames(counts), , drop = FALSE], min.cells = 0, min.features = 0)
ref <- SCTransform(ref, vst.flavor = "v2", variable.features.n = nfeat, verbose = FALSE)
ref <- RunPCA(ref, npcs = dims, verbose = FALSE, seed.use = 42)
saveRDS(ref, out)

jsonlite::write_json(list(
  built = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), export_dir = normalizePath(dir), out = normalizePath(out),
  n_cells = ncol(ref), n_genes = nrow(ref), dims = dims, variable_features = nfeat,
  label_columns = setdiff(colnames(obs), character(0)),
  versions = list(Seurat = as.character(packageVersion("Seurat")), R = paste(R.version$major, R.version$minor, sep = "."))
), sub("\\.rds$", "_build.json", out), auto_unbox = TRUE, pretty = TRUE)
message(sprintf("saved %s", out))
