# Custom Seurat reference

Builds a reference for step 05's `annotation.custom_reference` from an atlas `.h5ad`
(for example a CELLxGENE download). One-time and run by hand; not a pipeline step.

1. Export a label-stratified subsample with raw counts and gene symbols, using the
   bundled Python env:

   ```bash
   scripts/scvi/.venv/bin/python scripts/refs/export_h5ad.py \
     --h5ad atlas.h5ad --out-dir ref_export \
     --label-columns cell_type --keep-obs donor_id,assay \
     --max-per-label 2000 --max-cells 100000
   ```

   It refuses non-integer counts (a SCT reference needs raw counts), sums duplicated
   gene symbols, and writes `export_summary.json` with cells per label.

2. Build the reference (SCTransform v2 + PCA). Large atlases need a big-memory job:

   ```bash
   Rscript scripts/refs/build_seurat_reference.R --dir ref_export --out cellref.rds --dims 50
   ```

3. Point the config at it:

   ```yaml
   annotation:
     custom_reference:
       path: "/path/to/cellref.rds"
       label_columns: ["cell_type"]
       dims: 50
   ```
