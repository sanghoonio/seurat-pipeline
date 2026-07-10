#!/usr/bin/env bash
set -euo pipefail

# Pipeline runner. Executes each step in order, captures stdout/stderr per step,
# then collects every *_summary.json into a single pipeline_report.json that
# downstream tools (or LLMs) can read.
#
# Usage:
#   ./run.sh                # uses ./config.yaml
#   ./run.sh path/to/cfg.yaml

CONFIG="${1:-config.yaml}"

if [ ! -f "$CONFIG" ]; then
  echo "ERROR: config file not found: $CONFIG" >&2
  echo "Copy config.template.yaml to config.yaml and edit it for your data." >&2
  exit 1
fi

# Resolve repo dir (so the pipeline can be invoked from any working directory)
# and the absolute config path.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="$(cd "$(dirname "$CONFIG")" && pwd)/$(basename "$CONFIG")"

# Honor the config's paths.output_dir / paths.log_dir (resolved against
# paths.base_dir) so analysis artifacts land where the config says, not in the
# repo. Falls back to ./outputs and ./logs if python/yaml is unavailable.
read_cfg_path() {
  python3 - "$CONFIG" "$1" "$2" <<'PY' 2>/dev/null
import sys, os, yaml
cfg_path, key, default = sys.argv[1], sys.argv[2], sys.argv[3]
with open(cfg_path) as f:
    cfg = yaml.safe_load(f) or {}
paths = cfg.get('paths', {}) or {}
base = paths.get('base_dir', '.') or '.'
val = paths.get(key, default) or default
if not os.path.isabs(val):
    val = os.path.normpath(os.path.join(base, val))
print(val)
PY
}
OUTDIR="$(read_cfg_path output_dir outputs || true)"; OUTDIR="${OUTDIR:-outputs}"
LOGDIR="$(read_cfg_path log_dir logs || true)";       LOGDIR="${LOGDIR:-logs}"

mkdir -p "$LOGDIR" "$OUTDIR"

RUN_TS=$(date +%Y%m%d_%H%M%S)

echo "=== Seurat Pipeline Run: $RUN_TS ==="
echo "Config:  $CONFIG"
echo "Outputs: $OUTDIR"
echo "Logs:    $LOGDIR"
echo ""

STEPS=(
  "00_cell_calling"
  "01_load_qc_metrics"
  "02_soupx"
  "03_filter"
  "04_normalize_cluster"
  "05_annotate"
)

for step in "${STEPS[@]}"; do
  echo "[$(date +%H:%M:%S)] Starting $step ..."

  Rscript "${SCRIPT_DIR}/scripts/${step}.R" \
    --config "$CONFIG" \
    > "${LOGDIR}/${step}_${RUN_TS}.stdout.log" \
    2> "${LOGDIR}/${step}_${RUN_TS}.stderr.log"

  EXIT_CODE=$?

  if [ $EXIT_CODE -ne 0 ]; then
    echo "[$(date +%H:%M:%S)] FAILED: $step (exit code $EXIT_CODE)"
    echo "--- stderr tail ---"
    tail -20 "${LOGDIR}/${step}_${RUN_TS}.stderr.log"
    exit $EXIT_CODE
  fi

  echo "[$(date +%H:%M:%S)] Completed $step"
done

echo ""

# --- Collect all JSON summaries into a single pipeline report ---
echo "Collecting pipeline summary..."
python3 -c "
import json, glob, os

summaries = {}
for f in sorted(glob.glob('${OUTDIR}/*_summary.json') + glob.glob('${OUTDIR}/*summary.json')):
    name = os.path.basename(f).replace('_summary.json', '').replace('.json', '')
    with open(f) as fh:
        summaries[name] = json.load(fh)

# Add config for full reproducibility context
import yaml
with open('${CONFIG}') as fh:
    summaries['config'] = yaml.safe_load(fh)

with open('${OUTDIR}/pipeline_report.json', 'w') as fh:
    json.dump(summaries, fh, indent=2, default=str)

print(f'Wrote pipeline_report.json ({len(summaries)} sections)')
" 2>/dev/null || echo "(pipeline_report.json collection skipped — python3/pyyaml not available)"

echo "=== Pipeline complete ==="
echo "Outputs: $OUTDIR/"
echo "Logs:    $LOGDIR/"
echo ""
echo "Key files for downstream / AI consumption:"
echo "  $OUTDIR/pipeline_report.json    <- merged summary of all steps + config"
echo "  $OUTDIR/qc_metrics.csv          <- per-cell QC metrics (pre-filter)"
echo "  $OUTDIR/filter_summary.json     <- filtering decisions and cell counts"
echo "  $OUTDIR/cluster_summary.json    <- clustering results and composition"
echo "  $OUTDIR/cluster_annotations.csv <- per-cluster labels (Azimuth/SingleR/markers)"
echo "  $OUTDIR/umap_coordinates.csv    <- UMAP + cluster + sample per cell"
echo "  $OUTDIR/seurat_final.rds        <- final annotated Seurat object"
