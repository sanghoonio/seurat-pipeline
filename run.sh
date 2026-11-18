#!/usr/bin/env bash
set -euo pipefail

# Pipeline runner. Executes each step in order, captures stdout/stderr per step,
# then collects every *_summary.json into a single pipeline_report.json that
# downstream tools (or LLMs) can read.
#
# Usage:
#   ./run.sh                          # uses ./config.yaml
#   ./run.sh path/to/cfg.yaml
#   ./run.sh path/to/cfg.yaml --report-only   # just rebuild pipeline_report.json
#                                             # (e.g. after running steps one by one)

CONFIG="config.yaml"
REPORT_ONLY=0
for arg in "$@"; do
  case "$arg" in
    --report-only) REPORT_ONLY=1 ;;
    -*) echo "ERROR: unknown option: $arg (usage: ./run.sh [config.yaml] [--report-only])" >&2; exit 1 ;;
    *) CONFIG="$arg" ;;
  esac
done

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

# Step -> summary file. The report reads exactly these (a leftover file from an
# older pipeline version, e.g. soupx_summary.json, is ignored).
STEPS=(
  "00_cell_calling:cell_calling_summary.json"
  "01_load_qc_metrics:qc_summary.json"
  "02_ambient:ambient_summary.json"
  "03_filter:filter_summary.json"
  "04_normalize_cluster:cluster_summary.json"
  "05_annotate:annotation_summary.json"
)

# --- Collect all step summaries into a single pipeline report ---
build_report() {
  echo "Collecting pipeline summary..."
  python3 - "$OUTDIR" "$CONFIG" "${STEPS[@]}" <<'PY' || echo "WARNING: pipeline_report.json not written (see error above)" >&2
import json, os, sys, yaml
outdir, config, steps = sys.argv[1], sys.argv[2], sys.argv[3:]
summaries, index = {}, {}
for entry in steps:
    step, fname = entry.split(':', 1)
    section = fname[:-len('_summary.json')]
    path = os.path.join(outdir, fname)
    if not os.path.exists(path):
        index[section] = {'step': step, 'status': 'missing'}
        continue
    with open(path) as fh:
        summaries[section] = s = json.load(fh)
    # Per-section status/time so stale, failed or missing steps are visible at a glance
    index[section] = {'step': step, 'status': s.get('status'), 'timestamp': s.get('timestamp'),
                      'pipeline_commit': s.get('pipeline_commit'),
                      'n_warnings': len(s.get('warnings') or [])}
summaries['_index'] = index
with open(config) as fh:
    summaries['config'] = yaml.safe_load(fh)
with open(os.path.join(outdir, 'pipeline_report.json'), 'w') as fh:
    json.dump(summaries, fh, indent=2, default=str)
print(f"Wrote pipeline_report.json ({len(summaries) - 2} of {len(steps)} step summaries present)")
PY
}

if [ "$REPORT_ONLY" = 1 ]; then
  build_report
  exit 0
fi

echo "=== Seurat Pipeline Run: $RUN_TS ==="
echo "Config:  $CONFIG"
echo "Outputs: $OUTDIR"
echo "Logs:    $LOGDIR"
echo ""

# A full run starts clean: remove every step summary first, so none from an
# earlier run can survive into this run's report.
for entry in "${STEPS[@]}"; do rm -f "${OUTDIR}/${entry#*:}"; done

for entry in "${STEPS[@]}"; do
  step="${entry%%:*}"
  echo "[$(date +%H:%M:%S)] Starting $step ..."

  EXIT_CODE=0
  Rscript "${SCRIPT_DIR}/scripts/${step}.R" \
    --config "$CONFIG" \
    > "${LOGDIR}/${step}_${RUN_TS}.stdout.log" \
    2> "${LOGDIR}/${step}_${RUN_TS}.stderr.log" || EXIT_CODE=$?

  if [ $EXIT_CODE -ne 0 ]; then
    echo "[$(date +%H:%M:%S)] FAILED: $step (exit code $EXIT_CODE)"
    echo "--- stderr tail ---"
    tail -20 "${LOGDIR}/${step}_${RUN_TS}.stderr.log"
    build_report
    exit $EXIT_CODE
  fi

  echo "[$(date +%H:%M:%S)] Completed $step"
done

echo ""

build_report

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
