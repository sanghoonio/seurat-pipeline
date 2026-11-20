#!/usr/bin/env bash
set -euo pipefail

# Build the pipeline report: every run of every step (defined in the config or
# found on disk) with its status, written to <output_dir>/pipeline_report.json.
#
# Steps are run one at a time, each for a chosen run label:
#   Rscript scripts/02_ambient.R --config config.yaml --label decontx
#
# Usage:
#   ./run.sh [config.yaml] --report-only

CONFIG="config.yaml"
REPORT_ONLY=0
for arg in "$@"; do
  case "$arg" in
    --report-only) REPORT_ONLY=1 ;;
    -*) echo "ERROR: unknown option: $arg (usage: ./run.sh [config.yaml] --report-only)" >&2; exit 1 ;;
    *) CONFIG="$arg" ;;
  esac
done

if [ "$REPORT_ONLY" -ne 1 ]; then
  echo "ERROR: run.sh only builds the report (--report-only). Run steps one at a time:" >&2
  echo "  Rscript scripts/<step>.R --config $CONFIG [--label <label>] [--no-plots]" >&2
  exit 1
fi
if [ ! -f "$CONFIG" ]; then
  echo "ERROR: config file not found: $CONFIG" >&2
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
Rscript "$SCRIPT_DIR/scripts/report.R" --config "$CONFIG"
