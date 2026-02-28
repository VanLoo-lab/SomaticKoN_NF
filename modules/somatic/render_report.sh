#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
render_report.sh --sample-id ID --outdir DIR   [--metrics-dir DIR] [--plots-dir DIR] [--report-qmd FILE] [--format png|pdf|both]
EOF
}

SAMPLE=""
OUTDIR=""
METRICS_DIR=""
PLOTS_DIR=""
QMD=""
FORMAT="png"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample-id) SAMPLE="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --metrics-dir) METRICS_DIR="$2"; shift 2;;
    --plots-dir) PLOTS_DIR="$2"; shift 2;;
    --report-qmd) QMD="$2"; shift 2;;
    --format) FORMAT="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1" >&2; usage; exit 2;;
  esac
done

[[ -n "$SAMPLE" && -n "$OUTDIR" ]] || { echo "ERROR: --sample-id and --outdir required" >&2; usage; exit 2; }

METRICS_DIR="${METRICS_DIR:-$OUTDIR/metrics}"
PLOTS_DIR="${PLOTS_DIR:-$OUTDIR/plots}"
QMD="${QMD:-$(dirname "$0")/somatic_report.qmd}"

mkdir -p "$PLOTS_DIR"

Rscript "$(dirname "$0")/plot_somatic_qc.R" "$METRICS_DIR" "$PLOTS_DIR" --sample-id "$SAMPLE" --format "$FORMAT"

quarto render "$QMD" -P sample_id:"$SAMPLE" -P metrics_dir:"$METRICS_DIR" -P plots_dir:"$PLOTS_DIR" --output-dir "$OUTDIR"

echo "Report written to: $OUTDIR"
