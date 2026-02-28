#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Usage:
  mutect2_scatter.sh \
    --sample-id S1 --normal-id N1 \
    --tumor-bam T.bam --normal-bam N.bam \
    --ref ref.fa \
    --interval chr1|intervals.list|region.bed \
    --germline-resource gnomad.vcf.gz \
    --pon pon.vcf.gz \
    --out-prefix S1.chr1 \
    [--extra-args "..."]

Outputs (in cwd):
  <out-prefix>.unfiltered.vcf.gz
  <out-prefix>.unfiltered.vcf.gz.stats
  <out-prefix>.f1r2.tar.gz
  <out-prefix>.log
EOF
}

# defaults
EXTRA_ARGS=""

# args
SAMPLE_ID=""
NORMAL_ID=""
TBAM=""
NBAM=""
REF=""
INTERVAL=""
GERMLINE=""
PON=""
OUT_PREFIX=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample-id) SAMPLE_ID="$2"; shift 2;;
    --normal-id) NORMAL_ID="$2"; shift 2;;
    --tumor-bam) TBAM="$2"; shift 2;;
    --normal-bam) NBAM="$2"; shift 2;;
    --ref) REF="$2"; shift 2;;
    --interval) INTERVAL="$2"; shift 2;;
    --germline-resource) GERMLINE="$2"; shift 2;;
    --pon) PON="$2"; shift 2;;
    --out-prefix) OUT_PREFIX="$2"; shift 2;;
    --extra-args) EXTRA_ARGS="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 2;;
  esac
done

# validate
for f in "$TBAM" "$NBAM" "$REF" "${REF}.fai"; do
  [[ -f "$f" ]] || { echo "Missing required file: $f" >&2; exit 10; }
done
[[ -n "$SAMPLE_ID" && -n "$NORMAL_ID" && -n "$INTERVAL" && -n "$GERMLINE" && -n "$PON" && -n "$OUT_PREFIX" ]] || {
  echo "Missing required arguments" >&2; usage; exit 11;
}

LOG="${OUT_PREFIX}.log"
VCF="${OUT_PREFIX}.unfiltered.vcf.gz"
F1R2="${OUT_PREFIX}.f1r2.tar.gz"

echo "[Mutect2 scatter] sample=$SAMPLE_ID normal=$NORMAL_ID interval=$INTERVAL" | tee "$LOG"

# Run Mutect2
gatk Mutect2 \
  -R "$REF" \
  -I "$TBAM" \
  -I "$NBAM" \
  -normal "$NORMAL_ID" \
  -L "$INTERVAL" \
  --germline-resource "$GERMLINE" \
  --panel-of-normals "$PON" \
  --f1r2-tar-gz "$F1R2" \
  $EXTRA_ARGS \
  -O "$VCF" >>"$LOG" 2>&1

# Ensure stats exists (GATK writes <vcf>.stats)
[[ -f "${VCF}.stats" ]] || {
  echo "Expected stats file missing: ${VCF}.stats" | tee -a "$LOG"
  exit 12
}

echo "[Mutect2 scatter] done: $VCF" | tee -a "$LOG"