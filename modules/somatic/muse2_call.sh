#!/usr/bin/env bash
#set -euo pipefail

usage() {
  cat <<'EOF'
MuSE2 somatic calling (tumor/normal) - runner script

This script DOES NOT source resolve_refs.sh.
You should source resolve_refs.sh outside and pass resolved paths like:
  --ref "$REF_GENOME" --dbsnp "$DBSNP"

Required:
  --sample-id   STR
  --tumor-bam   PATH
  --normal-bam  PATH
  --ref         PATH to reference FASTA
  --dbsnp       PATH to dbSNP VCF(.gz) used by MuSE sump
  --outdir      PATH output directory

Optional:
  --threads     INT (default: 8)
  --exome       true|false (default: false) - use -E flag for WES data
  --keep-all    true|false (default: false) - keep unfiltered VCF in addition to PASS
  -h|--help

Outputs (in <outdir>/):
  <sample-id>.muse2.pass.vcf.gz       PASS variants only (default output)
  <sample-id>.muse2.pass.vcf.gz.tbi   Index
  
  (if --keep-all true):
  <sample-id>.muse2.all.vcf.gz        All variants (unfiltered)
  <sample-id>.muse2.all.vcf.gz.tbi    Index
EOF
}

THREADS=8
EXOME="false"
KEEP_ALL="false"

SAMPLE_ID=""
TBAM=""
NBAM=""
REF_FASTA=""
DBSNP_VCF=""
OUTDIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample-id) SAMPLE_ID="$2"; shift 2;;
    --tumor-bam) TBAM="$2"; shift 2;;
    --normal-bam) NBAM="$2"; shift 2;;
    --ref) REF_FASTA="$2"; shift 2;;
    --dbsnp) DBSNP_VCF="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --exome) EXOME="$2"; shift 2;;
    --keep-all) KEEP_ALL="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1" >&2; usage; exit 2;;
  esac
done

# validate
[[ -n "$SAMPLE_ID" && -n "$TBAM" && -n "$NBAM" && -n "$REF_FASTA" && -n "$DBSNP_VCF" && -n "$OUTDIR" ]] || {
  echo "ERROR: missing required arguments" >&2
  usage
  exit 2
}

for f in "$TBAM" "$NBAM"; do
  [[ -f "$f" ]] || { echo "ERROR: missing BAM: $f" >&2; exit 10; }
  # optional index check (MuSE needs BAM index)
  if [[ ! -f "${f}.bai" && ! -f "${f%.bam}.bai" ]]; then
    echo "WARNING: BAM index not found for $f (expected .bai next to BAM)" >&2
  fi
done

[[ -f "$REF_FASTA" ]] || { echo "ERROR: reference fasta not found: $REF_FASTA" >&2; exit 11; }
[[ -f "${REF_FASTA}.fai" ]] || echo "WARNING: missing fasta index: ${REF_FASTA}.fai" >&2

[[ -f "$DBSNP_VCF" ]] || { echo "ERROR: dbSNP VCF not found: $DBSNP_VCF" >&2; exit 12; }
if [[ "$DBSNP_VCF" == *.vcf.gz ]]; then
  [[ -f "${DBSNP_VCF}.tbi" ]] || echo "WARNING: missing dbSNP tabix index: ${DBSNP_VCF}.tbi" >&2
fi

command -v bcftools >/dev/null 2>&1 || { echo "ERROR: bcftools not found in PATH" >&2; exit 20; }
command -v tabix >/dev/null 2>&1 || { echo "ERROR: tabix not found in PATH" >&2; exit 20; }

# ---- setup dirs ----
mkdir -p "$OUTDIR"

# Temp prefix for intermediate files
WORK_PREFIX="${OUTDIR}/${SAMPLE_ID}.muse2.tmp"
MUSE_TXT="${WORK_PREFIX}.MuSE.txt"
MUSE_VCF_RAW="${WORK_PREFIX}.vcf"

# Final output files
PASS_VCF="${OUTDIR}/${SAMPLE_ID}.muse2.pass.vcf.gz"
ALL_VCF="${OUTDIR}/${SAMPLE_ID}.muse2.all.vcf.gz"

echo "[MuSE2] sample=${SAMPLE_ID}"
echo "[MuSE2] tumor_bam=${TBAM}"
echo "[MuSE2] normal_bam=${NBAM}"
echo "[MuSE2] ref=${REF_FASTA}"
echo "[MuSE2] dbsnp=${DBSNP_VCF}"
echo "[MuSE2] outdir=${OUTDIR}"
echo "[MuSE2] exome=${EXOME} keep_all=${KEEP_ALL}"

# 1) call -> <WORK_PREFIX>.MuSE.txt
echo "[MuSE2] Running MuSE call..."
MuSE call \
  -f "$REF_FASTA" \
  -O "$WORK_PREFIX" \
  -n "$THREADS" \
  "$TBAM" \
  "$NBAM"

[[ -f "$MUSE_TXT" ]] || { echo "ERROR: expected output not found: $MUSE_TXT" >&2; exit 21; }

# 2) sump -> VCF
echo "[MuSE2] Running MuSE sump..."

SUMP_ARGS=(
  -I "$MUSE_TXT"
  -O "$MUSE_VCF_RAW"
  -n "$THREADS" \
  -D "$DBSNP_VCF"
)

# Add -E for exome/WES data, -G for WGS data
if [[ "$EXOME" == "true" ]]; then
  SUMP_ARGS+=( -E )
else
  SUMP_ARGS+=( -G )
fi

MuSE sump "${SUMP_ARGS[@]}"

[[ -f "$MUSE_VCF_RAW" ]] || { echo "ERROR: expected output not found: $MUSE_VCF_RAW" >&2; exit 22; }

# Count raw variants
RAW_COUNT=$(grep -v "^#" "$MUSE_VCF_RAW" | wc -l)
echo "[MuSE2] Raw variants: $RAW_COUNT"

# 3) Filter PASS only
echo "[MuSE2] Filtering PASS variants..."
bcftools view -f PASS -Oz -o "$PASS_VCF" "$MUSE_VCF_RAW"
tabix -p vcf "$PASS_VCF"

PASS_COUNT=$(bcftools view -H "$PASS_VCF" | wc -l)
echo "[MuSE2] PASS variants: $PASS_COUNT"

# 4) Optionally keep all variants
if [[ "$KEEP_ALL" == "true" ]]; then
  echo "[MuSE2] Keeping unfiltered VCF..."
  bcftools view -Oz -o "$ALL_VCF" "$MUSE_VCF_RAW"
  tabix -p vcf "$ALL_VCF"
  echo "[MuSE2] All variants saved to: $ALL_VCF"
fi

# 5) Clean up intermediate files
echo "[MuSE2] Cleaning up intermediate files..."
rm -f "$MUSE_TXT" "$MUSE_VCF_RAW"

# ---- summary ----
echo ""
echo "[MuSE2] ========== Summary =========="
echo "[MuSE2] Sample: $SAMPLE_ID"
echo "[MuSE2] Raw variants: $RAW_COUNT"
echo "[MuSE2] PASS variants: $PASS_COUNT"
echo "[MuSE2]"
echo "[MuSE2] Output files:"
echo "[MuSE2]   $PASS_VCF"
if [[ "$KEEP_ALL" == "true" ]]; then
  echo "[MuSE2]   $ALL_VCF"
fi
echo "[MuSE2] =============================="
echo "[MuSE2] done"

