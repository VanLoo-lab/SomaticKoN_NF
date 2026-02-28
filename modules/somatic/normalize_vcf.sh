#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
normalize_vcf.sh

Normalize and standardize a VCF using bcftools:
  - bgzip + tabix if needed
  - left-align/normalize against reference
  - (optional) split multiallelics
  - (optional) emit "pre-split" normalized VCF as well

Required:
  --in-vcf   PATH   (vcf or vcf.gz)
  --ref      PATH   (reference fasta)

Optional:
  --out      PATH   (output split-normalized vcf.gz)
                     default: <input>.norm.vcf.gz
  --out-presplit PATH  (output pre-split normalized vcf.gz)
                     default: <input>.norm.presplit.vcf.gz
  --emit-presplit true|false   default: true
  --threads  INT    default: 4
  --regions  PATH   (BED or region list) optional speed-up
  --keep-standard-contigs true|false   default: false
  --contig-regex  REGEX  used if keep-standard-contigs=true
                  default (hg38 chr): '^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y|M)$'
  --workdir  PATH   directory for temp files (default: same as output)

Outputs:
  <out>.tbi tabix index
  <out-presplit>.tbi tabix index (if emit-presplit=true)
EOF
}

IN_VCF=""
REF=""
OUT=""
OUT_PRESPLIT=""
EMIT_PRESPLIT="true"
THREADS=4
REGIONS=""
KEEP_STD="false"
CONTIG_REGEX='^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y|M)$'
WORKDIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --in-vcf) IN_VCF="$2"; shift 2;;
    --ref) REF="$2"; shift 2;;
    --out) OUT="$2"; shift 2;;
    --out-presplit) OUT_PRESPLIT="$2"; shift 2;;
    --emit-presplit) EMIT_PRESPLIT="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --regions) REGIONS="$2"; shift 2;;
    --keep-standard-contigs) KEEP_STD="$2"; shift 2;;
    --contig-regex) CONTIG_REGEX="$2"; shift 2;;
    --workdir) WORKDIR="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1" >&2; usage; exit 2;;
  esac
done

[[ -n "$IN_VCF" && -n "$REF" ]] || { echo "ERROR: --in-vcf and --ref are required" >&2; usage; exit 2; }
[[ -f "$IN_VCF" ]] || { echo "ERROR: input not found: $IN_VCF" >&2; exit 2; }
[[ -f "$REF" ]] || { echo "ERROR: reference not found: $REF" >&2; exit 2; }
[[ -f "${REF}.fai" ]] || { echo "ERROR: missing fasta index: ${REF}.fai" >&2; exit 2; }

command -v bcftools >/dev/null 2>&1 || { echo "ERROR: bcftools not found in PATH" >&2; exit 127; }
command -v tabix    >/dev/null 2>&1 || { echo "ERROR: tabix not found in PATH" >&2; exit 127; }

# default output names
base="$(basename "$IN_VCF")"
base="${base%.vcf.gz}"
base="${base%.vcf}"

if [[ -z "$OUT" ]]; then
  OUT="${base}.norm.vcf.gz"
fi
if [[ -z "$OUT_PRESPLIT" ]]; then
  OUT_PRESPLIT="${base}.norm.presplit.vcf.gz"
fi

# Determine workdir - use output directory if not specified
if [[ -z "$WORKDIR" ]]; then
  WORKDIR="$(dirname "$OUT")"
fi

# Create workdir if it doesn't exist
mkdir -p "$WORKDIR"

# Create temp directory under workdir (same filesystem as output)
tmpdir="${WORKDIR}/.normalize_tmp.$$"
mkdir -p "$tmpdir"

# Cleanup function
cleanup() {
  if [[ -d "$tmpdir" ]]; then
    rm -rf "$tmpdir"
    echo "Cleaned up temp directory: $tmpdir"
  fi
}

# Trap to cleanup on exit (success or failure)
trap cleanup EXIT

echo "Using temp directory: $tmpdir"

# ensure bgzip + index for input
INPUT_GZ="$IN_VCF"
if [[ "$IN_VCF" != *.vcf.gz ]]; then
  INPUT_GZ="${tmpdir}/input.vcf.gz"
  bcftools view --threads "$THREADS" "$IN_VCF" -Oz -o "$INPUT_GZ"
  tabix -f -p vcf "$INPUT_GZ"
else
  if [[ ! -f "${IN_VCF}.tbi" && ! -f "${IN_VCF}.csi" ]]; then
    tabix -f -p vcf "$IN_VCF"
  fi
fi

# Optional region restriction
REGION_ARGS=()
if [[ -n "$REGIONS" ]]; then
  [[ -f "$REGIONS" ]] || { echo "ERROR: regions file not found: $REGIONS" >&2; exit 2; }
  REGION_ARGS=( -R "$REGIONS" )
fi

filter_std_contigs_or_copy() {
  local in_vcf="$1"
  local out_vcf="$2"

  # Ensure output directory exists
  mkdir -p "$(dirname "$out_vcf")"

  if [[ "$KEEP_STD" == "true" ]]; then
    bcftools view \
      --threads "$THREADS" \
      -i "REGEX(%CHROM,\"$CONTIG_REGEX\")" \
      -Oz -o "$out_vcf" \
      "$in_vcf"
    tabix -f -p vcf "$out_vcf"
  else
    # Use cp instead of mv to avoid cross-filesystem issues
    cp "$in_vcf" "$out_vcf"
    if [[ -f "${in_vcf}.tbi" ]]; then
      cp "${in_vcf}.tbi" "${out_vcf}.tbi"
    else
      tabix -f -p vcf "$out_vcf"
    fi
  fi
}

# --------------------------
# 1) PRE-SPLIT normalization
#    (normalize/left-align vs reference, do NOT split multiallelics)
# --------------------------
presplit_tmp="${tmpdir}/presplit.norm.vcf.gz"

echo "Step 1: Pre-split normalization..."
bcftools norm \
  "${REGION_ARGS[@]}" \
  -f "$REF" \
  --check-ref w \
  --threads "$THREADS" \
  "$INPUT_GZ" \
| bcftools view -Oz --threads "$THREADS" -o "$presplit_tmp"

tabix -f -p vcf "$presplit_tmp"

if [[ "$EMIT_PRESPLIT" == "true" ]]; then
  filter_std_contigs_or_copy "$presplit_tmp" "$OUT_PRESPLIT"
  echo "Wrote pre-split: $OUT_PRESPLIT"
  echo "Index: ${OUT_PRESPLIT}.tbi"
fi

# Use presplit_tmp for step 2 (it still exists in tmpdir)
presplit_for_split="$presplit_tmp"

# --------------------------
# 2) SPLIT multiallelics on top of pre-split normalized
# --------------------------
split_tmp="${tmpdir}/split.norm.vcf.gz"

echo "Step 2: Split multiallelics..."
bcftools norm \
  "${REGION_ARGS[@]}" \
  -m -any \
  -f "$REF" \
  --check-ref w \
  --threads "$THREADS" \
  "$presplit_for_split" \
| bcftools view -Oz --threads "$THREADS" -o "$split_tmp"

tabix -f -p vcf "$split_tmp"

filter_std_contigs_or_copy "$split_tmp" "$OUT"
echo "Wrote split: $OUT"
echo "Index: ${OUT}.tbi"

echo "normalize_vcf.sh completed successfully"