#!/usr/bin/env bash
#set -euo pipefail

usage() {
  cat <<'EOF'
Strelka2 somatic variant calling (tumor/normal) - runner script

This script DOES NOT source resolve_refs.sh.
Source resolve_refs.sh outside and pass resolved paths, e.g.:
  --ref "$REF_GENOME"

Required:
  --sample-id    STR
  --tumor-bam    PATH
  --normal-bam   PATH
  --ref          PATH to reference FASTA
  --outdir       PATH output directory (will create <outdir>/<sample-id>/)

Optional:
  --threads      INT  (default: 8)
  --exome        true|false (default: false)
  --call-regions PATH to BED.gz|BED (limits calling; great for chr22 tests)
  --run-mode     local|sge|pbs (default: local)  # usually "local" under Nextflow
  --keep-all     true|false (default: false) - keep all variants, not just PASS
  -h|--help

Outputs (in <outdir>/<sample-id>/):
  <sample-id>.strelka2.snvs.pass.vcf.gz       PASS SNVs only
  <sample-id>.strelka2.indels.pass.vcf.gz     PASS indels only
  <sample-id>.strelka2.pass.vcf.gz            Merged PASS SNVs + indels
  
  (if --keep-all true, also outputs unfiltered versions)
EOF
}

SAMPLE_ID=""
TBAM=""
NBAM=""
REF_FASTA=""
OUTDIR=""

THREADS=8
EXOME="false"
CALL_REGIONS=""
RUN_MODE="local"
KEEP_ALL="false"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample-id) SAMPLE_ID="$2"; shift 2;;
    --tumor-bam) TBAM="$2"; shift 2;;
    --normal-bam) NBAM="$2"; shift 2;;
    --ref) REF_FASTA="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;

    --threads) THREADS="$2"; shift 2;;
    --exome) EXOME="$2"; shift 2;;
    --call-regions) CALL_REGIONS="$2"; shift 2;;
    --run-mode) RUN_MODE="$2"; shift 2;;
    --keep-all) KEEP_ALL="$2"; shift 2;;

    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1" >&2; usage; exit 2;;
  esac
done

# ---- validate ----
[[ -n "$SAMPLE_ID" && -n "$TBAM" && -n "$NBAM" && -n "$REF_FASTA" && -n "$OUTDIR" ]] || {
  echo "ERROR: missing required arguments" >&2
  usage
  exit 2
}

[[ -f "$TBAM" ]] || { echo "ERROR: tumor BAM not found: $TBAM" >&2; exit 10; }
[[ -f "$NBAM" ]] || { echo "ERROR: normal BAM not found: $NBAM" >&2; exit 10; }
[[ -f "$REF_FASTA" ]] || { echo "ERROR: reference FASTA not found: $REF_FASTA" >&2; exit 11; }
[[ -f "${REF_FASTA}.fai" ]] || echo "WARNING: missing FASTA index: ${REF_FASTA}.fai" >&2

if [[ -n "$CALL_REGIONS" ]]; then
  [[ -f "$CALL_REGIONS" ]] || { echo "ERROR: call-regions not found: $CALL_REGIONS" >&2; exit 12; }
fi

#command -v configureStrelkaSomaticWorkflow.py >/dev/null 2>&1 || {
#  echo "ERROR: configureStrelkaSomaticWorkflow.py not found in PATH (strelka2 not activated?)" >&2
#  exit 20
#}

#command -v bcftools >/dev/null 2>&1 || {
#  echo "ERROR: bcftools not found in PATH" >&2
#  exit 20
#}

# ---- setup dirs ----
SAMPLE_OUT="${OUTDIR%/}/${SAMPLE_ID}"
WORKDIR="${SAMPLE_OUT}/strelka_work"
RUN_SCRIPT="${WORKDIR}/runWorkflow.py"

mkdir -p "$SAMPLE_OUT"
mkdir -p "$WORKDIR"

echo "[Strelka2] sample=${SAMPLE_ID}"
echo "[Strelka2] tumor_bam=${TBAM}"
echo "[Strelka2] normal_bam=${NBAM}"
echo "[Strelka2] ref=${REF_FASTA}"
echo "[Strelka2] out=${SAMPLE_OUT}"
echo "[Strelka2] threads=${THREADS} exome=${EXOME} run_mode=${RUN_MODE}"
[[ -n "$CALL_REGIONS" ]] && echo "[Strelka2] call_regions=${CALL_REGIONS}"

# ---- configure ----
CFG_ARGS=(
  --tumorBam "$TBAM"
  --normalBam "$NBAM"
  --referenceFasta "$REF_FASTA"
  --runDir "$WORKDIR"
)

if [[ "$EXOME" == "true" ]]; then
  CFG_ARGS+=( --exome )
fi

if [[ -n "$CALL_REGIONS" ]]; then
  # Strelka accepts a "callRegions" bed.gz (recommended). BED also often works.
  CFG_ARGS+=( --callRegions "$CALL_REGIONS" )
fi

configureStrelkaSomaticWorkflow.py "${CFG_ARGS[@]}"

[[ -f "$RUN_SCRIPT" ]] || { echo "ERROR: Strelka runWorkflow.py not found after configuration" >&2; exit 21; }

# ---- run ----
# Under Nextflow, "local" mode is typical: Strelka does its own parallelism using -j
python "$RUN_SCRIPT" -m "$RUN_MODE" -j "$THREADS"

# ---- check raw outputs ----
SNV_VCF_RAW="${WORKDIR}/results/variants/somatic.snvs.vcf.gz"
INDEL_VCF_RAW="${WORKDIR}/results/variants/somatic.indels.vcf.gz"

[[ -f "$SNV_VCF_RAW" ]] || { echo "ERROR: missing output: $SNV_VCF_RAW" >&2; exit 30; }
[[ -f "$INDEL_VCF_RAW" ]] || { echo "ERROR: missing output: $INDEL_VCF_RAW" >&2; exit 31; }

# ---- filter PASS only and create final outputs ----
echo "[Strelka2] Filtering PASS variants..."

# Output file paths
SNV_PASS="${SAMPLE_OUT}/${SAMPLE_ID}.strelka2.snvs.pass.vcf.gz"
INDEL_PASS="${SAMPLE_OUT}/${SAMPLE_ID}.strelka2.indels.pass.vcf.gz"
MERGED_PASS="${SAMPLE_OUT}/${SAMPLE_ID}.strelka2.pass.vcf.gz"

# Filter SNVs - PASS only
bcftools view -f PASS -Oz -o "$SNV_PASS" "$SNV_VCF_RAW"
bcftools index -t "$SNV_PASS"

SNV_PASS_COUNT=$(bcftools view -H "$SNV_PASS" | wc -l)
echo "[Strelka2] PASS SNVs: $SNV_PASS_COUNT"

# Filter indels - PASS only
bcftools view -f PASS -Oz -o "$INDEL_PASS" "$INDEL_VCF_RAW"
bcftools index -t "$INDEL_PASS"

INDEL_PASS_COUNT=$(bcftools view -H "$INDEL_PASS" | wc -l)
echo "[Strelka2] PASS indels: $INDEL_PASS_COUNT"

# Merge SNVs and indels into single VCF
echo "[Strelka2] Merging SNVs and indels..."
bcftools concat -a -Oz -o "$MERGED_PASS" "$SNV_PASS" "$INDEL_PASS"
bcftools sort -Oz -o "${MERGED_PASS}.tmp" "$MERGED_PASS"
mv "${MERGED_PASS}.tmp" "$MERGED_PASS"
bcftools index -t "$MERGED_PASS"

MERGED_PASS_COUNT=$(bcftools view -H "$MERGED_PASS" | wc -l)
echo "[Strelka2] Total PASS variants (merged): $MERGED_PASS_COUNT"

# ---- optionally keep unfiltered versions ----
if [[ "$KEEP_ALL" == "true" ]]; then
  echo "[Strelka2] Copying unfiltered VCFs..."
  
  SNV_ALL="${SAMPLE_OUT}/${SAMPLE_ID}.strelka2.snvs.all.vcf.gz"
  INDEL_ALL="${SAMPLE_OUT}/${SAMPLE_ID}.strelka2.indels.all.vcf.gz"
  MERGED_ALL="${SAMPLE_OUT}/${SAMPLE_ID}.strelka2.all.vcf.gz"
  
  cp "$SNV_VCF_RAW" "$SNV_ALL"
  cp "${SNV_VCF_RAW}.tbi" "${SNV_ALL}.tbi" 2>/dev/null || bcftools index -t "$SNV_ALL"
  
  cp "$INDEL_VCF_RAW" "$INDEL_ALL"
  cp "${INDEL_VCF_RAW}.tbi" "${INDEL_ALL}.tbi" 2>/dev/null || bcftools index -t "$INDEL_ALL"
  
  # Merge all variants
  bcftools concat -a -Oz -o "$MERGED_ALL" "$SNV_ALL" "$INDEL_ALL"
  bcftools sort -Oz -o "${MERGED_ALL}.tmp" "$MERGED_ALL"
  mv "${MERGED_ALL}.tmp" "$MERGED_ALL"
  bcftools index -t "$MERGED_ALL"
  
  SNV_ALL_COUNT=$(bcftools view -H "$SNV_ALL" | wc -l)
  INDEL_ALL_COUNT=$(bcftools view -H "$INDEL_ALL" | wc -l)
  MERGED_ALL_COUNT=$(bcftools view -H "$MERGED_ALL" | wc -l)
  
  echo "[Strelka2] All SNVs (unfiltered): $SNV_ALL_COUNT"
  echo "[Strelka2] All indels (unfiltered): $INDEL_ALL_COUNT"
  echo "[Strelka2] All variants (merged, unfiltered): $MERGED_ALL_COUNT"
fi

# ---- summary ----
echo ""
echo "[Strelka2] ========== Summary =========="
echo "[Strelka2] Sample: $SAMPLE_ID"
echo "[Strelka2] PASS SNVs: $SNV_PASS_COUNT"
echo "[Strelka2] PASS indels: $INDEL_PASS_COUNT"
echo "[Strelka2] Total PASS: $MERGED_PASS_COUNT"
echo "[Strelka2]"
echo "[Strelka2] Output files:"
echo "[Strelka2]   $SNV_PASS"
echo "[Strelka2]   $INDEL_PASS"
echo "[Strelka2]   $MERGED_PASS"
if [[ "$KEEP_ALL" == "true" ]]; then
  echo "[Strelka2]   $SNV_ALL"
  echo "[Strelka2]   $INDEL_ALL"
  echo "[Strelka2]   $MERGED_ALL"
fi
echo "[Strelka2] =============================="
echo "[Strelka2] done"
