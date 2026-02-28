#!/usr/bin/env bash
#set -euo pipefail

usage() {
  cat <<EOF
Usage:
  mutect2_gather.sh \
    --sample-id S1 --normal-id N1 \
    --tumor-bam T.bam --normal-bam N.bam \
    --ref ref.fa \
    --common-snps commonsnps.vcf.gz \
    --scatter-vcfs "S1.chr1.unfiltered.vcf.gz S1.chr2.unfiltered.vcf.gz ..." \
    --scatter-stats "S1.chr1.unfiltered.vcf.gz.stats S1.chr2.unfiltered.vcf.gz.stats ..." \
    --scatter-f1r2 "S1.chr1.f1r2.tar.gz S1.chr2.f1r2.tar.gz ..." \
    --out-prefix S1.mutect2 \
    [--pass-only true|false] \
    [--funcotator-datasources /path/to/ds] \
    [--funcotator-ref-version hg38|hg19] \
    [--keep-intermediates true|false]

Outputs:
  <out-prefix>.unfiltered.vcf.gz
  <out-prefix>.filtered.vcf.gz
  <out-prefix>.pass.vcf.gz
  <out-prefix>.contamination.table
  <out-prefix>.read-orientation-model.tar.gz
  <out-prefix>.metrics.json
  <out-prefix>.log
  optionally:
    <out-prefix>.pass.funcotator.vcf.gz
EOF
}

PASS_ONLY="true"
KEEP_INTERMEDIATES="false"
FUNCOTATOR_DS=""
FUNCOTATOR_REFVER=""

SAMPLE_ID=""
NORMAL_ID=""
TBAM=""
NBAM=""
REF=""
COMMONSNPS=""
SCATTER_VCFS=""
SCATTER_STATS=""
SCATTER_F1R2=""
OUT_PREFIX=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample-id) SAMPLE_ID="$2"; shift 2;;
    --normal-id) NORMAL_ID="$2"; shift 2;;
    --tumor-bam) TBAM="$2"; shift 2;;
    --normal-bam) NBAM="$2"; shift 2;;
    --ref) REF="$2"; shift 2;;
    --common-snps) COMMONSNPS="$2"; shift 2;;
    --scatter-vcfs) SCATTER_VCFS="$2"; shift 2;;
    --scatter-stats) SCATTER_STATS="$2"; shift 2;;
    --scatter-f1r2) SCATTER_F1R2="$2"; shift 2;;
    --out-prefix) OUT_PREFIX="$2"; shift 2;;
    --pass-only) PASS_ONLY="$2"; shift 2;;
    --keep-intermediates) KEEP_INTERMEDIATES="$2"; shift 2;;
    --funcotator-datasources) FUNCOTATOR_DS="$2"; shift 2;;
    --funcotator-ref-version) FUNCOTATOR_REFVER="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 2;;
  esac
done

LOG="${OUT_PREFIX}.log"
echo "[Mutect2 gather] sample=$SAMPLE_ID" | tee "$LOG"

# validate basics
for f in "$TBAM" "$NBAM" "$REF" "${REF}.fai" "$COMMONSNPS"; do
  [[ -f "$f" ]] || { echo "Missing required file: $f" >&2; exit 10; }
done

# If commonsnps is bgzipped VCF, require tabix index too
if [[ "$COMMONSNPS" == *.vcf.gz ]]; then
  [[ -f "${COMMONSNPS}.tbi" ]] || { echo "Missing required index: ${COMMONSNPS}.tbi" >&2; exit 10; }
fi

echo $SCATTER_VCFS

[[ -n "$SCATTER_VCFS" && -n "$SCATTER_STATS" && -n "$SCATTER_F1R2" && -n "$OUT_PREFIX" ]] || {
  echo "Missing required arguments" >&2; usage; exit 11;
}

UNFILTERED="${OUT_PREFIX}.unfiltered.vcf.gz"
FILTERED="${OUT_PREFIX}.filtered.vcf.gz"
PASSVCF="${OUT_PREFIX}.pass.vcf.gz"
ORIENT="${OUT_PREFIX}.read-orientation-model.tar.gz"
CONT_T="${OUT_PREFIX}.contamination.table"

# 1) concat VCFs (keep order deterministic; bcftools concat requires consistent contigs)
#    recommend: ensure your scatter outputs were generated from non-overlapping intervals.
echo "[Mutect2 gather] concat vcfs -> $UNFILTERED" | tee -a "$LOG"
# bcftools concat writes bgzip output with -Oz, then index
bcftools concat -a -Oz -o "$UNFILTERED" $SCATTER_VCFS >>"$LOG" 2>&1
tabix -p vcf "$UNFILTERED" >>"$LOG" 2>&1

# 2) learn orientation model
echo "[Mutect2 gather] LearnReadOrientationModel -> $ORIENT" | tee -a "$LOG"
gatk LearnReadOrientationModel $SCATTER_F1R2 -O "$ORIENT" >>"$LOG" 2>&1

# 3) contamination
echo "[Mutect2 gather] GetPileupSummaries / CalculateContamination" | tee -a "$LOG"
gatk GetPileupSummaries \
  -I "$TBAM" -V "$COMMONSNPS" -L "$COMMONSNPS" \
  -O "${OUT_PREFIX}.tumor.pileups.table" >>"$LOG" 2>&1

gatk GetPileupSummaries \
  -I "$NBAM" -V "$COMMONSNPS" -L "$COMMONSNPS" \
  -O "${OUT_PREFIX}.normal.pileups.table" >>"$LOG" 2>&1

gatk CalculateContamination \
  -I "${OUT_PREFIX}.tumor.pileups.table" \
  -matched "${OUT_PREFIX}.normal.pileups.table" \
  -O "$CONT_T" >>"$LOG" 2>&1

# 4) merge stats
echo "[Mutect2 gather] MergeMutectStats" | tee -a "$LOG"
gatk MergeMutectStats $SCATTER_STATS -O "${UNFILTERED}.stats" >>"$LOG" 2>&1

# 5) filter mutect calls
echo "[Mutect2 gather] FilterMutectCalls -> $FILTERED" | tee -a "$LOG"
gatk FilterMutectCalls \
  -R "$REF" \
  -V "$UNFILTERED" \
  --contamination-table "$CONT_T" \
  --stats "${UNFILTERED}.stats" \
  --ob-priors "$ORIENT" \
  -O "$FILTERED" >>"$LOG" 2>&1

# 6) select PASS (optional)
if [[ "$PASS_ONLY" == "true" ]]; then
  echo "[Mutect2 gather] SelectVariants PASS-only -> $PASSVCF" | tee -a "$LOG"
  gatk SelectVariants \
    -R "$REF" \
    -V "$FILTERED" \
    --exclude-filtered true \
    -O "$PASSVCF" >>"$LOG" 2>&1
else
  cp "$FILTERED" "$PASSVCF"
  [[ -f "${FILTERED}.tbi" ]] && cp "${FILTERED}.tbi" "${PASSVCF}.tbi" || true
fi

# ensure indexed
tabix -p vcf "$PASSVCF" >>"$LOG" 2>&1 || true

# 7) optional funcotator
FUNC_OUT=""
if [[ -n "$FUNCOTATOR_DS" ]]; then
  [[ -n "$FUNCOTATOR_REFVER" ]] || { echo "funcotator-ref-version required if datasources provided" | tee -a "$LOG"; exit 20; }
  
  # Write uncompressed VCF first, then bgzip+tabix (Funcotator does not reliably create a valid .vcf.gz directly)
  FUNC_VCF="${OUT_PREFIX}.pass.funcotator.vcf"
  FUNC_OUT="${FUNC_VCF}.gz"

  echo "[Mutect2 gather] Funcotator -> $FUNC_OUT" | tee -a "$LOG"
  gatk Funcotator \
    -R "$REF" \
    -V "$PASSVCF" \
    -O "$FUNC_VCF" \
    --output-file-format VCF \
    --data-sources-path "$FUNCOTATOR_DS" \
    --ref-version "$FUNCOTATOR_REFVER" >>"$LOG" 2>&1

  bgzip -f "$FUNC_VCF" >>"$LOG" 2>&1
  tabix -p vcf "$FUNC_OUT" >>"$LOG" 2>&1
fi

# 8) minimal metrics
echo "[Mutect2 gather] write metrics.json" | tee -a "$LOG"
unfiltered_n=$(bcftools view -H "$UNFILTERED" | wc -l | tr -d ' ')
pass_n=$(bcftools view -H "$PASSVCF" | wc -l | tr -d ' ')

cat > "${OUT_PREFIX}.metrics.json" <<EOF
{
  "sample_id": "${SAMPLE_ID}",
  "caller": "mutect2",
  "unfiltered_variants": ${unfiltered_n},
  "pass_variants": ${pass_n},
  "pass_only": ${PASS_ONLY},
  "funcotator_enabled": $( [[ -n "$FUNCOTATOR_DS" ]] && echo "true" || echo "false" )
}
EOF

if [[ "$KEEP_INTERMEDIATES" != "true" ]]; then
  rm -f "${OUT_PREFIX}.tumor.pileups.table" "${OUT_PREFIX}.normal.pileups.table" || true
fi

echo "[Mutect2 gather] done" | tee -a "$LOG"
