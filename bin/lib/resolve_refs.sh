#!/usr/bin/env bash
#set -euo pipefail

# resolve_refs.sh
# Usage:
#   source resolve_refs.sh --ref hg38 --refdir ./refs [--seq WGS|EXOME]
#
# Exports:
#   REF_GENOME, GNOMAD, PON, DBSNP, INDELS, INDEL_1000G, INDEL_KNOWN, COMMONSNPS
#   contigs (bash array), chr (prefix string)

usage_resolve_refs() {
  cat <<'EOF'
resolve_refs.sh: map ref build to reference files under a user-provided directory.

Required:
  --ref    hg38|hg19|mm10|mm39
  --refdir PATH   (directory containing reference assets)

Optional:
  --seq    WGS|EXOME   (only used where PoN differs; default WGS)

Expected refdir layouts (recommended):
  refs/
    hg38/
      Homo_sapiens_assembly38.fasta
      Homo_sapiens_assembly38.fasta.fai
      Homo_sapiens_assembly38.dict
      somatic-hg38-af-only-gnomad.hg38.vcf.gz
      somatic-hg38-af-only-gnomad.hg38.vcf.gz.tbi
      somatic-hg38-1000g_pon.hg38.vcf.gz
      somatic-hg38-1000g_pon.hg38.vcf.gz.tbi
      somatic-hg38-small_exac_common_3.hg38.vcf.gz
      somatic-hg38-small_exac_common_3.hg38.vcf.gz.tbi
      (optional for BQSR etc)
      Homo_sapiens_assembly38.dbsnp138.vcf.gz
      Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
      Homo_sapiens_assembly38.known_indels.vcf.gz

    hg19/
      somatic-b37-Homo_sapiens_assembly19.fasta
      somatic-b37-Homo_sapiens_assembly19.fasta.fai
      somatic-b37-Homo_sapiens_assembly19.dict
      somatic-b37-af-only-gnomad.raw.sites.vcf.gz (or .vcf)
      somatic-b37-af-only-gnomad.raw.sites.vcf.gz.tbi (if gz)
      somatic-b37-small_exac_common_3.vcf.gz (or .vcf)
      somatic-b37-small_exac_common_3.vcf.gz.tbi (if gz)
      somatic-b37-Mutect2-WGS-panel-b37.vcf.gz (or .vcf)
      somatic-b37-Mutect2-exome-panel.vcf.gz (or .vcf)

    mm10/
      mm10.fa
      mm10.fa.fai
      mm10.dict
      mgp.v3.snps.rsIDdbSNPv137.mm10.vcf.gz (or .vcf)
      mgp.v3.indels.rsIDdbSNPv137.mm10.vcf.gz (or .vcf)

    mm39/
      mm39.fa
      mm39.fa.fai
      mm39.dict
      mgp_REL2021_snps.rsID.renamed.vcf.gz
      mgp_REL2021_indels.rsID.renamed.vcf.gz
      mm39_PON.vcf.gz

EOF
}

# -------- parse args --------
REF=""
REFDIR=""
SEQ="WGS"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --ref) REF="$2"; shift 2;;
    --refdir) REFDIR="$2"; shift 2;;
    --seq) SEQ="$2"; shift 2;;
    -h|--help) usage_resolve_refs; return 0;;
    *) echo "resolve_refs.sh: Unknown arg: $1" >&2; usage_resolve_refs; return 2;;
  esac
done

[[ -n "$REF" && -n "$REFDIR" ]] || { echo "resolve_refs.sh: --ref and --refdir are required" >&2; usage_resolve_refs; return 2; }
[[ -d "$REFDIR" ]] || { echo "resolve_refs.sh: refdir not found: $REFDIR" >&2; return 2; }

# -------- helpers --------
require_file() {
  local f="$1"
  [[ -f "$f" ]] || { echo "Missing required reference file: $f" >&2; return 10; }
}

# Allow either .vcf or .vcf.gz; prefer .vcf.gz if present.
pick_vcf() {
  local base="$1"
  if [[ -f "${base}.vcf.gz" ]]; then
    echo "${base}.vcf.gz"
  elif [[ -f "${base}.vcf" ]]; then
    echo "${base}.vcf"
  else
    # caller will error later with require_file on returned path
    echo "${base}.vcf.gz"
  fi
}

# If gz, require tbi too (Mutect2 expects indexed resources typically)
maybe_require_index() {
  local f="$1"
  if [[ "$f" == *.vcf.gz ]]; then
    require_file "${f}.tbi"
  fi
}

# -------- resolve --------
REFROOT="${REFDIR%/}/${REF}"

unset GNOMAD REF_GENOME PON DBSNP INDEL_1000G INDEL_KNOWN INDELS COMMONSNPS chrlist contigs chr

case "$REF" in
  hg38)
    REF_GENOME="${REFROOT}/Homo_sapiens_assembly38.fasta"
    GNOMAD="$(pick_vcf "${REFROOT}/somatic-hg38-af-only-gnomad.hg38")"
    PON="$(pick_vcf "${REFROOT}/somatic-hg38-1000g_pon.hg38")"
    COMMONSNPS="$(pick_vcf "${REFROOT}/somatic-hg38-small_exac_common_3.hg38")"

    # optional extras (don’t hard-require unless your pipeline needs them)
    DBSNP="$(pick_vcf "${REFROOT}/Homo_sapiens_assembly38.dbsnp138")"
    INDEL_1000G="$(pick_vcf "${REFROOT}/Mills_and_1000G_gold_standard.indels.hg38")"
    INDEL_KNOWN="$(pick_vcf "${REFROOT}/Homo_sapiens_assembly38.known_indels")"

    chr="chr"
    contigs=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y")
    ;;

  hg19)
    REF_GENOME="${REFROOT}/somatic-b37-Homo_sapiens_assembly19.fasta"
    GNOMAD="$(pick_vcf "${REFROOT}/somatic-b37-af-only-gnomad.raw.sites")"
    COMMONSNPS="$(pick_vcf "${REFROOT}/somatic-b37-small_exac_common_3")"

    # PoN depends on seq
    if [[ "$SEQ" == "WGS" ]]; then
      PON="$(pick_vcf "${REFROOT}/somatic-b37-Mutect2-WGS-panel-b37")"
    else
      PON="$(pick_vcf "${REFROOT}/somatic-b37-Mutect2-exome-panel")"
    fi

    DBSNP="$(pick_vcf "${REFROOT}/Homo_sapiens_assembly19.dbsnp138")"
    INDEL_1000G="$(pick_vcf "${REFROOT}/Mills_and_1000G_gold_standard.indels.b37")"
    INDEL_KNOWN="$(pick_vcf "${REFROOT}/Homo_sapiens_assembly19.known_indels")"

    chr=""
    contigs=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y")
    ;;

  mm10)
    REF_GENOME="${REFROOT}/mm10.fa"
    COMMONSNPS="$(pick_vcf "${REFROOT}/mgp.v3.snps.rsIDdbSNPv137.mm10")"
    INDELS="$(pick_vcf "${REFROOT}/mgp.v3.indels.rsIDdbSNPv137.mm10")"
    chr=""
    contigs=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "X" "Y")
    ;;

  mm39)
    REF_GENOME="${REFROOT}/mm39.fa"
    COMMONSNPS="$(pick_vcf "${REFROOT}/mgp_REL2021_snps.rsID.renamed")"
    INDELS="$(pick_vcf "${REFROOT}/mgp_REL2021_indels.rsID.renamed")"
    PON="$(pick_vcf "${REFROOT}/mm39_PON")"
    chr="chr"
    contigs=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "X" "Y")
    ;;

  *)
    echo "resolve_refs.sh: ref not found: $REF" >&2
    usage_resolve_refs
    return 2
    ;;
esac

# ---- require the essentials for Mutect2 consensus workflow ----
require_file "$REF_GENOME"
require_file "${REF_GENOME}.fai"
# dict name can vary; require either <ref>.dict or explicit common name if present
if [[ -f "${REF_GENOME%.fasta}.dict" ]]; then
  : # ok
elif [[ -f "${REF_GENOME}.dict" ]]; then
  : # ok
elif [[ -f "${REFROOT}/$(basename "${REF_GENOME%.*}").dict" ]]; then
  : # ok
else
  echo "Missing reference dict for: $REF_GENOME (expected .dict in same folder)" >&2
  return 10
fi

# Mutect2 essentials for human/mouse:
# - germline resource + PoN for human; for mouse mm39 you also set PoN.
# - common snps required for contamination calculation in your workflow.
if [[ "$REF" == "hg38" || "$REF" == "hg19" ]]; then
  require_file "$GNOMAD"; maybe_require_index "$GNOMAD"
  require_file "$PON"; maybe_require_index "$PON"
fi

require_file "$COMMONSNPS"; maybe_require_index "$COMMONSNPS"
if [[ -n "${INDELS:-}" && -f "$INDELS" ]]; then
  maybe_require_index "$INDELS"
fi

export REF_GENOME GNOMAD PON DBSNP INDEL_1000G INDEL_KNOWN INDELS COMMONSNPS chr contigs
