#!/usr/bin/env bash
set -euo pipefail

# split_mnv.sh
# Split MNVs (same-length REF/ALT where length>1) into per-base SNV records.
#
# Notes:
# - This is a "consensus-friendly" split: it preserves FORMAT/sample columns as-is
#   (which may not be biologically exact after splitting). This is usually fine if
#   you use the split VCF for K-of-N voting, then later emit a merged/clean VCF.
# - Indels (REF/ALT different length) are not split.
# - Multiallelic records (ALT contains comma) are NOT expanded here; use bcftools norm -m -any before this script.

usage() {
  cat <<'EOF'
split_mnv.sh

Split MNVs (e.g., CA->TG) into per-base SNVs while preserving header and producing bgzip VCF.

Required:
  --in-vcf   PATH   input VCF/VCF.GZ
  --out      PATH   output .vcf.gz

Optional:
  --workdir  PATH   keep intermediate files here (default: <out>.work)
  --keep-workdir true|false  default: true
  --add-orig-tags true|false default: true
  --threads  INT    default: 4

Output:
  <out> (bgzip VCF) + <out>.tbi

Tags added (if --add-orig-tags=true):
  INFO/ORIG_EVENT = CHROM:POS:REF>ALT
  INFO/ORIG_POS   = original POS
  INFO/ORIG_REF   = original REF
  INFO/ORIG_ALT   = original ALT
  INFO/ORIG_LEN   = original REF length
  INFO/SPLIT_MNV  = 1

EOF
}

IN_VCF=""
OUT=""
WORKDIR=""
KEEP_WORKDIR="true"
ADD_ORIG="true"
THREADS=4

while [[ $# -gt 0 ]]; do
  case "$1" in
    --in-vcf) IN_VCF="$2"; shift 2;;
    --out) OUT="$2"; shift 2;;
    --workdir) WORKDIR="$2"; shift 2;;
    --keep-workdir) KEEP_WORKDIR="$2"; shift 2;;
    --add-orig-tags) ADD_ORIG="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "split_mnv.sh: Unknown arg: $1" >&2; usage; exit 2;;
  esac
done

[[ -n "$IN_VCF" && -n "$OUT" ]] || { echo "ERROR: --in-vcf and --out are required" >&2; usage; exit 2; }
[[ -f "$IN_VCF" ]] || { echo "ERROR: input not found: $IN_VCF" >&2; exit 2; }

if [[ -z "$WORKDIR" ]]; then
  WORKDIR="${OUT}.work"
fi
mkdir -p "$WORKDIR"

cleanup() {
  if [[ "$KEEP_WORKDIR" != "true" ]]; then
    rm -rf "$WORKDIR"
  fi
}
trap cleanup EXIT

LOG="${WORKDIR}/split_mnv.log"
echo "[split_mnv] IN:  $IN_VCF" | tee "$LOG"
echo "[split_mnv] OUT: $OUT" | tee -a "$LOG"
echo "[split_mnv] WORKDIR: $WORKDIR" | tee -a "$LOG"

# Ensure bgzip input + index (bcftools view can read plain VCF, but we prefer .vcf.gz)
INPUT_GZ="$IN_VCF"
if [[ "$IN_VCF" != *.vcf.gz ]]; then
  INPUT_GZ="${WORKDIR}/input.vcf.gz"
  bcftools view --threads "$THREADS" -Oz -o "$INPUT_GZ" "$IN_VCF" >>"$LOG" 2>&1
  bcftools index -f -t "$INPUT_GZ" >>"$LOG" 2>&1
else
  if [[ ! -f "${IN_VCF}.tbi" && ! -f "${IN_VCF}.csi" ]]; then
    bcftools index -f -t "$IN_VCF" >>"$LOG" 2>&1
  fi
fi

HEADER="${WORKDIR}/header.txt"
BODY="${WORKDIR}/body.tsv"
OUT_TMP="${WORKDIR}/split.tmp.vcf"

# Extract header + body
bcftools view -h "$INPUT_GZ" > "$HEADER"
bcftools view -H "$INPUT_GZ" > "$BODY"

# Add INFO header lines if requested (insert before #CHROM)
if [[ "$ADD_ORIG" == "true" ]]; then
  awk '
    BEGIN{added=0}
    /^#CHROM/ && added==0{
      print "##INFO=<ID=ORIG_EVENT,Number=1,Type=String,Description=\"Original unsplit event CHROM:POS:REF>ALT\">"
      print "##INFO=<ID=ORIG_POS,Number=1,Type=Integer,Description=\"Original POS before splitting\">"
      print "##INFO=<ID=ORIG_REF,Number=1,Type=String,Description=\"Original REF before splitting\">"
      print "##INFO=<ID=ORIG_ALT,Number=1,Type=String,Description=\"Original ALT before splitting\">"
      print "##INFO=<ID=ORIG_LEN,Number=1,Type=Integer,Description=\"Length of original REF for split MNV\">"
      print "##INFO=<ID=SPLIT_MNV,Number=0,Type=Flag,Description=\"Record produced by splitting an MNV\">"
      added=1
    }
    {print}
  ' "$HEADER" > "${WORKDIR}/header.withinfo.txt"
  mv "${WORKDIR}/header.withinfo.txt" "$HEADER"
fi

# Split MNVs in body
# VCF columns:
# 1 CHROM 2 POS 3 ID 4 REF 5 ALT 6 QUAL 7 FILTER 8 INFO 9 FORMAT 10+ samples...
#
# Rules:
# - If ALT contains comma -> leave unchanged (assume user ran bcftools norm -m -any beforehand)
# - If length(REF)==length(ALT)>1 -> split into L lines at POS+i with REF[i], ALT[i]
# - Else -> pass-through

awk -v add_orig="$ADD_ORIG" '
BEGIN{FS=OFS="\t"}
{
  chrom=$1; pos=$2; id=$3; ref=$4; alt=$5;
  qual=$6; filt=$7; info=$8;

  # leave multi-allelic as-is
  if (index(alt,",")>0) { print; next }

  lref=length(ref); lalt=length(alt);

  # split only same-length MNVs (not indels), length>1
  if (lref==lalt && lref>1) {
    orig_event = chrom ":" pos ":" ref ">" alt
    for (i=1; i<=lref; i++) {
      npos = pos + (i-1)
      r = substr(ref,i,1)
      a = substr(alt,i,1)

      ninfo = info
      if (add_orig=="true") {
        # append with semicolon if INFO not empty
        if (ninfo=="." || ninfo=="") ninfo=""
        else ninfo=ninfo ";"
        ninfo = ninfo "ORIG_EVENT=" orig_event ";ORIG_POS=" pos ";ORIG_REF=" ref ";ORIG_ALT=" alt ";ORIG_LEN=" lref ";SPLIT_MNV"
        if (ninfo=="") ninfo="."
      }

      # print updated VCF line, keep remaining columns as-is
      $2=npos; $4=r; $5=a; $8=(ninfo==""?".":ninfo)
      print
    }
    next
  }

  # pass-through
  print
}
' "$BODY" > "${WORKDIR}/body.split.tsv"

# Rebuild VCF
cat "$HEADER" "${WORKDIR}/body.split.tsv" > "$OUT_TMP"

# bgzip + index
bcftools view --threads "$THREADS" -Oz -o "$OUT" "$OUT_TMP" >>"$LOG" 2>&1
bcftools index -f -t "$OUT" >>"$LOG" 2>&1

echo "[split_mnv] Wrote: $OUT" | tee -a "$LOG"
echo "[split_mnv] Index: ${OUT}.tbi" | tee -a "$LOG"
echo "[split_mnv] Log:   $LOG" | tee -a "$LOG"
