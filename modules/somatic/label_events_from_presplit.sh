#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
label_events_from_presplit.sh

Compare pre-split VCF (std.norm.vcf.gz) vs final split VCF (std.norm.split.mnvsplit.vcf.gz),
detect both MULTIALLELIC and MNV-derived atomic records, and add EVENT tags to final VCF.

Required:
  --presplit   PATH   pre-split VCF (bgzip vcf.gz preferred)
  --final      PATH   final VCF after split+mnvsplit (bgzip vcf.gz preferred)
  --out        PATH   output annotated VCF.gz

Optional:
  --threads INT       default: 4
  --workdir PATH      default: <out>.workdir (kept, NOT deleted)

Outputs:
  <out> and <out>.tbi
  Workdir contains intermediate TSVs for debugging.

Event INFO added (if matched):
  EVENT_ID, EVENT_TYPE, EVENT_PART, EVENT_SIZE, EVENT_STATUS (FULL/PARTIAL)

Matching key:
  CHROM, POS, REF, ALT  (ALT assumed biallelic in final; still robust if comma exists)

EOF
}

PRESPLIT=""
FINAL=""
OUT=""
THREADS=4
WORKDIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --presplit) PRESPLIT="$2"; shift 2;;
    --final) FINAL="$2"; shift 2;;
    --out) OUT="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --workdir) WORKDIR="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1" >&2; usage; exit 2;;
  esac
done

[[ -n "$PRESPLIT" && -n "$FINAL" && -n "$OUT" ]] || { echo "ERROR: missing required args" >&2; usage; exit 2; }

command -v bcftools >/dev/null 2>&1 || { echo "ERROR: bcftools not found" >&2; exit 127; }
command -v bgzip   >/dev/null 2>&1 || { echo "ERROR: bgzip not found" >&2; exit 127; }
command -v tabix   >/dev/null 2>&1 || { echo "ERROR: tabix not found" >&2; exit 127; }

[[ -f "$PRESPLIT" ]] || { echo "ERROR: presplit not found: $PRESPLIT" >&2; exit 2; }
[[ -f "$FINAL" ]] || { echo "ERROR: final not found: $FINAL" >&2; exit 2; }

if [[ -z "$WORKDIR" ]]; then
  WORKDIR="${OUT}.workdir"
fi
mkdir -p "$WORKDIR"

log="$WORKDIR/label_events.log"
ts() { date +"%F %T"; }
echo "[$(ts)] label_events_from_presplit starting" | tee -a "$log" >&2
echo "[$(ts)] presplit=$PRESPLIT" | tee -a "$log" >&2
echo "[$(ts)] final=$FINAL" | tee -a "$log" >&2
echo "[$(ts)] out=$OUT" | tee -a "$log" >&2
echo "[$(ts)] workdir=$WORKDIR" | tee -a "$log" >&2

# Ensure FINAL is bgz + indexed (needed for annotation output sanity)
if [[ "$FINAL" != *.vcf.gz ]]; then
  echo "ERROR: --final must be bgzipped .vcf.gz for robust downstream; got: $FINAL" | tee -a "$log" >&2
  exit 2
fi
if [[ ! -f "${FINAL}.tbi" && ! -f "${FINAL}.csi" ]]; then
  echo "[$(ts)] Indexing FINAL VCF" | tee -a "$log" >&2
  tabix -f -p vcf "$FINAL" >>"$log" 2>&1
fi

# -----------------------
# 1) Extract FINAL keys
# -----------------------
echo "[$(ts)] Extracting FINAL variant keys" | tee -a "$log" >&2
final_keys="$WORKDIR/final.keys.tsv"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$FINAL" \
| awk 'BEGIN{FS=OFS="\t"} NF==4 && $4!="."{
    # split ALT just in case (should already be biallelic)
    n=split($4,alts,",");
    for(i=1;i<=n;i++) print $1,$2,$3,alts[i]
  }' > "$final_keys"

final_count=$(wc -l < "$final_keys")
echo "[$(ts)] FINAL has $final_count variant keys" | tee -a "$log" >&2

if [[ ! -s "$final_keys" ]]; then
  echo "[$(ts)] WARNING: final has no variants; writing header-only output" | tee -a "$log" >&2
  bcftools view -h "$FINAL" -Ov | bgzip -c > "$OUT"
  tabix -f -p vcf "$OUT"
  exit 0
fi

# ---------------------------------------
# 2) Expand events from PRESPLIT to atomic
#    Output columns:
#    CHROM POS REF ALT EVENT_ID EVENT_TYPE EVENT_PART EVENT_SIZE
# ---------------------------------------
echo "[$(ts)] Expanding PRESPLIT events (MULTIALLELIC and MNV)" | tee -a "$log" >&2
expanded="$WORKDIR/presplit.expanded.tsv"

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$PRESPLIT" \
| awk '
BEGIN{FS=OFS="\t"}
function is_acgtn(s,   i,c){
  for(i=1;i<=length(s);i++){
    c=substr(s,i,1)
    if(c!="A" && c!="C" && c!="G" && c!="T" && c!="N") return 0
  }
  return 1
}
NF==4 && $4!="."{
  chrom=$1; pos=$2; ref=$3; alt=$4

  # MULTIALLELIC: ALT has commas. Expand each allele as its own atomic key.
  if(index(alt,",")>0){
    n=split(alt,alts,",")
    event_id=chrom ":" pos ":" ref ">" alt
    for(i=1;i<=n;i++){
      print chrom,pos,ref,alts[i],event_id,"MULTIALLELIC",i,n
    }
    next
  }

  # MNV candidate: same length, length>1, pure bases
  if(length(ref)>1 && length(ref)==length(alt) && is_acgtn(ref) && is_acgtn(alt)){
    L=length(ref)
    event_id=chrom ":" pos ":" ref ">" alt
    for(i=1;i<=L;i++){
      r=substr(ref,i,1)
      a=substr(alt,i,1)
      # skip "no-change" positions if any (rare, but makes mapping cleaner)
      if(r==a) continue
      print chrom, (pos+i-1), r, a, event_id, "MNV", i, L
    }
    next
  }

  # Otherwise: no event expansion emitted
}
' > "$expanded"

expanded_count=$(wc -l < "$expanded")
echo "[$(ts)] Expanded to $expanded_count atomic event parts" | tee -a "$log" >&2

# expanded can be empty if there are no multiallelics/MNV in presplit — that's fine.

# -------------------------------------------------------
# 3) Intersect expanded events with FINAL keys,
#    compute EVENT_STATUS (FULL vs PARTIAL per EVENT_ID),
#    and build an annotation TSV for bcftools annotate.
#
#    Annotation file columns (CORRECTED):
#      CHROM POS REF ALT EVENT_ID EVENT_TYPE EVENT_PART EVENT_SIZE EVENT_STATUS
# -------------------------------------------------------
echo "[$(ts)] Computing EVENT_STATUS and building annotation file" | tee -a "$log" >&2
anno="$WORKDIR/event_anno.tsv"
anno_gz="${anno}.gz"

awk '
BEGIN{FS=OFS="\t"}
FNR==NR{
  # final keys
  k=$1 OFS $2 OFS $3 OFS $4
  final[k]=1
  next
}
NF==8{
  chrom=$1; pos=$2; ref=$3; alt=$4
  event_id=$5; event_type=$6; part=$7; size=$8

  k=chrom OFS pos OFS ref OFS alt
  # keep only expansions that exist in final
  if(!(k in final)) next

  # store record lines to output later, but first compute counts per event
  line_count++
  rec_k[line_count]=k
  rec_event[line_count]=event_id
  rec_type[line_count]=event_type
  rec_part[line_count]=part
  rec_size[line_count]=size

  # size should be consistent per event_id; keep max just in case
  if(!(event_id in esize) || size>esize[event_id]) esize[event_id]=size
  epresent[event_id]++
}
END{
  # Output annotation lines for each matched record
  # CORRECTED: Each INFO field as separate column
  for(i=1;i<=line_count;i++){
    k=rec_k[i]
    split(k,a,OFS)
    event_id=rec_event[i]
    size=rec_size[i]
    present=epresent[event_id]
    status=(present==esize[event_id] ? "FULL" : "PARTIAL")

    # Output: CHROM POS REF ALT EVENT_ID EVENT_TYPE EVENT_PART EVENT_SIZE EVENT_STATUS
    print a[1], a[2], a[3], a[4], event_id, rec_type[i], rec_part[i], esize[event_id], status
  }
}
' "$final_keys" "$expanded" > "$anno"

anno_count=$(wc -l < "$anno")
echo "[$(ts)] Matched $anno_count variants with EVENT annotations" | tee -a "$log" >&2

# ---------------------------------------
# 4) Annotate FINAL with EVENT tags
# ---------------------------------------
header_add="$WORKDIR/event_header.hdr"
cat > "$header_add" <<'HDR'
##INFO=<ID=EVENT_ID,Number=1,Type=String,Description="Event identifier for grouped variants (MNV or MULTIALLELIC)">
##INFO=<ID=EVENT_TYPE,Number=1,Type=String,Description="Event type: MNV or MULTIALLELIC">
##INFO=<ID=EVENT_PART,Number=1,Type=Integer,Description="1-based part index within the event">
##INFO=<ID=EVENT_SIZE,Number=1,Type=Integer,Description="Total number of parts in the event">
##INFO=<ID=EVENT_STATUS,Number=1,Type=String,Description="FULL if all parts are present in final VCF, else PARTIAL">
HDR

# Sort annotation file
echo "[$(ts)] Sorting annotation file" | tee -a "$log" >&2
sort -k1,1V -k2,2n "$anno" > "${anno}.sorted"
mv "${anno}.sorted" "$anno"

# Compress and index
bgzip -f -c "$anno" > "$anno_gz"
tabix -f -s 1 -b 2 -e 2 "$anno_gz" >>"$log" 2>&1

if [[ -s "$anno" ]]; then
  echo "[$(ts)] Annotating FINAL with EVENT tags (matched lines: $anno_count)" | tee -a "$log" >&2
  # CORRECTED: Column specification matches the actual TSV format
  bcftools annotate \
    -a "$anno_gz" \
    -c CHROM,POS,REF,ALT,INFO/EVENT_ID,INFO/EVENT_TYPE,INFO/EVENT_PART,INFO/EVENT_SIZE,INFO/EVENT_STATUS \
    -h "$header_add" \
    --threads "$THREADS" \
    -Oz -o "$OUT" \
    "$FINAL" >>"$log" 2>&1
else
  echo "[$(ts)] No event matches found; writing FINAL through with header additions only." | tee -a "$log" >&2
  # Add header lines even if no records match
  bcftools annotate \
    -h "$header_add" \
    --threads "$THREADS" \
    -Oz -o "$OUT" \
    "$FINAL" >>"$log" 2>&1
fi

echo "[$(ts)] Indexing output VCF" | tee -a "$log" >&2
tabix -f -p vcf "$OUT" >>"$log" 2>&1

echo "[$(ts)] Wrote: $OUT" | tee -a "$log" >&2
echo "[$(ts)] Index: ${OUT}.tbi" | tee -a "$log" >&2
echo "[$(ts)] Summary:" | tee -a "$log" >&2
echo "[$(ts)]   FINAL variants: $final_count" | tee -a "$log" >&2
echo "[$(ts)]   Expanded event parts: $expanded_count" | tee -a "$log" >&2
echo "[$(ts)]   Matched annotations: $anno_count" | tee -a "$log" >&2
echo "[$(ts)] Done." | tee -a "$log" >&2