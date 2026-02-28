#!/usr/bin/env bash
set -euo pipefail

# make_somatic_metrics.sh
#
# Generate metric TSVs for plotting:
# - Variant type counts per stage/caller
# - K-of-N summary from k_pass files
# - Caller agreement matrix
# - Variant type breakdown by K threshold
# - Genomic distribution (by chr bins)

usage() {
  cat <<'EOF'
make_somatic_metrics.sh

Required:
  --outdir DIR               output directory for metrics TSVs
  --sample-id STR            sample name (used in tables)

Provide VCFs (repeatable by stage + caller):
  --vcf stage=STAGE caller=CALLER path=VCF.gz

Where STAGE is one of:
  raw | norm | split | consensus | merged | event
CALLER is one of:
  mutect2 | muse2 | strelka2 | <anything>

K-of-N inputs (from consensus_k_of_n.sh --mode compare):
  --kpass k=INT path=FILE    k_pass TSV file for a specific K threshold (repeatable)
                             Example: --kpass "k=1 path=workdir/k1.pass.tsv"
                                      --kpass "k=2 path=workdir/k2.pass.tsv"
                                      --kpass "k=3 path=workdir/k3.pass.tsv"

Optional:
  --bin-size INT             default 1000000 (1Mb)
  --keep-temp                keep intermediate files (default: delete)

Anchor VCF depth summary (adds DP stats to report):
  --anchor-vcf stage=STAGE caller=CALLER path=VCF.gz   anchor VCF to summarize depth from (once)
  --dp-sample STR                                     optional: sample name inside VCF for DP extraction
                                                     (default: first sample / bcftools default)

Outputs (in outdir):
  variant_counts_by_stage.tsv       counts by sample/stage/caller/vartype
  variant_counts_wide.tsv           wide matrix: rows=stage:caller, cols=vartype
  upset_matrix.tsv                  wide UpSet matrix (0/1 per caller) from k1.pass
  caller_agreement_matrix.tsv       pairwise caller overlap counts
  caller_contribution.tsv           variants unique/shared per caller
  k_curve.tsv                       K threshold vs variant counts
  k_curve_by_vartype.tsv            K threshold vs counts by variant type
  genomic_bins.tsv                  variant density by genomic bin
  genomic_summary.tsv               per-chromosome totals
EOF
}

OUTDIR=""
SAMPLE=""
BIN_SIZE=1000000
KEEP_TEMP="false"
ANCHOR_SPEC=""
DP_SAMPLE=""

declare -a VCF_SPECS=()
declare -a KPASS_SPECS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2;;
    --sample-id) SAMPLE="$2"; shift 2;;
    --vcf) VCF_SPECS+=("$2"); shift 2;;
    --kpass) KPASS_SPECS+=("$2"); shift 2;;
    --bin-size) BIN_SIZE="$2"; shift 2;;
    --anchor-vcf) ANCHOR_SPEC="$2"; shift 2;;
    --dp-sample) DP_SAMPLE="$2"; shift 2;;
    --keep-temp) KEEP_TEMP="true"; shift;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1" >&2; usage; exit 2;;
  esac
done

[[ -n "$OUTDIR" && -n "$SAMPLE" ]] || { 
  echo "ERROR: --outdir and --sample-id required" >&2
  usage
  exit 2
}

mkdir -p "$OUTDIR"
TMPDIR="$OUTDIR/tmp.$$"
mkdir -p "$TMPDIR"

if [[ "$KEEP_TEMP" != "true" ]]; then
  trap "rm -rf '$TMPDIR'" EXIT
fi

command -v bcftools >/dev/null 2>&1 || { echo "ERROR: bcftools not in PATH" >&2; exit 127; }
command -v awk >/dev/null 2>&1 || { echo "ERROR: awk not in PATH" >&2; exit 127; }

ts() { date +"%F %T"; }
log() { echo "[$(ts)] $*" >&2; }

log "Starting metrics generation for sample: $SAMPLE"

# ----------------------------
# Helper function: parse key=value specs
# ----------------------------
parse_spec() {
  local spec="$1"
  local key="$2"
  echo "$spec" | tr ' ' '\n' | grep "^${key}=" | cut -d= -f2-
}

# ----------------------------
# Initialize output files
# ----------------------------
counts_long="$OUTDIR/variant_counts_by_stage.tsv"
counts_wide="$OUTDIR/variant_counts_wide.tsv"
upset_matrix="$OUTDIR/upset_matrix.tsv"
caller_agreement="$OUTDIR/caller_agreement_matrix.tsv"
caller_contrib="$OUTDIR/caller_contribution.tsv"
genome_bins="$OUTDIR/genomic_bins.tsv"
genome_summary="$OUTDIR/genomic_summary.tsv"
event_by_type="$OUTDIR/event_status_by_type.tsv"
event_summary="$OUTDIR/event_status_summary.tsv"
vcf_dp_summary="$OUTDIR/vcf_dp_summary.tsv"

# K-curve outputs
kcurve="$OUTDIR/k_curve.tsv"
kcurve_vartype="$OUTDIR/k_curve_by_vartype.tsv"

##################################################
### Help function for Counting events of MNV ##############
###################################
count_event_status_from_vcf() {
  local vcf="$1"
  local stage="$2"
  local caller="$3"

  bcftools query -f '%INFO/EVENT_TYPE\t%INFO/MERGED_STATUS\n' "$vcf" 2>/dev/null \
  | awk -v sample="$SAMPLE" -v stage="$stage" -v caller="$caller" -v OFS="\t" '
    {
      t=$1; s=$2
      if(t=="." || s==".") next
      if(t!="MNV" && t!="MULTIALLELIC") next
      if(s!="FULL" && s!="PARTIAL") next
      c[t,s]++
      tot[t]++
      tot_status[s]++
      grand++
    }
    END{
      # write detailed
      for(k in c){
        split(k,a,SUBSEP)
        print sample,stage,caller,a[1],a[2],c[k] >> "'"$event_by_type"'"
      }

      # write summary
      if(grand==0){
        print sample,stage,caller,"event_total",0 >> "'"$event_summary"'"
        print sample,stage,caller,"event_full",0 >> "'"$event_summary"'"
        print sample,stage,caller,"event_partial",0 >> "'"$event_summary"'"
        print sample,stage,caller,"mnv_total",0 >> "'"$event_summary"'"
        print sample,stage,caller,"multiallelic_total",0 >> "'"$event_summary"'"
      } else {
        print sample,stage,caller,"event_total",grand >> "'"$event_summary"'"
        print sample,stage,caller,"event_full",tot_status["FULL"]+0 >> "'"$event_summary"'"
        print sample,stage,caller,"event_partial",tot_status["PARTIAL"]+0 >> "'"$event_summary"'"
        print sample,stage,caller,"mnv_total",tot["MNV"]+0 >> "'"$event_summary"'"
        print sample,stage,caller,"multiallelic_total",tot["MULTIALLELIC"]+0 >> "'"$event_summary"'"
      }
    }'
}

echo -e "sample\tstage\tcaller\tEVENT_TYPE\tEVENT_STATUS\tcount" > "$event_by_type"
echo -e "sample\tstage\tcaller\tmetric\tvalue" > "$event_summary"


# ----------------------------
# Process VCFs: extract variants, classify, count
# ----------------------------
echo -e "sample\tstage\tcaller\tvartype\tcount" > "$counts_long"
echo -e "sample\tstage\tcaller\tCHROM\tbin_start\tbin_end\tcount" > "$genome_bins"

log "Processing ${#VCF_SPECS[@]} VCF specifications"

for spec in "${VCF_SPECS[@]}"; do
  stage="$(parse_spec "$spec" "stage")"
  caller="$(parse_spec "$spec" "caller")"
  path="$(parse_spec "$spec" "path")"

  [[ -n "$stage" && -n "$caller" && -n "$path" ]] || {
    echo "ERROR: malformed --vcf spec: $spec" >&2
    echo "Expected format: stage=STAGE caller=CALLER path=PATH" >&2
    exit 2
  }
  
  [[ -f "$path" ]] || {
    echo "ERROR: VCF not found: $path" >&2
    exit 2
  }

  log "Processing: stage=$stage caller=$caller"

  tmp="${TMPDIR}/${stage}.${caller}.sites.tsv"
  
  # Extract CHROM POS REF ALT
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$path" 2>/dev/null \
    | awk 'BEGIN{FS=OFS="\t"} NF==4 && $4!="." && $4!="" {print}' \
    > "$tmp"

  if [[ ! -s "$tmp" ]]; then
    log "  WARNING: No variants found in $path"
    continue
  fi

  # NEW: event-level FULL/PARTIAL counts (only meaningful for merged/event stage)
  if [[ "$stage" == "merged" || "$stage" == "event" ]]; then
    log "  Event-level metrics (EVENT_TYPE/EVENT_STATUS) from $path"
    count_event_status_from_vcf "$path" "$stage" "$caller" || true
  fi

  # Expand multiallelic and classify
  tmp_typed="${tmp}.typed"
  awk 'BEGIN{FS=OFS="\t"}
  {
    n=split($4,alts,",");
    for(i=1;i<=n;i++){
      ref=$3; alt=alts[i];
      if (alt ~ /,/) vt="MULTIALLELIC";
      else if (length(ref)==1 && length(alt)==1) vt="SNP";
      else if (length(ref)==length(alt) && length(ref)>1) vt="MNV";
      else vt="INDEL";
      print $1,$2,$3,alts[i],vt;
    }
  }' "$tmp" > "$tmp_typed"

  if [[ ! -s "$tmp_typed" ]]; then
    log "  WARNING: No valid variants after expansion"
    continue
  fi

  # Count by vartype
  awk -v sample="$SAMPLE" -v stage="$stage" -v caller="$caller" \
    'BEGIN{FS=OFS="\t"} {c[$5]++} END{for(v in c) print sample,stage,caller,v,c[v]}' \
    "$tmp_typed" >> "$counts_long"

  # Genomic bins
  awk -v sample="$SAMPLE" -v stage="$stage" -v caller="$caller" -v B="$BIN_SIZE" \
    'BEGIN{FS=OFS="\t"}
    {
      pos=$2+0
      bin_start = int((pos-1)/B)*B + 1
      bin_end = bin_start + B - 1
      k=$1 OFS bin_start OFS bin_end
      n[k]++
    }
    END{for(k in n){split(k,a,OFS); print sample,stage,caller,a[1],a[2],a[3],n[k]}}' \
    "$tmp_typed" >> "$genome_bins"

  log "  Processed $(wc -l < "$tmp_typed") variants"
done

# Sort genomic bins
sort -t$'\t' -k1,1 -k2,2 -k3,3 -k4,4V -k5,5n "$genome_bins" -o "$genome_bins"

# ----------------------------
# Generate wide counts matrix
# ----------------------------
log "Generating wide variant counts matrix"

awk 'BEGIN{FS=OFS="\t"}
NR>1 {
  key=$1 OFS $2 OFS $3
  keys[key]=1
  vartypes[$4]=1
  data[key,$4]=$5
}
END{
  # Header
  printf "sample\tstage\tcaller"
  n=asorti(vartypes, vt_sorted)
  for(i=1;i<=n;i++) printf "\t%s", vt_sorted[i]
  printf "\ttotal\n"
  
  # Data rows
  for(key in keys){
    split(key,a,OFS)
    printf "%s\t%s\t%s", a[1],a[2],a[3]
    total=0
    for(i=1;i<=n;i++){
      v=data[key,vt_sorted[i]]+0
      printf "\t%d", v
      total+=v
    }
    printf "\t%d\n", total
  }
}' "$counts_long" | sort -t$'\t' -k1,1 -k2,2 -k3,3 > "$counts_wide"

# ----------------------------
# Generate genomic summary (per chromosome)
# ----------------------------
log "Generating genomic summary"

echo -e "sample\tstage\tcaller\tCHROM\ttotal_variants" > "$genome_summary"
awk 'BEGIN{FS=OFS="\t"}
NR>1 {
  key=$1 OFS $2 OFS $3 OFS $4
  n[key]+=$7
}
END{
  for(key in n){
    print key, n[key]
  }
}' "$genome_bins" | sort -t$'\t' -k1,1 -k2,2 -k3,3 -k4,4V >> "$genome_summary"

# ----------------------------
# Process K-pass files
# ----------------------------
# k_pass TSV format: CHROM POS REF ALT callers n_callers pass
# Column 5 = callers (comma-separated list)
# Column 6 = n_callers (count)
# Column 7 = pass (1 if n_callers >= K)

# Additional outputs for k-pass analysis
kcurve_exact="$OUTDIR/k_curve_exact.tsv"
kcurve_exact_vartype="$OUTDIR/k_curve_exact_by_vartype.tsv"
kcurve_exact_caller="$OUTDIR/k_curve_exact_by_caller.tsv"

if [[ ${#KPASS_SPECS[@]} -gt 0 ]]; then
  log "Processing ${#KPASS_SPECS[@]} k-pass specifications"
  
  # Initialize K-curve outputs
  # Cumulative: variants with >= K callers
  echo -e "sample\tK\tpass_variants\ttotal_variants\tpercent_pass" > "$kcurve"
  echo -e "sample\tK\tvartype\tcount" > "$kcurve_vartype"
  
  # Exact: variants with exactly K callers
  echo -e "sample\tn_callers\tcount\tpercent" > "$kcurve_exact"
  echo -e "sample\tn_callers\tvartype\tcount" > "$kcurve_exact_vartype"
  # For n_callers=1, show which caller; for n_callers>1, show "shared"
  echo -e "sample\tn_callers\tcaller_or_shared\tvartype\tcount" > "$kcurve_exact_caller"
  
  # Find the k=1 file for full union (used for all analyses)
  k1_file=""
  
  # First pass: collect all K values and find k1 file
  declare -A K_FILES=()
  for spec in "${KPASS_SPECS[@]}"; do
    k_val="$(parse_spec "$spec" "k")"
    k_path="$(parse_spec "$spec" "path")"
    
    [[ -n "$k_val" && -n "$k_path" ]] || {
      echo "ERROR: malformed --kpass spec: $spec" >&2
      echo "Expected format: k=INT path=FILE" >&2
      exit 2
    }
    
    [[ -f "$k_path" ]] || {
      echo "ERROR: k-pass file not found: $k_path" >&2
      exit 2
    }
    
    K_FILES["$k_val"]="$k_path"
    
    if [[ "$k_val" == "1" ]]; then
      k1_file="$k_path"
    fi
  done
  
  # Process each K value for cumulative K-curve (>= K callers)
  for k_val in "${!K_FILES[@]}"; do
    k_path="${K_FILES[$k_val]}"
    log "Processing k=$k_val from $k_path"
    
    # Count total and passing variants (cumulative: >= K)
    stats=$(awk 'BEGIN{FS="\t"; total=0; pass=0}
      NF>=7 {
        total++
        if($7==1) pass++
      }
      END{
        pct = (total>0) ? 100*pass/total : 0
        printf "%d\t%d\t%.2f", pass, total, pct
      }' "$k_path")
    
    pass_n=$(echo "$stats" | cut -f1)
    total_n=$(echo "$stats" | cut -f2)
    pct=$(echo "$stats" | cut -f3)
    
    echo -e "${SAMPLE}\t${k_val}\t${pass_n}\t${total_n}\t${pct}" >> "$kcurve"
    
    # Count by variant type for passing variants (cumulative)
    awk -v sample="$SAMPLE" -v k="$k_val" 'BEGIN{FS=OFS="\t"}
    NF>=7 && $7==1 {
      ref=$3; alt=$4
      if (alt ~ /,/) vt="MULTIALLELIC"
      else if (length(ref)==1 && length(alt)==1) vt="SNP"
      else if (length(ref)==length(alt) && length(ref)>1) vt="MNV"
      else vt="INDEL"
      counts[vt]++
    }
    END{
      for(vt in counts) print sample, k, vt, counts[vt]
    }' "$k_path" >> "$kcurve_vartype"
    
    log "  k>=$k_val: $pass_n / $total_n variants ($pct%)"
  done
  
  # Sort cumulative k_curve files by K value
  head -1 "$kcurve" > "${kcurve}.tmp"
  tail -n +2 "$kcurve" | sort -t$'\t' -k2,2n >> "${kcurve}.tmp"
  mv "${kcurve}.tmp" "$kcurve"
  
  head -1 "$kcurve_vartype" > "${kcurve_vartype}.tmp"
  tail -n +2 "$kcurve_vartype" | sort -t$'\t' -k2,2n -k3,3 >> "${kcurve_vartype}.tmp"
  mv "${kcurve_vartype}.tmp" "$kcurve_vartype"
  
  # ----------------------------
  # Generate EXACT counts from k1 file (variants with exactly N callers)
  # ----------------------------
  if [[ -n "$k1_file" && -f "$k1_file" ]]; then
    log "Generating exact K counts and caller breakdown from k=1 file"
    
    # Get total variants
    total_variants=$(wc -l < "$k1_file")
    
    # Exact counts: variants with exactly n_callers = 1, 2, 3, ...
    awk -v sample="$SAMPLE" -v total="$total_variants" 'BEGIN{FS=OFS="\t"}
    NF>=6 {
      n = $6+0
      counts[n]++
    }
    END{
      for(n in counts){
        pct = (total>0) ? 100*counts[n]/total : 0
        printf "%s\t%d\t%d\t%.2f\n", sample, n, counts[n], pct
      }
    }' "$k1_file" | sort -t$'\t' -k2,2n >> "$kcurve_exact"
    
    # Exact counts by variant type
    awk -v sample="$SAMPLE" 'BEGIN{FS=OFS="\t"}
    NF>=6 {
      n = $6+0
      ref=$3; alt=$4
      if (alt ~ /,/) vt="MULTIALLELIC"
      else if (length(ref)==1 && length(alt)==1) vt="SNP"
      else if (length(ref)==length(alt) && length(ref)>1) vt="MNV"
      else vt="INDEL"
      counts[n,vt]++
    }
    END{
      for(key in counts){
        split(key, a, SUBSEP)
        print sample, a[1], a[2], counts[key]
      }
    }' "$k1_file" | sort -t$'\t' -k2,2n -k3,3 >> "$kcurve_exact_vartype"
    
    # Exact counts by caller (for unique) or "shared_N" (for N>1)
    # This shows: for variants called by exactly 1 caller, which caller
    #             for variants called by exactly 2 callers, which pair
    #             for variants called by all callers, "all_callers"
    awk -v sample="$SAMPLE" 'BEGIN{FS=OFS="\t"}
    NF>=6 {
      n = $6+0
      callers = $5
      ref=$3; alt=$4
      
      if (alt ~ /,/) vt="MULTIALLELIC"
      else if (length(ref)==1 && length(alt)==1) vt="SNP"
      else if (length(ref)==length(alt) && length(ref)>1) vt="MNV"
      else vt="INDEL"
      
      if(n == 1){
        # Unique to one caller - use caller name
        caller_label = callers
      } else {
        # Shared - use sorted caller combination
        # Sort the callers for consistent labeling
        split(callers, arr, ",")
        n_arr = 0
        for(i in arr) n_arr++
        
        # Simple bubble sort for small arrays
        for(i=1; i<=n_arr; i++){
          for(j=i+1; j<=n_arr; j++){
            if(arr[i] > arr[j]){
              tmp = arr[i]; arr[i] = arr[j]; arr[j] = tmp
            }
          }
        }
        
        caller_label = arr[1]
        for(i=2; i<=n_arr; i++) caller_label = caller_label "+" arr[i]
      }
      
      counts[n, caller_label, vt]++
    }
    END{
      for(key in counts){
        split(key, a, SUBSEP)
        print sample, a[1], a[2], a[3], counts[key]
      }
    }' "$k1_file" | sort -t$'\t' -k2,2n -k3,3 -k4,4 >> "$kcurve_exact_caller"
    
    log "  Exact K counts generated"
    
    # ----------------------------
    # Get the list of callers from the data
    # ----------------------------
    callers_list=$(awk 'BEGIN{FS="\t"}
      NF>=5 {
        split($5, arr, ",")
        for(i in arr) callers[arr[i]]=1
      }
      END{
        # Sort callers alphabetically
        n = asorti(callers, sorted)
        for(i=1; i<=n; i++){
          if(i>1) printf ","
          printf "%s", sorted[i]
        }
      }' "$k1_file")
    
    log "  Detected callers: $callers_list"
    
    # Generate UpSet matrix (wide format with 0/1 per caller)
    awk -v sample="$SAMPLE" -v callers="$callers_list" 'BEGIN{
      FS=OFS="\t"
      n_callers = split(callers, caller_arr, ",")
      for(i=1; i<=n_callers; i++) caller_idx[caller_arr[i]] = i
      
      # Print header
      printf "sample\tCHROM\tPOS\tREF\tALT\tn_callers\tcallers"
      for(i=1; i<=n_callers; i++) printf "\t%s", caller_arr[i]
      printf "\tvartype\n"
    }
    NF>=6 {
      ref=$3; alt=$4
      if (alt ~ /,/) vt="MULTIALLELIC"
      else if (length(ref)==1 && length(alt)==1) vt="SNP"
      else if (length(ref)==length(alt) && length(ref)>1) vt="MNV"
      else vt="INDEL"
      
      # Initialize presence array
      for(i=1; i<=n_callers; i++) present[i] = 0
      
      # Mark callers present for this variant
      split($5, var_callers, ",")
      for(i in var_callers) {
        c = var_callers[i]
        if(c in caller_idx) present[caller_idx[c]] = 1
      }
      
      # Output row
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s", sample, $1, $2, $3, $4, $6, $5
      for(i=1; i<=n_callers; i++) printf "\t%d", present[i]
      printf "\t%s\n", vt
    }' "$k1_file" > "$upset_matrix"
    
    log "  UpSet matrix: $(( $(wc -l < "$upset_matrix") - 1 )) variants"
    
    # Generate caller agreement matrix (pairwise overlaps)
    awk -v callers="$callers_list" 'BEGIN{
      FS=OFS="\t"
      n_callers = split(callers, caller_arr, ",")
      for(i=1; i<=n_callers; i++) caller_idx[caller_arr[i]] = i
    }
    NF>=5 {
      # Parse callers for this variant
      split($5, var_callers, ",")
      
      # Mark presence
      for(i=1; i<=n_callers; i++) present[i] = 0
      for(i in var_callers) {
        c = var_callers[i]
        if(c in caller_idx) present[caller_idx[c]] = 1
      }
      
      # Count single and pairwise
      for(i=1; i<=n_callers; i++){
        if(present[i]==1){
          single[i]++
          for(j=i; j<=n_callers; j++){
            if(present[j]==1) pair[i,j]++
          }
        }
      }
    }
    END{
      # Print header
      printf "caller"
      for(i=1; i<=n_callers; i++) printf "\t%s", caller_arr[i]
      printf "\n"
      
      # Print matrix
      for(i=1; i<=n_callers; i++){
        printf "%s", caller_arr[i]
        for(j=1; j<=n_callers; j++){
          if(i==j) v=single[i]+0
          else if(i<j) v=pair[i,j]+0
          else v=pair[j,i]+0
          printf "\t%d", v
        }
        printf "\n"
      }
    }' "$k1_file" > "$caller_agreement"
    
    log "  Caller agreement matrix generated"
    
    # Generate caller contribution summary
    echo -e "sample\tcaller\tunique\tshared_2\tshared_3\tshared_all\ttotal" > "$caller_contrib"
    
    awk -v sample="$SAMPLE" -v callers="$callers_list" 'BEGIN{
      FS=OFS="\t"
      n_callers = split(callers, caller_arr, ",")
      for(i=1; i<=n_callers; i++) caller_idx[caller_arr[i]] = i
    }
    NF>=6 {
      n = $6+0  # n_callers for this variant
      
      # Parse callers
      split($5, var_callers, ",")
      for(i=1; i<=n_callers; i++) present[i] = 0
      for(i in var_callers) {
        c = var_callers[i]
        if(c in caller_idx) present[caller_idx[c]] = 1
      }
      
      # Categorize by agreement level
      for(i=1; i<=n_callers; i++){
        if(present[i]==1){
          total[i]++
          if(n==1) unique[i]++
          else if(n==2) shared2[i]++
          else if(n==3) shared3[i]++
          else shared_all[i]++
        }
      }
    }
    END{
      for(i=1; i<=n_callers; i++){
        printf "%s\t%s\t%d\t%d\t%d\t%d\t%d\n",
          sample, caller_arr[i], unique[i]+0, shared2[i]+0, shared3[i]+0, shared_all[i]+0, total[i]+0
      }
    }' "$k1_file" >> "$caller_contrib"
    
    log "  Caller contribution generated"
    
  else
    log "WARNING: No k=1 file found, skipping detailed analysis"
    echo -e "sample\tn_callers\tcount\tpercent" > "$kcurve_exact"
    echo -e "sample\tn_callers\tvartype\tcount" > "$kcurve_exact_vartype"
    echo -e "sample\tn_callers\tcaller_or_shared\tvartype\tcount" > "$kcurve_exact_caller"
    echo -e "sample\tCHROM\tPOS\tREF\tALT\tn_callers\tcallers\tvartype" > "$upset_matrix"
    echo -e "caller" > "$caller_agreement"
    echo -e "sample\tcaller\tunique\tshared_2\tshared_3\tshared_all\ttotal" > "$caller_contrib"
  fi

else
  log "No --kpass files provided, skipping K-curve and caller stats"
  echo -e "sample\tK\tpass_variants\ttotal_variants\tpercent_pass" > "$kcurve"
  echo -e "sample\tK\tvartype\tcount" > "$kcurve_vartype"
  echo -e "sample\tn_callers\tcount\tpercent" > "$kcurve_exact"
  echo -e "sample\tn_callers\tvartype\tcount" > "$kcurve_exact_vartype"
  echo -e "sample\tn_callers\tcaller_or_shared\tvartype\tcount" > "$kcurve_exact_caller"
  echo -e "sample\tCHROM\tPOS\tREF\tALT\tn_callers\tcallers\tvartype" > "$upset_matrix"
  echo -e "caller" > "$caller_agreement"
  echo -e "sample\tcaller\tunique\tshared_2\tshared_3\tshared_all\ttotal" > "$caller_contrib"
fi


# ----------------------------
# Summary statistics
# ----------------------------
log "Generating summary statistics"

summary_out="$OUTDIR/summary_stats.tsv"
echo -e "sample\tmetric\tvalue" > "$summary_out"

# Total variants per stage
awk -v sample="$SAMPLE" 'BEGIN{FS=OFS="\t"}
NR>1 {
  key=$2 OFS $3
  n[key]+=$5
}
END{
  for(key in n){
    print sample, "total_variants:" key, n[key]
  }
}' "$counts_long" >> "$summary_out"

# Variants by type (total across all)
awk -v sample="$SAMPLE" 'BEGIN{FS=OFS="\t"}
NR>1 { n[$4]+=$5 }
END{
  for(vt in n) print sample, "vartype:" vt, n[vt]
}' "$counts_long" >> "$summary_out"

# K-curve summary if available
if [[ -s "$kcurve" ]] && [[ $(wc -l < "$kcurve") -gt 1 ]]; then
  awk -v sample="$SAMPLE" 'BEGIN{FS=OFS="\t"}
  NR>1 {
    print sample, "k" $2 "_pass", $3
    print sample, "k" $2 "_pct", $5
  }' "$kcurve" >> "$summary_out"
fi

# ----------------------------
# Anchor VCF DP summary
# ----------------------------

echo -e "sample\tstage\tcaller\tdp_source\tmean_dp\tmedian_dp\tmin_dp\tmax_dp\tn_sites" > "$vcf_dp_summary"

# Parse spec helper
parse_spec() {
  local spec="$1"
  local key="$2"
  # Use { grep ... || true; } to prevent set -e from exiting
  { echo "$spec" | tr ' ' '\n' | grep "^${key}=" | cut -d= -f2-; } || true
}

summarize_vcf_dp_from_anchor() {
  local spec="$1"

  local stage caller vcf
  stage="$(parse_spec "$spec" "stage")"
  caller="$(parse_spec "$spec" "caller")"
  vcf="$(parse_spec "$spec" "path")"

  [[ -n "$stage" && -n "$caller" && -n "$vcf" ]] || {
    echo "ERROR: malformed --anchor-vcf spec: $spec" >&2
    echo "Expected: stage=STAGE caller=CALLER path=VCF.gz" >&2
    exit 2
  }
  [[ -f "$vcf" ]] || { echo "ERROR: anchor VCF not found: $vcf" >&2; exit 2; }

  log "Anchor VCF DP summary: stage=$stage caller=$caller vcf=$vcf"

  local hdr dp_source dp_file sorted
  
  # Get header and remove any carriage returns (Windows line endings)
  hdr="$(bcftools view -h "$vcf" 2>/dev/null | tr -d '\r')" || hdr=""

  dp_file="$TMPDIR/anchor.dp.txt"
  sorted="$TMPDIR/anchor.dp.sorted.txt"
  : > "$dp_file"

  # Helper to write one line summary (or NA if no sites)
  _emit_dp_summary() {
    local src="$1"
    local f="$2"
    local n mean med min max mid m1 m2

    if [[ ! -s "$f" ]]; then
      echo -e "${SAMPLE}\t${stage}\t${caller}\t${src}\tNA\tNA\tNA\tNA\t0" >> "$vcf_dp_summary"
      return
    fi

    # Sort numeric
    sort -n "$f" > "$sorted"
    n=$(wc -l < "$sorted")
    n="${n// /}"  # Remove any whitespace

    # Calculate mean / min / max in one pass using awk
    local stats
    stats=$(awk 'BEGIN{s=0; min=""; max=""}
         {
           x=$1+0
           s+=x
           if(min=="" || x<min) min=x
           if(max=="" || x>max) max=x
         }
         END{
           if(NR==0){print "NA NA NA"; exit}
           printf "%.6f %s %s\n", s/NR, min, max
         }' "$sorted")
    
    mean=$(echo "$stats" | awk '{print $1}')
    min=$(echo "$stats" | awk '{print $2}')
    max=$(echo "$stats" | awk '{print $3}')

    # Calculate median
    if (( n % 2 == 1 )); then
      mid=$(( (n+1)/2 ))
      med=$(awk -v m="$mid" 'NR==m{print $1; exit}' "$sorted")
    else
      m1=$(( n/2 ))
      m2=$(( m1+1 ))
      med=$(awk -v a="$m1" -v b="$m2" 'NR==a{v=$1} NR==b{print (v+$1)/2; exit}' "$sorted")
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\n" \
      "$SAMPLE" "$stage" "$caller" "$src" "$mean" "$med" "$min" "$max" "$n" \
      >> "$vcf_dp_summary"
  }

  # Build sample argument - only add if DP_SAMPLE is set
  local sample_args=""
  if [[ -n "${DP_SAMPLE:-}" ]]; then
    sample_args="-s ${DP_SAMPLE}"
  fi

  # Check for DP fields using grep -c (returns count, always succeeds)
  # Then check if count > 0
  local format_dp_count format_ad_count info_dp_count
  format_dp_count=$(echo "$hdr" | grep -c '^##FORMAT=<ID=DP[,>]' 2>/dev/null || echo "0")
  format_ad_count=$(echo "$hdr" | grep -c '^##FORMAT=<ID=AD[,>]' 2>/dev/null || echo "0")
  info_dp_count=$(echo "$hdr" | grep -c '^##INFO=<ID=DP[,>]' 2>/dev/null || echo "0")

  # Priority 1: FORMAT/DP
  if [[ "$format_dp_count" -gt 0 ]]; then
    dp_source="FORMAT/DP"
    log "  Using $dp_source"
    if [[ -n "$sample_args" ]]; then
      { bcftools query $sample_args -f '[%DP]\n' "$vcf" 2>/dev/null \
        | awk '$1!="." && $1!="" && $1~/^[0-9]+$/ {print $1}' > "$dp_file"; } || true
    else
      { bcftools query -f '[%DP]\n' "$vcf" 2>/dev/null \
        | awk '$1!="." && $1!="" && $1~/^[0-9]+$/ {print $1}' > "$dp_file"; } || true
    fi
    _emit_dp_summary "$dp_source" "$dp_file"
    return
  fi

  # Priority 2: FORMAT/AD -> sum(AD) as DP proxy
  if [[ "$format_ad_count" -gt 0 ]]; then
    dp_source="FORMAT/AD_sum"
    log "  Using $dp_source"
    if [[ -n "$sample_args" ]]; then
      { bcftools query $sample_args -f '[%AD]\n' "$vcf" 2>/dev/null \
        | awk 'BEGIN{FS=","} 
               $1!="." && $1!="" {
                 dp=0
                 for(i=1;i<=NF;i++) if($i~/^[0-9]+$/) dp+=$i
                 if(dp>0) print dp
               }' > "$dp_file"; } || true
    else
      { bcftools query -f '[%AD]\n' "$vcf" 2>/dev/null \
        | awk 'BEGIN{FS=","} 
               $1!="." && $1!="" {
                 dp=0
                 for(i=1;i<=NF;i++) if($i~/^[0-9]+$/) dp+=$i
                 if(dp>0) print dp
               }' > "$dp_file"; } || true
    fi
    _emit_dp_summary "$dp_source" "$dp_file"
    return
  fi

  # Priority 3: INFO/DP
  if [[ "$info_dp_count" -gt 0 ]]; then
    dp_source="INFO/DP"
    log "  Using $dp_source"
    { bcftools query -f '%INFO/DP\n' "$vcf" 2>/dev/null \
      | awk '$1!="." && $1!="" && $1~/^[0-9]+$/ {print $1}' > "$dp_file"; } || true
    _emit_dp_summary "$dp_source" "$dp_file"
    return
  fi

  # Nothing found
  log "  WARNING: No DP-like fields found in anchor VCF header (FORMAT/DP, FORMAT/AD, INFO/DP). Writing NA."
  echo -e "${SAMPLE}\t${stage}\t${caller}\tNONE\tNA\tNA\tNA\tNA\t0" >> "$vcf_dp_summary"
}


if [[ -n "${ANCHOR_SPEC}" ]]; then
  summarize_vcf_dp_from_anchor "$ANCHOR_SPEC"
else
  log "No --anchor-vcf provided; skipping VCF DP summary"
fi
# Anchor VCF DP summary into summary_stats.tsv (if present)
if [[ -s "$vcf_dp_summary" ]] && [[ $(wc -l < "$vcf_dp_summary") -gt 1 ]]; then
  awk -v sample="$SAMPLE" 'BEGIN{FS=OFS="\t"}
    NR==2 {
      stage=$2; caller=$3; src=$4
      print sample, "anchor_vcf_dp_source:" stage ":" caller, src
      print sample, "anchor_vcf_mean_dp:" stage ":" caller, $5
      print sample, "anchor_vcf_median_dp:" stage ":" caller, $6
      print sample, "anchor_vcf_min_dp:" stage ":" caller, $7
      print sample, "anchor_vcf_max_dp:" stage ":" caller, $8
      print sample, "anchor_vcf_n_sites:" stage ":" caller, $9
    }' "$vcf_dp_summary" >> "$summary_out"
fi

# ----------------------------
# summary
# ----------------------------

log "Done. Outputs in: $OUTDIR"
log "Files generated:"
for f in "$counts_long" "$counts_wide" "$upset_matrix" \
         "$caller_agreement" "$caller_contrib" \
         "$kcurve" "$kcurve_vartype" \
         "$kcurve_exact" "$kcurve_exact_vartype" "$kcurve_exact_caller" \
         "$genome_bins" "$genome_summary" "$summary_out"; do
  if [[ -f "$f" ]]; then
    lines=$(wc -l < "$f")
    log "  $(basename "$f"): $lines lines"
  fi
done
