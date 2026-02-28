#!/usr/bin/env bash
set -euo pipefail

# consensus_k_of_n_enhanced.sh
#
# Build a K-of-N consensus across multiple caller VCFs.
# Enhanced with --k-range to output mutations at all K thresholds.

usage() {
  cat <<'EOF'
consensus_k_of_n_enhanced.sh

Modes:
  --mode union     Build union.calls.tsv + union.anno.tsv only (no pass/emit)
  --mode kpass     Derive pass TSV for a given --k from an existing union.anno.tsv
  --mode emit      Emit consensus VCF using an existing pass TSV (no union recompute)
  --mode all       union -> kpass -> emit (default)
  --mode compare   Build union, then output pass TSVs + VCFs for all K in range

Required (mode=all|union|compare):
  --ref PATH
  --vcf caller=PATH            Repeatable. Example:
                                 --vcf mutect2=mutect.vcf.gz
                                 --vcf muse2=muse.vcf.gz
                                 --vcf strelka2=strelka.vcf.gz
  --workdir PATH               Persistent work directory (NOT deleted)

Required (mode=all|kpass|emit):
  --k INT                      K threshold (e.g. 2 for 2-of-3)

Required (mode=compare):
  --k-range MIN:MAX            Range of K values to test (e.g., 1:3 for K=1,2,3)
                               Default for compare mode: 1:N where N = number of callers

Required (mode=all|emit|compare):
  --out-prefix STR             Output prefix (filename prefix; can include path)

Optional:
  --normalize true|false       Default: false
  --normalize-script PATH      Default: ./normalize_vcf.sh (or in PATH)
  --threads INT                Default: 4
  --anchor CALLER              Header source + preference seed (default: first --vcf)
  --emit-vcf true|false        Default: true
  --keep-all-annot true|false  Default: true
  --log PATH                   Default: <workdir>/consensus_k_of_n.log

Inputs/Outputs wiring (for Nextflow fan-out):
  --union-anno PATH            Use this union.anno.tsv instead of <workdir>/union.anno.tsv
  --pass-in PATH               Use this pass TSV for mode=emit (default: <workdir>/consensus.pass.tsv)
  --pass-out PATH              Write pass TSV here for mode=kpass (default: <workdir>/k<k>.pass.tsv)

Workdir outputs:
  union.calls.tsv              CHROM POS REF ALT CALLER
  union.anno.tsv               CHROM POS REF ALT callers n_callers
  consensus.pass.tsv           subset where pass=1 for current k (mode=all)
  k<k>.pass.tsv                subset for selected k (mode=kpass/compare)
  comparison.summary.tsv       variant counts at each K (mode=compare)
  comparison.detail.tsv        per-variant K support detail (mode=compare)

Main outputs:
  <out-prefix>.consensus.vcf        consensus calls as VCF (if emit-vcf=true)
  <out-prefix>.k<K>.vcf             consensus VCF for each K (mode=compare)

EOF
}

# ---------------- args ----------------
MODE="all"

K=""
K_RANGE=""
REF=""
OUT_PREFIX=""
WORKDIR=""
NORMALIZE="false"
NORMALIZE_SCRIPT="normalize_vcf.sh"
THREADS=4
ANCHOR=""
EMIT_VCF="true"
KEEP_ALL="true"
LOG=""

UNION_ANNO_IN=""
PASS_IN=""
PASS_OUT=""

declare -a VCF_SPECS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode) MODE="$2"; shift 2;;

    --k) K="$2"; shift 2;;
    --k-range) K_RANGE="$2"; shift 2;;
    --ref) REF="$2"; shift 2;;
    --vcf) VCF_SPECS+=("$2"); shift 2;;
    --out-prefix) OUT_PREFIX="$2"; shift 2;;
    --workdir) WORKDIR="$2"; shift 2;;

    --normalize) NORMALIZE="$2"; shift 2;;
    --normalize-script) NORMALIZE_SCRIPT="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --anchor) ANCHOR="$2"; shift 2;;
    --emit-vcf) EMIT_VCF="$2"; shift 2;;
    --keep-all-annot) KEEP_ALL="$2"; shift 2;;
    --log) LOG="$2"; shift 2;;

    --union-anno) UNION_ANNO_IN="$2"; shift 2;;
    --pass-in) PASS_IN="$2"; shift 2;;
    --pass-out) PASS_OUT="$2"; shift 2;;

    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1" >&2; usage; exit 2;;
  esac
done

# ---------------- helpers ----------------
ts() { date +"%F %T"; }
log() { echo "[$(ts)] $*" | tee -a "$LOG" >&2; }

need() { [[ -n "${1:-}" ]] || { echo "ERROR: $2" >&2; usage; exit 2; }; }

# ---------------- basic checks by mode ----------------
need "$WORKDIR" "missing required: --workdir"

mkdir -p "$WORKDIR"
if [[ -z "$LOG" ]]; then
  LOG="$WORKDIR/consensus_k_of_n.log"
fi

log "Starting consensus_k_of_n (mode=$MODE)"
log "threads=$THREADS  normalize=$NORMALIZE"
log "workdir=$WORKDIR"

command -v bcftools >/dev/null 2>&1 || { echo "ERROR: bcftools not found in PATH" >&2; exit 127; }
command -v tabix   >/dev/null 2>&1 || { echo "ERROR: tabix not found in PATH" >&2; exit 127; }
command -v awk     >/dev/null 2>&1 || { echo "ERROR: awk not found in PATH" >&2; exit 127; }
command -v bgzip   >/dev/null 2>&1 || { echo "ERROR: bgzip not found in PATH" >&2; exit 127; }

# Validate ref needed for union/emit/all/compare (or if normalize true)
if [[ "$MODE" == "union" || "$MODE" == "emit" || "$MODE" == "all" || "$MODE" == "compare" || "$NORMALIZE" == "true" ]]; then
  need "$REF" "missing required: --ref"
  [[ -f "$REF" ]] || { echo "ERROR: ref fasta not found: $REF" >&2; exit 2; }
  [[ -f "${REF}.fai" ]] || { echo "ERROR: missing fasta index: ${REF}.fai" >&2; exit 2; }
  log "ref=$REF"
fi

# Validate K for kpass/emit/all
if [[ "$MODE" == "kpass" || "$MODE" == "emit" || "$MODE" == "all" ]]; then
  need "$K" "missing required: --k"
  [[ "$K" =~ ^[0-9]+$ ]] || { echo "ERROR: --k must be a positive integer, got: $K" >&2; exit 2; }
  [[ "$K" -gt 0 ]] || { echo "ERROR: --k must be > 0" >&2; exit 2; }
  log "K=$K"
fi

# out-prefix required for emit/all/compare
if [[ "$MODE" == "emit" || "$MODE" == "all" || "$MODE" == "compare" ]]; then
  need "$OUT_PREFIX" "missing required: --out-prefix (mode=$MODE)"
  log "out-prefix=$OUT_PREFIX"
fi

# normalize script check only if used
if [[ "$NORMALIZE" == "true" ]]; then
  if command -v "$NORMALIZE_SCRIPT" >/dev/null 2>&1; then
    NORMALIZE_SCRIPT="$(command -v "$NORMALIZE_SCRIPT")"
  fi
  [[ -x "$NORMALIZE_SCRIPT" ]] || { echo "ERROR: normalize-script not executable: $NORMALIZE_SCRIPT" >&2; exit 2; }
  log "normalize-script=$NORMALIZE_SCRIPT"
fi

# Default file locations
union_tsv="$WORKDIR/union.calls.tsv"
anno_tsv_default="$WORKDIR/union.anno.tsv"
anno_tsv="${UNION_ANNO_IN:-$anno_tsv_default}"

pass_tsv_default="$WORKDIR/consensus.pass.tsv"
pass_tsv_k="$WORKDIR/k${K}.pass.tsv"

# For wiring
if [[ -z "$PASS_IN" ]]; then PASS_IN="$pass_tsv_default"; fi
if [[ -z "$PASS_OUT" ]]; then PASS_OUT="$pass_tsv_k"; fi

# ---------------- parse VCF specs (only needed for union/emit/all/compare) ----------------
declare -A CALLER2VCF=()
declare -A CALLER2WORKVCF=()
declare -a CALLERS=()
declare -a PREF_CALLERS=()

parse_vcfs() {
  [[ ${#VCF_SPECS[@]} -ge 2 ]] || { echo "ERROR: provide at least 2 --vcf caller=path specs" >&2; exit 2; }

  for spec in "${VCF_SPECS[@]}"; do
    [[ "$spec" == *=* ]] || { echo "ERROR: --vcf must be caller=path, got: $spec" >&2; exit 2; }
    local caller="${spec%%=*}"
    local vcf="${spec#*=}"
    [[ -n "$caller" && -n "$vcf" ]] || { echo "ERROR: malformed --vcf $spec" >&2; exit 2; }
    [[ -f "$vcf" ]] || { echo "ERROR: VCF not found: $vcf" >&2; exit 2; }
    CALLER2VCF["$caller"]="$vcf"
    CALLERS+=("$caller")
  done

  # default anchor: first caller
  if [[ -z "$ANCHOR" ]]; then
    ANCHOR="${CALLERS[0]}"
  fi
  [[ -n "${CALLER2VCF[$ANCHOR]:-}" ]] || { echo "ERROR: anchor caller not found among --vcf: $ANCHOR" >&2; exit 2; }

  # Build preference list anchor first
  PREF_CALLERS=("$ANCHOR")
  for c in "${CALLERS[@]}"; do
    [[ "$c" == "$ANCHOR" ]] && continue
    PREF_CALLERS+=("$c")
  done

  log "Callers: ${CALLERS[*]}"
  log "Anchor (header + preference only): $ANCHOR"
  log "Preference order for duplicate keys: ${PREF_CALLERS[*]}"
}

prepare_per_caller_vcfs() {
  for caller in "${CALLERS[@]}"; do
    local in_vcf="${CALLER2VCF[$caller]}"
    local out_vcf="$WORKDIR/${caller}.input.vcf.gz"

    if [[ "$in_vcf" == *.vcf.gz ]]; then
      cp -f "$in_vcf" "$out_vcf"
      if [[ -f "${in_vcf}.tbi" ]]; then
        cp -f "${in_vcf}.tbi" "${out_vcf}.tbi"
      elif [[ -f "${in_vcf}.csi" ]]; then
        tabix -f -p vcf "$out_vcf" >>"$LOG" 2>&1
      else
        tabix -f -p vcf "$out_vcf" >>"$LOG" 2>&1
      fi
    else
      bcftools view --threads "$THREADS" -Oz -o "$out_vcf" "$in_vcf" >>"$LOG" 2>&1
      tabix -f -p vcf "$out_vcf" >>"$LOG" 2>&1
    fi

    if [[ "$NORMALIZE" == "true" ]]; then
      local norm_out="$WORKDIR/${caller}.norm.vcf.gz"
      log "Normalizing $caller: $in_vcf -> $norm_out"
      "$NORMALIZE_SCRIPT" \
        --in-vcf "$out_vcf" \
        --ref "$REF" \
        --out "$norm_out" \
        --threads "$THREADS" >>"$LOG" 2>&1
      CALLER2WORKVCF["$caller"]="$norm_out"
    else
      CALLER2WORKVCF["$caller"]="$out_vcf"
    fi
  done
}

build_union_and_anno() {
  : > "$union_tsv"
  log "Building union.calls.tsv -> $union_tsv"

  for caller in "${CALLERS[@]}"; do
    local vcf="${CALLER2WORKVCF[$caller]}"
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$vcf" \
      | awk -v c="$caller" 'BEGIN{FS=OFS="\t"} NF==4 && $4!="." { print $1,$2,$3,$4,c }' \
      >> "$union_tsv"
  done

  [[ -s "$union_tsv" ]] || { log "ERROR: union.calls.tsv is empty. Upstream VCFs may contain no variants."; exit 1; }

  log "union.calls.tsv lines: $(wc -l < "$union_tsv")"
  log "union.calls.tsv callers: $(cut -f5 "$union_tsv" | sort | uniq -c | tr '\n' '; ')"

  log "Annotating K-of-N (union.anno.tsv) -> $anno_tsv_default"

  awk '
  BEGIN { FS=OFS="\t" }
  NF>=5 {
    key = $1 FS $2 FS $3 FS $4
    caller = $5
    pair = key FS caller

    if (!(pair in seen)) {
      seen[pair] = 1
      n[key]++
      if (callers_list[key] == "")
        callers_list[key] = caller
      else
        callers_list[key] = callers_list[key] "," caller
    }
  }
  END {
    for (key in n) {
      split(key, a, FS)
      print a[1], a[2], a[3], a[4], callers_list[key], n[key]
    }
  }
  ' "$union_tsv" | sort -k1,1V -k2,2n > "$anno_tsv_default"

  [[ -s "$anno_tsv_default" ]] || { log "ERROR: union.anno.tsv is empty; annotation failed."; exit 1; }
  log "union.anno.tsv lines: $(wc -l < "$anno_tsv_default")"
}

kpass_from_union_anno() {
  local k="$1"
  local in_anno="$2"
  local out_pass="$3"

  [[ -s "$in_anno" ]] || { log "ERROR: union.anno.tsv missing/empty: $in_anno"; exit 2; }

  log "Deriving pass TSV from union.anno.tsv (k=$k) -> $out_pass"
  awk -v K="$k" 'BEGIN{FS=OFS="\t"}
    NF>=6{
      n=$6+0
      pass=(n>=K?1:0)
      print $1,$2,$3,$4,$5,$6,pass
    }
  ' "$in_anno" > "$out_pass"

  local pass_n
  pass_n=$(awk 'BEGIN{FS="\t"} $7==1{n++} END{print n+0}' "$out_pass")
  log "Pass variants (k=$k): $pass_n"
  echo "$pass_n"
}

emit_consensus_vcf_from_pass() {
  local pass_file="$1"
  local out_vcf="$2"

  log "Emitting minimal consensus VCF -> $out_vcf"
  local anchor_vcf="${CALLER2WORKVCF[$ANCHOR]}"

  local pass_count
  pass_count=$(awk 'BEGIN{FS="\t"} NF>=7 && $7==1{n++} END{print n+0}' "$pass_file")

  if [[ "$pass_count" -eq 0 ]]; then
    log "WARNING: No variants passed threshold. Writing header-only VCF."
    bcftools view -h "$anchor_vcf" -Ov > "$out_vcf"
    log "Wrote empty consensus VCF: $out_vcf"
    return 0
  fi

  local header_vcf="$WORKDIR/consensus.header.vcf"
  local body_vcf="$WORKDIR/consensus.body.vcf"
  local merged_unsorted="$WORKDIR/consensus.merged.unsorted.vcf"

  # Merge headers from all callers to capture all INFO/FORMAT definitions
  log "Merging headers from all callers to capture all field definitions"
  local all_headers="$WORKDIR/all.headers.txt"
  : > "$all_headers"
  for caller in "${CALLERS[@]}"; do
    bcftools view -h "${CALLER2WORKVCF[$caller]}" >> "$all_headers"
  done

  local header_base="$WORKDIR/header.base.vcf"
  local header_extra="$WORKDIR/header.extra.txt"
  bcftools view -h "$anchor_vcf" > "$header_base"

  : > "$header_extra"
  grep '^##INFO='   "$all_headers" | sort -u >> "$header_extra" || true
  grep '^##FORMAT=' "$all_headers" | sort -u >> "$header_extra" || true
  grep '^##FILTER=' "$all_headers" | sort -u >> "$header_extra" || true
  grep '^##contig=' "$all_headers" | sort -u >> "$header_extra" || true

  awk '
    BEGIN {
      while ((getline line < "'"$header_extra"'") > 0) {
        if (line ~ /^##INFO=<ID=/) {
          match(line, /ID=([^,>]+)/, arr); info_defs[arr[1]] = line
        } else if (line ~ /^##FORMAT=<ID=/) {
          match(line, /ID=([^,>]+)/, arr); format_defs[arr[1]] = line
        } else if (line ~ /^##FILTER=<ID=/) {
          match(line, /ID=([^,>]+)/, arr); filter_defs[arr[1]] = line
        } else if (line ~ /^##contig=/) {
          match(line, /ID=([^,>]+)/, arr); contig_defs[arr[1]] = line
        }
      }
      close("'"$header_extra"'")
    }
    /^##fileformat=/ { print; next }
    /^##FILTER=<ID=/ { match($0,/ID=([^,>]+)/,arr); filter_defs[arr[1]]=$0; next }
    /^##INFO=<ID=/   { match($0,/ID=([^,>]+)/,arr); info_defs[arr[1]]=$0; next }
    /^##FORMAT=<ID=/ { match($0,/ID=([^,>]+)/,arr); format_defs[arr[1]]=$0; next }
    /^##contig=/     { match($0,/ID=([^,>]+)/,arr); contig_defs[arr[1]]=$0; next }
    /^#CHROM/ {
      for (id in filter_defs) print filter_defs[id]
      for (id in info_defs)   print info_defs[id]
      for (id in format_defs) print format_defs[id]
      for (id in contig_defs) print contig_defs[id]
      print "##INFO=<ID=CALLERS,Number=.,Type=String,Description=\"List of callers that called this variant\">"
      print "##INFO=<ID=NCALLERS,Number=1,Type=Integer,Description=\"Number of callers that called this variant\">"
      print "##INFO=<ID=SOURCE_CALLER,Number=1,Type=String,Description=\"Caller from which this variant record was taken\">"
      print
      next
    }
    { print }
  ' "$header_base" > "$header_vcf"

  : > "$body_vcf"

  log "Extracting passing variants from callers (preference order: ${PREF_CALLERS[*]})"
  for caller in "${PREF_CALLERS[@]}"; do
    bcftools view -H "${CALLER2WORKVCF[$caller]}" -Ov --threads "$THREADS" \
      | awk -v c="$caller" 'BEGIN{OFS="\t"} {print c,$0}' > "$WORKDIR/${caller}.body.tsv"
  done

  # preference string for awk
  local pref_str=""
  for i in "${!PREF_CALLERS[@]}"; do
    pref_str+="${PREF_CALLERS[$i]}"
    [[ $i -lt $((${#PREF_CALLERS[@]} - 1)) ]] && pref_str+=","
  done

  local all_bodies="$WORKDIR/all.bodies.tsv"
  : > "$all_bodies"
  for caller in "${PREF_CALLERS[@]}"; do
    cat "$WORKDIR/${caller}.body.tsv" >> "$all_bodies"
  done

  log "Building consensus VCF body by matching pass variants to caller records"

  awk -v pref="$pref_str" '
    BEGIN {
      FS = OFS = "\t"
      n_pref = split(pref, pref_arr, ",")
      for (i = 1; i <= n_pref; i++) pref_rank[pref_arr[i]] = i
    }
    # load pass variants
    FNR == NR {
      if (NF>=7 && $7 == 1) {
        key = $1 OFS $2 OFS $3 OFS $4
        pass[key] = 1
        pass_callers[key] = $5
        pass_ncallers[key] = $6
      }
      next
    }
    # process bodies
    {
      caller = $1
      key = $2 OFS $3 OFS $5 OFS $6
      if (!(key in pass)) next

      has_event = ($9 ~ /(^|;)EVENT_ID=/) ? 1 : 0

      if (!(key in best_caller)) {
        best_caller[key] = caller
        best_line[key] = $0
        best_has_event[key] = has_event
      } else {
        curr_rank = pref_rank[caller]
        best_rank = pref_rank[best_caller[key]]
        replace = 0
        if (has_event && !best_has_event[key]) replace = 1
        else if (has_event == best_has_event[key]) {
          if (curr_rank < best_rank) replace = 1
        }
        if (replace) {
          best_caller[key] = caller
          best_line[key] = $0
          best_has_event[key] = has_event
        }
      }
    }
    END {
      for (key in best_caller) {
        line = best_line[key]
        caller = best_caller[key]
        split(line, f, OFS)
        info = f[9]
        callers_list = pass_callers[key]
        ncall = pass_ncallers[key]
        if (info == "." || info == "") info = "CALLERS=" callers_list ";NCALLERS=" ncall ";SOURCE_CALLER=" caller
        else info = info ";CALLERS=" callers_list ";NCALLERS=" ncall ";SOURCE_CALLER=" caller
        f[9] = info

        out = f[2]
        for (i = 3; i <= length(f); i++) out = out OFS f[i]
        print out
      }
    }
  ' "$pass_file" "$all_bodies" > "$body_vcf"

  cat "$header_vcf" "$body_vcf" > "$merged_unsorted"

  local body_lines
  body_lines=$(wc -l < "$body_vcf" || echo 0)

  if [[ "$body_lines" -eq 0 ]]; then
    log "WARNING: No matched variants; writing header-only VCF."
    cp "$header_vcf" "$out_vcf"
    return 0
  fi

  grep '^#' "$merged_unsorted" > "$out_vcf"
  grep -v '^#' "$merged_unsorted" | sort -k1,1V -k2,2n >> "$out_vcf"

  log "Wrote: $out_vcf (${body_lines} variants)"

  for caller in "${PREF_CALLERS[@]}"; do
    rm -f "$WORKDIR/${caller}.body.tsv"
  done
}

# ---------------- compare mode functions ----------------
run_compare_mode() {
  # Must parse VCFs first to know number of callers
  parse_vcfs
  prepare_per_caller_vcfs

  local n_callers=${#CALLERS[@]}
  local k_min k_max

  # Parse K_RANGE or default to 1:N
  if [[ -n "$K_RANGE" ]]; then
    if [[ "$K_RANGE" =~ ^([0-9]+):([0-9]+)$ ]]; then
      k_min="${BASH_REMATCH[1]}"
      k_max="${BASH_REMATCH[2]}"
    else
      log "ERROR: --k-range must be MIN:MAX (e.g., 1:3), got: $K_RANGE"
      exit 2
    fi
  else
    k_min=1
    k_max=$n_callers
  fi

  [[ "$k_min" -ge 1 ]] || { log "ERROR: k_min must be >= 1"; exit 2; }
  [[ "$k_max" -le "$n_callers" ]] || { log "ERROR: k_max ($k_max) cannot exceed number of callers ($n_callers)"; exit 2; }
  [[ "$k_min" -le "$k_max" ]] || { log "ERROR: k_min ($k_min) must be <= k_max ($k_max)"; exit 2; }

  log "Compare mode: K range = $k_min to $k_max (from $n_callers callers)"

  # Build union
  build_union_and_anno

  # Summary file
  local summary_tsv="$WORKDIR/comparison.summary.tsv"
  echo -e "K\tpass_variants\ttotal_variants\tpercent_pass" > "$summary_tsv"

  local total_variants
  total_variants=$(wc -l < "$anno_tsv_default")

  # Generate pass TSV and VCF for each K (from high to low for logical ordering)
  for (( k=k_max; k>=k_min; k-- )); do
    local pass_tsv="$WORKDIR/k${k}.pass.tsv"
    local pass_count
    pass_count=$(kpass_from_union_anno "$k" "$anno_tsv_default" "$pass_tsv")

    local pct
    if [[ "$total_variants" -gt 0 ]]; then
      pct=$(awk "BEGIN {printf \"%.2f\", 100*$pass_count/$total_variants}")
    else
      pct="0.00"
    fi
    echo -e "${k}\t${pass_count}\t${total_variants}\t${pct}" >> "$summary_tsv"

    if [[ "$EMIT_VCF" == "true" ]]; then
      local out_vcf="${OUT_PREFIX}.k${k}.vcf"
      emit_consensus_vcf_from_pass "$pass_tsv" "$out_vcf"
    fi
  done

  log "Comparison summary written to: $summary_tsv"
  cat "$summary_tsv" | tee -a "$LOG"

  # Generate detail TSV showing which callers called each variant
  local detail_tsv="$WORKDIR/comparison.detail.tsv"
  log "Generating detailed comparison: $detail_tsv"

  # Header: CHROM POS REF ALT n_callers callers caller1 caller2 ... callerN
  local caller_header=""
  for c in "${CALLERS[@]}"; do
    caller_header+="\t${c}"
  done
  echo -e "CHROM\tPOS\tREF\tALT\tn_callers\tcallers${caller_header}" > "$detail_tsv"

  # Build caller presence columns
  local caller_list_str
  caller_list_str=$(printf '%s,' "${CALLERS[@]}")
  caller_list_str="${caller_list_str%,}"

  awk -v callers="$caller_list_str" '
    BEGIN {
      FS = OFS = "\t"
      n = split(callers, caller_arr, ",")
      for (i = 1; i <= n; i++) caller_idx[caller_arr[i]] = i
    }
    NF >= 6 {
      # Initialize all callers to 0
      for (i = 1; i <= n; i++) present[i] = 0

      # Parse callers list
      split($5, called, ",")
      for (i in called) {
        c = called[i]
        if (c in caller_idx) present[caller_idx[c]] = 1
      }

      # Output row
      printf "%s\t%s\t%s\t%s\t%s\t%s", $1, $2, $3, $4, $6, $5
      for (i = 1; i <= n; i++) printf "\t%d", present[i]
      printf "\n"
    }
  ' "$anno_tsv_default" >> "$detail_tsv"

  log "Compare mode complete."
}

# ---------------- mode dispatcher ----------------
case "$MODE" in
  union)
    parse_vcfs
    prepare_per_caller_vcfs
    build_union_and_anno
    log "Done (mode=union)."
    exit 0
    ;;
  kpass)
    # derive from union.anno.tsv only
    [[ -s "$anno_tsv" ]] || { log "ERROR: union.anno.tsv missing/empty: $anno_tsv"; exit 2; }
    kpass_from_union_anno "$K" "$anno_tsv" "$PASS_OUT" > /dev/null
    log "Done (mode=kpass). Wrote: $PASS_OUT"
    exit 0
    ;;
  emit)
    parse_vcfs
    prepare_per_caller_vcfs
    [[ -s "$PASS_IN" ]] || { log "ERROR: pass TSV missing/empty: $PASS_IN"; exit 2; }
    if [[ "$EMIT_VCF" == "true" ]]; then
      emit_consensus_vcf_from_pass "$PASS_IN" "${OUT_PREFIX}.consensus.vcf"
    else
      log "emit-vcf=false; skipping emit."
    fi
    log "Done (mode=emit)."
    exit 0
    ;;
  compare)
    run_compare_mode
    log "Done (mode=compare)."
    exit 0
    ;;
  all)
    parse_vcfs
    prepare_per_caller_vcfs
    build_union_and_anno

    # Create default pass file for this K
    kpass_from_union_anno "$K" "$anno_tsv_default" "$pass_tsv_default" > /dev/null
    pass_count=$(awk 'BEGIN{FS="\t"} $7==1{n++} END{print n+0}' "$pass_tsv_default")
    log "consensus.pass.tsv lines (pass=1): $pass_count"

    if [[ "$EMIT_VCF" == "true" ]]; then
      emit_consensus_vcf_from_pass "$pass_tsv_default" "${OUT_PREFIX}.consensus.vcf"
    else
      log "emit-vcf=false; skipping emit."
    fi

    log "Done (mode=all)."
    ;;
  *)
    log "ERROR: Unknown --mode '$MODE' (use union|kpass|emit|compare|all)"
    exit 2
    ;;
esac

if [[ "$KEEP_ALL" != "true" ]]; then
  log "keep-all-annot=false requested, but workdir is persistent by design (nothing deleted)."
fi

log "Summary:"
if [[ -f "$anno_tsv_default" ]]; then
  log "  Total unique variants (union.anno.tsv): $(wc -l < "$anno_tsv_default")"
fi
if [[ -f "$pass_tsv_default" ]]; then
  log "  Passing K=$K (consensus.pass.tsv): $(awk 'BEGIN{FS="\t"} $7==1{n++} END{print n+0}' "$pass_tsv_default")"
fi
