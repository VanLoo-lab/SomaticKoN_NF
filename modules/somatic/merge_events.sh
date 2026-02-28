#!/usr/bin/env bash
set -euo pipefail

# merge_events.sh
#
# Merge MNV and MULTIALLELIC variants back together based on EVENT_ID labels.
# Mark as FULL if all parts present, PARTIAL if only some parts present.
# Keep non-event variants as-is.

usage() {
  cat <<'EOF'
merge_events.sh

Merge split MNV/MULTIALLELIC variants back to their original form based on EVENT_ID.
Mark merged events as FULL (all parts present) or PARTIAL (some parts missing).

Required:
  --in VCF           Input VCF with EVENT_ID, EVENT_TYPE, EVENT_PART, EVENT_SIZE, EVENT_STATUS
  --out VCF          Output VCF with merged events
  --workdir PATH     Working directory for intermediate files

Optional:
  --threads INT      Default: 4
  --keep-split       If true, keep both split and merged variants (default: false)
  --log PATH         Log file (default: <workdir>/merge_events.log)

Output VCF INFO fields:
  - Merged events: MERGED_EVENT=1, MERGED_STATUS=FULL|PARTIAL, MERGED_FROM=<parts>
  - Non-event variants: unchanged
  - If --keep-split: split variants get SPLIT_OF_EVENT=<event_id>

Notes:
  - MNV merging: Reconstructs original REF>ALT at original position
  - MULTIALLELIC merging: Combines ALTs back to comma-separated list
  - FULL status: All EVENT_SIZE parts are present in input
  - PARTIAL status: Some parts are missing (gaps in EVENT_PART sequence)

EOF
}

IN_VCF=""
OUT_VCF=""
WORKDIR=""
THREADS=4
KEEP_SPLIT="false"
LOG=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --in) IN_VCF="$2"; shift 2;;
    --out) OUT_VCF="$2"; shift 2;;
    --workdir) WORKDIR="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --keep-split) KEEP_SPLIT="$2"; shift 2;;
    --log) LOG="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "ERROR: Unknown arg: $1" >&2; usage; exit 2;;
  esac
done

[[ -n "$IN_VCF" && -n "$OUT_VCF" && -n "$WORKDIR" ]] || {
  echo "ERROR: Missing required arguments" >&2
  usage
  exit 2
}

command -v bcftools >/dev/null 2>&1 || { echo "ERROR: bcftools not found" >&2; exit 127; }
command -v awk >/dev/null 2>&1 || { echo "ERROR: awk not found" >&2; exit 127; }

[[ -f "$IN_VCF" ]] || { echo "ERROR: Input VCF not found: $IN_VCF" >&2; exit 2; }

mkdir -p "$WORKDIR"

if [[ -z "$LOG" ]]; then
  LOG="$WORKDIR/merge_events.log"
fi

ts() { date +"%F %T"; }
log() { echo "[$(ts)] $*" | tee -a "$LOG" >&2; }

log "Starting merge_events"
log "Input: $IN_VCF"
log "Output: $OUT_VCF"
log "Workdir: $WORKDIR"
log "Keep split variants: $KEEP_SPLIT"

# Ensure input is uncompressed for processing
input_vcf="$WORKDIR/input.vcf"
if [[ "$IN_VCF" == *.gz ]]; then
  log "Decompressing input VCF"
  bcftools view "$IN_VCF" -Ov > "$input_vcf"
else
  cp "$IN_VCF" "$input_vcf"
fi

# Extract header
header="$WORKDIR/header.vcf"
grep '^#' "$input_vcf" > "$header"

# Add new INFO fields to header
new_header="$WORKDIR/header.new.vcf"
awk '
  /^#CHROM/ {
    print "##INFO=<ID=MERGED_EVENT,Number=0,Type=Flag,Description=\"This variant is a merged MNV or MULTIALLELIC event\">"
    print "##INFO=<ID=MERGED_STATUS,Number=1,Type=String,Description=\"FULL if all parts present, PARTIAL if some missing\">"
    print "##INFO=<ID=MERGED_FROM,Number=.,Type=String,Description=\"Comma-separated list of CHROM:POS:REF:ALT parts that were merged\">"
    if ("'"$KEEP_SPLIT"'" == "true") {
      print "##INFO=<ID=SPLIT_OF_EVENT,Number=1,Type=String,Description=\"EVENT_ID this split variant belongs to\">"
    }
  }
  {print}
' "$header" > "$new_header"

# Extract variants (body lines only)
variants="$WORKDIR/variants.tsv"
grep -v '^#' "$input_vcf" > "$variants"

log "Extracted $(wc -l < "$variants") variants"

# Process variants: identify events, merge them, output all
merged_body="$WORKDIR/merged.body.vcf"

awk -v keep_split="$KEEP_SPLIT" '
BEGIN {
  FS = OFS = "\t"
}

{
  chrom = $1
  pos = $2
  id = $3
  ref = $4
  alt = $5
  qual = $6
  filter = $7
  info = $8
  format_and_samples = ""
  for (i = 9; i <= NF; i++) {
    format_and_samples = format_and_samples (i > 9 ? OFS : "") $i
  }
  
  # Check if this variant has EVENT_ID
  event_id = ""
  event_type = ""
  event_part = ""
  event_size = ""
  event_status = ""
  
  if (match(info, /EVENT_ID=([^;]+)/, arr)) event_id = arr[1]
  if (match(info, /EVENT_TYPE=([^;]+)/, arr)) event_type = arr[1]
  if (match(info, /EVENT_PART=([^;]+)/, arr)) event_part = arr[1]
  if (match(info, /EVENT_SIZE=([^;]+)/, arr)) event_size = arr[1]
  if (match(info, /EVENT_STATUS=([^;]+)/, arr)) event_status = arr[1]
  
  if (event_id == "") {
    # Non-event variant: output as-is
    print chrom, pos, id, ref, alt, qual, filter, info, format_and_samples
    next
  }
  
  # Event variant: store for merging
  key = event_id
  
  # Store part info
  part_key = key SUBSEP event_part
  parts[key][event_part]["chrom"] = chrom
  parts[key][event_part]["pos"] = pos
  parts[key][event_part]["id"] = id
  parts[key][event_part]["ref"] = ref
  parts[key][event_part]["alt"] = alt
  parts[key][event_part]["qual"] = qual
  parts[key][event_part]["filter"] = filter
  parts[key][event_part]["info"] = info
  parts[key][event_part]["format_samples"] = format_and_samples
  
  # Track event metadata
  events[key]["type"] = event_type
  events[key]["size"] = event_size
  events[key]["status"] = event_status
  
  # Track which parts we have
  if (!(key in part_count)) part_count[key] = 0
  part_count[key]++
  
  if (!(key in min_part) || event_part < min_part[key]) min_part[key] = event_part
  if (!(key in max_part) || event_part > max_part[key]) max_part[key] = event_part
}

END {
  # Process each event
  for (event_id in events) {
    n_parts = part_count[event_id]
    expected_size = events[event_id]["size"]
    event_type = events[event_id]["type"]
    
    # Determine merged status: FULL if all parts present
    merged_status = (n_parts == expected_size) ? "FULL" : "PARTIAL"
    
    if (event_type == "MNV") {
      merge_mnv(event_id, n_parts, expected_size, merged_status)
    } else if (event_type == "MULTIALLELIC") {
      merge_multiallelic(event_id, n_parts, expected_size, merged_status)
    }
  }
}

function merge_mnv(event_id, n_parts, expected_size, merged_status,    
                   i, min_pos, ref_merged, alt_merged, chrom, first_part,
                   qual, filter, info, format_samples, merged_from, merged_info) {
  
  # Find the first part to get base position
  first_part = min_part[event_id]
  chrom = parts[event_id][first_part]["chrom"]
  min_pos = parts[event_id][first_part]["pos"]
  
  # Reconstruct MNV by concatenating bases in order
  ref_merged = ""
  alt_merged = ""
  merged_from = ""
  
  for (i = 1; i <= expected_size; i++) {
    if (i in parts[event_id]) {
      ref_merged = ref_merged parts[event_id][i]["ref"]
      alt_merged = alt_merged parts[event_id][i]["alt"]
      
      part_key = parts[event_id][i]["chrom"] ":" parts[event_id][i]["pos"] ":" \
                 parts[event_id][i]["ref"] ":" parts[event_id][i]["alt"]
      merged_from = merged_from (merged_from == "" ? "" : ",") part_key
    } else {
      # Missing part - use N for both
      ref_merged = ref_merged "N"
      alt_merged = alt_merged "N"
    }
  }
  
  # Use first part for QUAL, FILTER, FORMAT
  qual = parts[event_id][first_part]["qual"]
  filter = parts[event_id][first_part]["filter"]
  format_samples = parts[event_id][first_part]["format_samples"]
  
  # Merge INFO from all available parts
  merged_info = merge_info_fields(event_id, n_parts, expected_size)
  
  # Add MERGED_* tags
  merged_tags = "MERGED_EVENT;MERGED_STATUS=" merged_status ";MERGED_FROM=" merged_from
  
  if (merged_info == "" || merged_info == ".")
    info = merged_tags
  else
    info = merged_info ";" merged_tags
  
  # Output merged variant
  print chrom, min_pos, ".", ref_merged, alt_merged, qual, filter, info, format_samples
  
  # Optionally output split variants with original EVENT tags intact
  if (keep_split == "true") {
    for (i = 1; i <= expected_size; i++) {
      if (i in parts[event_id]) {
        split_info = parts[event_id][i]["info"] ";SPLIT_OF_EVENT=" event_id
        print parts[event_id][i]["chrom"], parts[event_id][i]["pos"], \
              parts[event_id][i]["id"], parts[event_id][i]["ref"], \
              parts[event_id][i]["alt"], parts[event_id][i]["qual"], \
              parts[event_id][i]["filter"], split_info, \
              parts[event_id][i]["format_samples"]
      }
    }
  }
}

function merge_info_fields(event_id, n_parts, expected_size,    
                            i, info_hash, key, val, tag, result, first_part) {
  # Collect all INFO tags from all parts
  # For conflicting tags, use first part value
  # For unique tags, keep all
  
  first_part = min_part[event_id]
  
  for (i = 1; i <= expected_size; i++) {
    if (!(i in parts[event_id])) continue
    
    info_str = parts[event_id][i]["info"]
    if (info_str == "." || info_str == "") continue
    
    # Parse INFO field
    n_tags = split(info_str, tags, ";")
    for (j = 1; j <= n_tags; j++) {
      tag = tags[j]
      if (tag == "") continue
      
      # Check if it is key=value or flag
      if (index(tag, "=") > 0) {
        split(tag, kv, "=")
        key = kv[1]
        val = kv[2]
        
        # Special handling for EVENT_PART - keep only from first part
        if (key == "EVENT_PART") {
          if (i == first_part) info_hash[key] = val
          continue
        }
        
        # For other tags, prefer first part if conflict
        if (!(key in info_hash)) {
          info_hash[key] = val
        }
      } else {
        # Flag (no value)
        info_hash[tag] = ""
      }
    }
  }
  
  # Reconstruct INFO string
  result = ""
  for (key in info_hash) {
    if (info_hash[key] == "") {
      # Flag
      result = result (result == "" ? "" : ";") key
    } else {
      # Key=value
      result = result (result == "" ? "" : ";") key "=" info_hash[key]
    }
  }
  
  return result
}

function merge_multiallelic(event_id, n_parts, expected_size, merged_status,
                             i, chrom, pos, ref, alt_list, qual, filter, info,
                             format_samples, merged_from, first_part, merged_info) {
  
  # For multiallelic: all parts should have same CHROM, POS, REF
  first_part = min_part[event_id]
  chrom = parts[event_id][first_part]["chrom"]
  pos = parts[event_id][first_part]["pos"]
  ref = parts[event_id][first_part]["ref"]
  
  # Collect all ALT alleles in order
  alt_list = ""
  merged_from = ""
  
  for (i = 1; i <= expected_size; i++) {
    if (i in parts[event_id]) {
      alt_list = alt_list (alt_list == "" ? "" : ",") parts[event_id][i]["alt"]
      
      part_key = parts[event_id][i]["chrom"] ":" parts[event_id][i]["pos"] ":" \
                 parts[event_id][i]["ref"] ":" parts[event_id][i]["alt"]
      merged_from = merged_from (merged_from == "" ? "" : ",") part_key
    } else {
      # Missing allele - use . placeholder
      alt_list = alt_list (alt_list == "" ? "" : ",") "."
    }
  }
  
  # Use first part for QUAL, FILTER, FORMAT
  qual = parts[event_id][first_part]["qual"]
  filter = parts[event_id][first_part]["filter"]
  format_samples = parts[event_id][first_part]["format_samples"]
  
  # Merge INFO from all available parts
  merged_info = merge_info_fields(event_id, n_parts, expected_size)
  
  # Add MERGED_* tags
  merged_tags = "MERGED_EVENT;MERGED_STATUS=" merged_status ";MERGED_FROM=" merged_from
  
  if (merged_info == "" || merged_info == ".")
    info = merged_tags
  else
    info = merged_info ";" merged_tags
  
  # Output merged variant
  print chrom, pos, ".", ref, alt_list, qual, filter, info, format_samples
  
  # Optionally output split variants with original EVENT tags intact
  if (keep_split == "true") {
    for (i = 1; i <= expected_size; i++) {
      if (i in parts[event_id]) {
        split_info = parts[event_id][i]["info"] ";SPLIT_OF_EVENT=" event_id
        print parts[event_id][i]["chrom"], parts[event_id][i]["pos"], \
              parts[event_id][i]["id"], parts[event_id][i]["ref"], \
              parts[event_id][i]["alt"], parts[event_id][i]["qual"], \
              parts[event_id][i]["filter"], split_info, \
              parts[event_id][i]["format_samples"]
      }
    }
  }
}
' "$variants" > "$merged_body"

log "Processed events, generated $(wc -l < "$merged_body") output variants"

# Combine header and body
cat "$new_header" "$merged_body" > "$OUT_VCF.tmp"

# Sort the output
log "Sorting output VCF"
grep '^#' "$OUT_VCF.tmp" > "$OUT_VCF"
grep -v '^#' "$OUT_VCF.tmp" | sort -k1,1V -k2,2n >> "$OUT_VCF"

rm -f "$OUT_VCF.tmp"

log "Wrote: $OUT_VCF"

# Summary statistics
total_variants=$(wc -l < "$merged_body")
merged_count=$(grep -c 'MERGED_EVENT' "$merged_body" || echo "0")
full_count=$(grep -c 'MERGED_STATUS=FULL' "$merged_body" || echo "0")
partial_count=$(grep -c 'MERGED_STATUS=PARTIAL' "$merged_body" || echo "0")

log "Summary:"
log "  Total output variants: $total_variants"
log "  Merged events: $merged_count"
log "  - FULL (all parts present): $full_count"
log "  - PARTIAL (some parts missing): $partial_count"
log "  Non-event variants: $((total_variants - merged_count))"

log "Done."