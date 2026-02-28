#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
standardize_vcf_samples.sh

Auto-detect tumor/normal sample IDs from a 2-sample VCF and standardize:
  - reorder to Tumor,Normal (or user-chosen names)
  - rename samples to standardized names (default: Tumor / Normal)

Required:
  --in-vcf   PATH   (vcf or vcf.gz)
  --out-vcf  PATH   (must end with .vcf.gz)

Optional (recommended):
  --tumor-hint   STRING   (substring match, case-insensitive; e.g. "MK545-A" or "tumor")
  --normal-hint  STRING   (substring match, case-insensitive; e.g. "Control" or "normal")

Optional (more flexible):
  --tumor-regex  REGEX    (regex match, case-sensitive by awk; you can embed (?i) in some regex engines but awk is basic)
  --normal-regex REGEX

Optional naming:
  --tumor-name   STRING   default: Tumor
  --normal-name  STRING   default: Normal

Other:
  --threads  INT  default: 2

Heuristics if no hints are provided:
  - tries to infer tumor/normal from common tokens in sample names:
    tumor:  tumor|case|primary|met|lesion|TUM
    normal: normal|ctrl|control|blood|germline|NORM

If it cannot confidently decide, it exits with an error and prints the sample IDs.
EOF
}

IN=""
OUT=""
THREADS=2
TUMOR_HINT=""
NORMAL_HINT=""
TUMOR_REGEX=""
NORMAL_REGEX=""
TUMOR_NAME="Tumor"
NORMAL_NAME="Normal"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --in-vcf) IN="$2"; shift 2;;
    --out-vcf) OUT="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --tumor-hint) TUMOR_HINT="$2"; shift 2;;
    --normal-hint) NORMAL_HINT="$2"; shift 2;;
    --tumor-regex) TUMOR_REGEX="$2"; shift 2;;
    --normal-regex) NORMAL_REGEX="$2"; shift 2;;
    --tumor-name) TUMOR_NAME="$2"; shift 2;;
    --normal-name) NORMAL_NAME="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1" >&2; usage; exit 2;;
  esac
done

[[ -n "$IN" && -n "$OUT" ]] || { echo "ERROR: --in-vcf and --out-vcf are required" >&2; usage; exit 2; }
[[ -f "$IN" ]] || { echo "ERROR: input not found: $IN" >&2; exit 2; }
[[ "$OUT" == *.vcf.gz ]] || { echo "ERROR: --out-vcf must end with .vcf.gz" >&2; exit 2; }

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

# Ensure bgzip+index for input (bcftools query -l works on plain VCF too, but reheader wants bgzipped reliably)
IN_GZ="$IN"
if [[ "$IN" != *.vcf.gz ]]; then
  IN_GZ="${tmpdir}/in.vcf.gz"
  bcftools view "$IN" -Oz --threads "$THREADS" -o "$IN_GZ"
  tabix -f -p vcf "$IN_GZ"
else
  # create index if missing
  if [[ ! -f "${IN}.tbi" && ! -f "${IN}.csi" ]]; then
    tabix -f -p vcf "$IN"
  fi
fi

mapfile -t SAMPLES < <(bcftools query -l "$IN_GZ")
if [[ "${#SAMPLES[@]}" -ne 2 ]]; then
  echo "ERROR: expected exactly 2 samples in VCF, found ${#SAMPLES[@]}:" >&2
  printf '  - %s\n' "${SAMPLES[@]}" >&2
  exit 2
fi

s1="${SAMPLES[0]}"
s2="${SAMPLES[1]}"

# Helper: case-insensitive substring match using bash
ci_contains() {
  local hay="$1" needle="$2"
  shopt -s nocasematch
  [[ "$hay" == *"$needle"* ]]
  local rc=$?
  shopt -u nocasematch
  return $rc
}

pick_by_hint_or_regex() {
  local a="$1" b="$2"
  local tumor="" normal=""

  # 1) regex mode (if provided)
  if [[ -n "$TUMOR_REGEX" || -n "$NORMAL_REGEX" ]]; then
    if [[ -n "$TUMOR_REGEX" ]]; then
      if echo "$a" | awk -v r="$TUMOR_REGEX" '$0 ~ r {exit 0} {exit 1}'; then tumor="$a"; fi
      if echo "$b" | awk -v r="$TUMOR_REGEX" '$0 ~ r {exit 0} {exit 1}'; then tumor="${tumor:-$b}"; fi
    fi
    if [[ -n "$NORMAL_REGEX" ]]; then
      if echo "$a" | awk -v r="$NORMAL_REGEX" '$0 ~ r {exit 0} {exit 1}'; then normal="$a"; fi
      if echo "$b" | awk -v r="$NORMAL_REGEX" '$0 ~ r {exit 0} {exit 1}'; then normal="${normal:-$b}"; fi
    fi
  fi

  # 2) hint mode (if provided)
  if [[ -z "$tumor" && -n "$TUMOR_HINT" ]]; then
    if ci_contains "$a" "$TUMOR_HINT"; then tumor="$a"; fi
    if ci_contains "$b" "$TUMOR_HINT"; then tumor="${tumor:-$b}"; fi
  fi
  if [[ -z "$normal" && -n "$NORMAL_HINT" ]]; then
    if ci_contains "$a" "$NORMAL_HINT"; then normal="$a"; fi
    if ci_contains "$b" "$NORMAL_HINT"; then normal="${normal:-$b}"; fi
  fi

  echo "$tumor|$normal"
}

# Heuristic fallback if no hints match
pick_by_heuristic() {
  local a="$1" b="$2"
  local tumor="" normal=""

  # tokens
  local tumor_tokens=("tumor" "case" "primary" "met" "lesion" "tum" "t")
  local normal_tokens=("normal" "ctrl" "control" "blood" "germline" "norm" "n")

  # score function
  score_name() {
    local name="$1"
    local score_t=0 score_n=0
    for tok in "${tumor_tokens[@]}"; do
      if ci_contains "$name" "$tok"; then score_t=$((score_t+1)); fi
    done
    for tok in "${normal_tokens[@]}"; do
      if ci_contains "$name" "$tok"; then score_n=$((score_n+1)); fi
    done
    echo "$score_t,$score_n"
  }

  IFS=',' read -r a_t a_n <<<"$(score_name "$a")"
  IFS=',' read -r b_t b_n <<<"$(score_name "$b")"

  # Decide only if clearly separable
  if [[ "$a_t" -gt "$a_n" && "$b_n" -gt "$b_t" ]]; then
    tumor="$a"; normal="$b"
  elif [[ "$b_t" -gt "$b_n" && "$a_n" -gt "$a_t" ]]; then
    tumor="$b"; normal="$a"
  fi

  echo "$tumor|$normal"
}

tumor=""
normal=""

# First try hints/regex if user provided anything
if [[ -n "$TUMOR_HINT" || -n "$NORMAL_HINT" || -n "$TUMOR_REGEX" || -n "$NORMAL_REGEX" ]]; then
  IFS='|' read -r tumor normal <<<"$(pick_by_hint_or_regex "$s1" "$s2")"
fi

# If still not decided, try heuristic
if [[ -z "$tumor" || -z "$normal" ]]; then
  IFS='|' read -r ht hn <<<"$(pick_by_heuristic "$s1" "$s2")"
  tumor="${tumor:-$ht}"
  normal="${normal:-$hn}"
fi

# If still ambiguous, fail loudly with guidance
if [[ -z "$tumor" || -z "$normal" || "$tumor" == "$normal" ]]; then
  echo "ERROR: Could not confidently determine tumor/normal samples from VCF header." >&2
  echo "Found sample IDs:" >&2
  printf '  - %s\n' "${SAMPLES[@]}" >&2
  echo "" >&2
  echo "Fix by providing hints or regex, e.g.:" >&2
  echo "  --tumor-hint 'MK545-A' --normal-hint 'Control'" >&2
  echo "or:" >&2
  echo "  --tumor-regex 'A$' --normal-regex 'Control$'" >&2
  exit 2
fi

echo ${IN_GZ}
# Reorder (tumor first, then normal)
bcftools view -s "${tumor},${normal}" -Oz --threads "$THREADS" \
  -o "${tmpdir}/reordered.vcf.gz" "$IN_GZ"
tabix -f -p vcf "${tmpdir}/reordered.vcf.gz"

# Rename to standardized names
cat > "${tmpdir}/rename.txt" <<EOF
${tumor}	${TUMOR_NAME}
${normal}	${NORMAL_NAME}
EOF

bcftools reheader -s "${tmpdir}/rename.txt" \
  -o "$OUT" "${tmpdir}/reordered.vcf.gz"

bcftools index -f -t "$OUT"

echo "Detected:"
echo "  tumor : $tumor  -> $TUMOR_NAME"
echo "  normal: $normal -> $NORMAL_NAME"
echo "Wrote:"
echo "  $OUT"
echo "  ${OUT}.tbi"
