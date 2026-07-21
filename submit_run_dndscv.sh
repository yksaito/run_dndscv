#!/bin/bash

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash submit_run_dndscv.sh [category ...]

Examples:
  bash submit_run_dndscv.sh
  bash submit_run_dndscv.sh CODE
  bash submit_run_dndscv.sh CODE level1 level2
  MAX_CONCURRENT=2 bash submit_run_dndscv.sh CODE
  REFDB=RefCDS_human_GRCh38.p12_dNdScv.0.1.0.rda bash submit_run_dndscv.sh CODE

Categories:
  CODE level1 level2 level3 level4 level5

Environment variables:
  MAX_CONCURRENT   Maximum simultaneously running array tasks (default: 2)
  REFDB            dndscv reference database passed to run_dndscv.R (default: hg38)
USAGE
}

if [[ $# -gt 0 && ( "$1" == "-h" || "$1" == "--help" ) ]]; then
  usage
  exit 0
fi

# submit_run_dndscv.sh縺ｯ繝励Ο繧ｸ繧ｧ繧ｯ繝医Ν繝ｼ繝医↓鄂ｮ縺�
# R繧ｹ繧ｯ繝ｪ繝励ヨ: src/R/run_dndscv.R
# array worker: src/run_dndscv.sh
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_DIR=$(pwd)
R_SCRIPT="${SCRIPT_DIR}/src/R/run_dndscv.R"
ARRAY_JOB_SCRIPT="${SCRIPT_DIR}/src/run_dndscv.sh"
INPUT_DIR="${PROJECT_DIR}/output/merged_mutation"
OUTPUT_DIR="${PROJECT_DIR}/output/dndscv"
MAX_CONCURRENT="${MAX_CONCURRENT:-2}"
REFDB="${REFDB:-hg38}"

if [[ ! -f "$R_SCRIPT" ]]; then
  echo "Error: R script not found: $R_SCRIPT" >&2
  exit 1
fi

if [[ ! -f "$ARRAY_JOB_SCRIPT" ]]; then
  echo "Error: Array job script not found: $ARRAY_JOB_SCRIPT" >&2
  exit 1
fi

if [[ ! -d "$INPUT_DIR" ]]; then
  echo "Error: Input directory not found: $INPUT_DIR" >&2
  exit 1
fi

if ! [[ "$MAX_CONCURRENT" =~ ^[1-9][0-9]*$ ]]; then
  echo "Error: MAX_CONCURRENT must be a positive integer." >&2
  exit 1
fi

if [[ $# -gt 0 ]]; then
  CATEGORIES=("$@")
else
  CATEGORIES=(CODE level1 level2 level3 level4 level5)
fi

VALID_CATEGORIES=" CODE level1 level2 level3 level4 level5 "
for category in "${CATEGORIES[@]}"; do
  if [[ "$VALID_CATEGORIES" != *" ${category} "* ]]; then
    echo "Error: Invalid category: $category" >&2
    usage >&2
    exit 1
  fi
done

LIST_DIR="${PROJECT_DIR}/output/dndscv_lists"
LOG_DIR="${PROJECT_DIR}/log/dndscv"
mkdir -p "$LIST_DIR" "$LOG_DIR"

TIMESTAMP=$(date '+%Y%m%d_%H%M%S')
CATEGORY_LABEL=$(IFS=-; echo "${CATEGORIES[*]}")
INPUT_LIST="${LIST_DIR}/dndscv_${CATEGORY_LABEL}_${TIMESTAMP}.txt"
CANDIDATE_LIST="${INPUT_LIST}.candidates"

trap 'rm -f -- "$CANDIDATE_LIST"' EXIT

# dNdScv逕ｨ縺ｫ邨ｱ蜷域ｸ医∩縺ｮ繝輔ぃ繧､繝ｫ繧堤ｵｶ蟇ｾ繝代せ縺ｧ蛻玲嫌縺吶ｋ
{
  for category in "${CATEGORIES[@]}"; do
    find "$INPUT_DIR" -maxdepth 1 -type f \
      -name "${category}_*_merged_mutation_for_dndscv.txt" -print
  done
} | sort -u > "$CANDIDATE_LIST"

N_CANDIDATES=$(wc -l < "$CANDIDATE_LIST" | tr -d ' ')

if [[ "$N_CANDIDATES" -eq 0 ]]; then
  echo "Error: No merged mutation files matched the requested categories." >&2
  echo "Searched directory: $INPUT_DIR" >&2
  exit 1
fi

# 螳御ｺ�ｸ医∩cohort繧帝勁縺阪∵悴螳溯｡慶ohort縺�縺代ｒarray input list縺ｸ譖ｸ縺榊�縺吶�
# sel_cv縺ｨ蜈ｨ蜃ｺ蜉妣xcel縺ｮ迚�婿縺�縺代′縺ゅｋ蝣ｴ蜷医ｄ0 byte縺ｮ蝣ｴ蜷医�荳榊ｮ悟�蜃ｺ蜉帙→縺励※蛛懈ｭ｢縺吶ｋ縲�
: > "$INPUT_LIST"
N_SKIPPED=0

while IFS= read -r input_file; do
  input_basename=$(basename "$input_file")
  analysis_id=${input_basename%_merged_mutation_for_dndscv.txt}

  sel_cv_output="${OUTPUT_DIR}/sel_cv_${analysis_id}.csv"
  all_output="${OUTPUT_DIR}/dndsCV_${analysis_id}.xlsx"

  if [[ -s "$sel_cv_output" && -s "$all_output" ]]; then
    echo "Skip completed cohort: $analysis_id"
    N_SKIPPED=$((N_SKIPPED + 1))
  elif [[ -e "$sel_cv_output" || -e "$all_output" ]]; then
    echo "Error: Incomplete dNdScv output for cohort: $analysis_id" >&2
    echo "Expected two non-empty files:" >&2
    echo "  $sel_cv_output" >&2
    echo "  $all_output" >&2
    rm -f "$INPUT_LIST"
    exit 1
  else
    printf '%s\n' "$input_file" >> "$INPUT_LIST"
  fi
done < "$CANDIDATE_LIST"

rm -f "$CANDIDATE_LIST"
trap - EXIT

N_TASKS=$(wc -l < "$INPUT_LIST" | tr -d ' ')

if [[ "$N_TASKS" -eq 0 ]]; then
  echo "All matched cohorts already have complete dNdScv outputs."
  echo "Matched cohorts : $N_CANDIDATES"
  echo "Skipped cohorts : $N_SKIPPED"
  rm -f "$INPUT_LIST"
  exit 0
fi

echo "Project directory       : $PROJECT_DIR"
echo "Categories              : ${CATEGORIES[*]}"
echo "Reference database      : $REFDB"
echo "Input list              : $INPUT_LIST"
echo "Matched cohorts         : $N_CANDIDATES"
echo "Skipped completed       : $N_SKIPPED"
echo "Number of array tasks   : $N_TASKS"
echo "Maximum concurrent tasks: $MAX_CONCURRENT"

qsub \
  -N "run_dndscv" \
  -t "1-${N_TASKS}" \
  -tc "$MAX_CONCURRENT" \
  -o "$LOG_DIR" \
  -e "$LOG_DIR" \
  "$ARRAY_JOB_SCRIPT" \
  "$R_SCRIPT" \
  "$INPUT_LIST" \
  "$PROJECT_DIR" \
  "$REFDB"
  