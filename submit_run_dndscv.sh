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
  MAX_CONCURRENT=10 bash submit_run_dndscv.sh CODE
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

# submit_run_dndscv.shはプロジェクトルートに置く
# Rスクリプト: src/R/run_dndscv.R
# array worker: src/run_dndscv.sh
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_DIR=$(pwd)
R_SCRIPT="${SCRIPT_DIR}/src/R/run_dndscv.R"
ARRAY_JOB_SCRIPT="${SCRIPT_DIR}/src/run_dndscv.sh"
INPUT_DIR="${PROJECT_DIR}/output/merged_mutation"
MAX_CONCURRENT="${MAX_CONCURRENT:-10}"
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

# dNdScv用に統合済みのファイルを絶対パスで列挙する
{
  for category in "${CATEGORIES[@]}"; do
    find "$INPUT_DIR" -maxdepth 1 -type f \
      -name "${category}_*_merged_mutation_for_dndscv.txt" -print
  done
} | sort -u > "$INPUT_LIST"

N_TASKS=$(wc -l < "$INPUT_LIST" | tr -d ' ')

if [[ "$N_TASKS" -eq 0 ]]; then
  echo "Error: No merged mutation files matched the requested categories." >&2
  echo "Searched directory: $INPUT_DIR" >&2
  rm -f "$INPUT_LIST"
  exit 1
fi

echo "Project directory       : $PROJECT_DIR"
echo "Categories              : ${CATEGORIES[*]}"
echo "Reference database      : $REFDB"
echo "Input list              : $INPUT_LIST"
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
