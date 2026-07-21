#!/bin/bash

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash submit_merge_mutation_drivers.sh [category ...]

Examples:
  bash submit_merge_mutation_drivers.sh CODE
  bash submit_merge_mutation_drivers.sh CODE level1 level2
  MAX_CONCURRENT=5 bash submit_merge_mutation_drivers.sh CODE

Categories:
  CODE level1 level2 level3 level4 level5

Environment variable:
  MAX_CONCURRENT   Maximum simultaneously running array tasks (default: 20)
USAGE
}

if [[ $# -gt 0 && ( "$1" == "-h" || "$1" == "--help" ) ]]; then
  usage
  exit 0
fi

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_DIR=$(pwd)
R_SCRIPT="${SCRIPT_DIR}/src/R/merge_mutation_drivers.R"
ARRAY_JOB_SCRIPT="${SCRIPT_DIR}/src/merge_mutation.sh"
MANIFEST_DIR="${PROJECT_DIR}/output/manifest"
LIST_DIR="${PROJECT_DIR}/output/manifest_lists"
LOG_DIR="${PROJECT_DIR}/log/merge_drv"
MAX_CONCURRENT="${MAX_CONCURRENT:-20}"

if [[ ! -f "$R_SCRIPT" ]]; then
  echo "Error: R script not found: $R_SCRIPT" >&2
  exit 1
fi

if [[ ! -f "$ARRAY_JOB_SCRIPT" ]]; then
  echo "Error: Array job script not found: $ARRAY_JOB_SCRIPT" >&2
  exit 1
fi

if [[ ! -d "$MANIFEST_DIR" ]]; then
  echo "Error: Manifest directory not found: $MANIFEST_DIR" >&2
  exit 1
fi

if [[ ! -d "$LOG_DIR" ]]; then
  echo "Error: Log directory not found: $LOG_DIR" >&2
  echo "Create it before submission." >&2
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

mkdir -p "$LIST_DIR"

TIMESTAMP=$(date '+%Y%m%d_%H%M%S')
CATEGORY_LABEL=$(IFS=-; echo "${CATEGORIES[*]}")
MANIFEST_LIST="${LIST_DIR}/drivers_${CATEGORY_LABEL}_${TIMESTAMP}.txt"

{
  for category in "${CATEGORIES[@]}"; do
    find "$MANIFEST_DIR" -maxdepth 1 -type f \
      -name "${category}_*_manifest.txt" -print
  done
} | sort -u > "$MANIFEST_LIST"

N_TASKS=$(wc -l < "$MANIFEST_LIST" | tr -d ' ')

if [[ "$N_TASKS" -eq 0 ]]; then
  echo "Error: No manifest files matched the requested categories." >&2
  rm -f "$MANIFEST_LIST"
  exit 1
fi

echo "Project directory: $PROJECT_DIR"
echo "Categories: ${CATEGORIES[*]}"
echo "Manifest list: $MANIFEST_LIST"
echo "Number of tasks: $N_TASKS"
echo "Maximum concurrent tasks: $MAX_CONCURRENT"

qsub \
  -N "merge_drv" \
  -t "1-${N_TASKS}" \
  -tc "$MAX_CONCURRENT" \
  -o "$LOG_DIR" \
  -e "$LOG_DIR" \
  "$ARRAY_JOB_SCRIPT" \
  "$R_SCRIPT" \
  "$MANIFEST_LIST" \
  "$PROJECT_DIR"
