#!/bin/bash

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash submit_merge_mutation.sh <merge|dndscv> [category ...]

Examples:
  bash submit_merge_mutation.sh merge
  bash submit_merge_mutation.sh merge CODE
  bash submit_merge_mutation.sh merge CODE level1 level2
  MAX_CONCURRENT=3 bash submit_merge_mutation.sh dndscv CODE

Categories:
  CODE level1 level2 level3 level4 level5

Environment variable:
  MAX_CONCURRENT   Maximum simultaneously running array tasks (default: 5)
USAGE
}

if [[ $# -lt 1 ]]; then
  usage >&2
  exit 1
fi

MODE="$1"
shift

case "$MODE" in
  merge)
    R_SCRIPT_NAME="merge_mutation.R"
    ;;
  dndscv)
    R_SCRIPT_NAME="merge_mutation_for_dndscv.R"
    ;;
  *)
    echo "Error: mode must be 'merge' or 'dndscv'." >&2
    usage >&2
    exit 1
    ;;
esac

# submit_merge_mutation.shはプロジェクトルートに置く
# Rスクリプト: src/R/
# array worker: src/merge_mutation.sh
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_DIR=$(pwd)
R_SCRIPT="${SCRIPT_DIR}/src/R/${R_SCRIPT_NAME}"
ARRAY_JOB_SCRIPT="${SCRIPT_DIR}/src/merge_mutation.sh"
MANIFEST_DIR="${PROJECT_DIR}/output/manifest"
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
    exit 1
  fi
done

LIST_DIR="${PROJECT_DIR}/output/manifest_lists"
mkdir -p "$LIST_DIR"

TIMESTAMP=$(date '+%Y%m%d_%H%M%S')
CATEGORY_LABEL=$(IFS=-; echo "${CATEGORIES[*]}")
MANIFEST_LIST="${LIST_DIR}/${MODE}_${CATEGORY_LABEL}_${TIMESTAMP}.txt"

# manifestを絶対パスで列挙し、同じファイルが複数回入らないようにする
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

echo "Mode: $MODE"
echo "Project directory: $PROJECT_DIR"
echo "Categories: ${CATEGORIES[*]}"
echo "Manifest list: $MANIFEST_LIST"
echo "Number of tasks: $N_TASKS"
echo "Maximum concurrent tasks: $MAX_CONCURRENT"

qsub \
  -t "1-${N_TASKS}" \
  -tc "$MAX_CONCURRENT" \
  "$ARRAY_JOB_SCRIPT" \
  "$R_SCRIPT" \
  "$MANIFEST_LIST" \
  "$PROJECT_DIR"
