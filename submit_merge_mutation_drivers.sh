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

For every selected cohort, the normal and _MSIHexcluded manifests are
submitted as separate array tasks. A complete three-file output bundle is
skipped independently for each task.

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
DND_DIR="${PROJECT_DIR}/output/dndscv"
OUT_DIR="${PROJECT_DIR}/output/merged_mutation_drivers"
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
ALL_MANIFESTS=$(mktemp "${LIST_DIR}/drivers_${CATEGORY_LABEL}_${TIMESTAMP}.all.XXXXXX")
TASK_LIST=$(mktemp "${LIST_DIR}/drivers_${CATEGORY_LABEL}_${TIMESTAMP}.XXXXXX.txt")

cleanup() {
  rm -f -- "$ALL_MANIFESTS"
}
trap cleanup EXIT

{
  for category in "${CATEGORIES[@]}"; do
    find "$MANIFEST_DIR" -maxdepth 1 -type f \
      -name "${category}_*_manifest.txt" -print
  done
} | sort -u > "$ALL_MANIFESTS"

if [[ ! -s "$ALL_MANIFESTS" ]]; then
  echo "Error: No manifest files matched the requested categories." >&2
  exit 1
fi

: > "$TASK_LIST"
TOTAL_MANIFESTS=0
SKIPPED_TASKS=0
ERRORS=0

while IFS= read -r manifest_file; do
  [[ -n "$manifest_file" ]] || continue
  TOTAL_MANIFESTS=$((TOTAL_MANIFESTS + 1))

  manifest_id=$(basename "$manifest_file")
  manifest_id=${manifest_id%_manifest.txt}

  if [[ "$manifest_id" == *_MSIHexcluded ]]; then
    cohort=${manifest_id%_MSIHexcluded}
    mutation_out="${OUT_DIR}/${cohort}_merged_mutation_drivers_MSIHexcluded.txt"
    gene_out="${OUT_DIR}/${cohort}_significant_genes_q_driver_lt_0_1_MSIHexcluded.tsv"
    missing_gene_out="${OUT_DIR}/${cohort}_significant_genes_not_in_ensGene_MSIHexcluded.tsv"
  else
    cohort="$manifest_id"
    mutation_out="${OUT_DIR}/${cohort}_merged_mutation_drivers.txt"
    gene_out="${OUT_DIR}/${cohort}_significant_genes_q_driver_lt_0_1.tsv"
    missing_gene_out="${OUT_DIR}/${cohort}_significant_genes_not_in_ensGene.tsv"
  fi

  dndscv_file="${DND_DIR}/sel_cv_${manifest_id}.csv"
  if [[ ! -s "$dndscv_file" ]]; then
    echo "Error: Missing or empty dNdScv result: $dndscv_file" >&2
    ERRORS=$((ERRORS + 1))
    continue
  fi

  complete_count=0
  existing_count=0
  for output_file in "$mutation_out" "$gene_out" "$missing_gene_out"; do
    if [[ -s "$output_file" ]]; then
      complete_count=$((complete_count + 1))
    fi
    if [[ -e "$output_file" ]]; then
      existing_count=$((existing_count + 1))
    fi
  done

  if [[ "$complete_count" -eq 3 ]]; then
    echo "Skip complete task: $manifest_id"
    SKIPPED_TASKS=$((SKIPPED_TASKS + 1))
  elif [[ "$existing_count" -eq 0 ]]; then
    printf '%s\n' "$manifest_file" >> "$TASK_LIST"
  else
    echo "Error: Incomplete output bundle for $manifest_id" >&2
    echo "Expected either no output files or three non-empty files:" >&2
    echo "  $mutation_out" >&2
    echo "  $gene_out" >&2
    echo "  $missing_gene_out" >&2
    ERRORS=$((ERRORS + 1))
  fi
done < "$ALL_MANIFESTS"

if [[ "$ERRORS" -gt 0 ]]; then
  rm -f -- "$TASK_LIST"
  echo "Error: Found $ERRORS incomplete prerequisite/output set(s); qsub was not called." >&2
  exit 1
fi

N_TASKS=$(wc -l < "$TASK_LIST" | tr -d ' ')

echo "Project directory: $PROJECT_DIR"
echo "Categories: ${CATEGORIES[*]}"
echo "Selected normal/MSI-H-excluded manifests: $TOTAL_MANIFESTS"
echo "Skipped complete tasks: $SKIPPED_TASKS"
echo "Array tasks to submit: $N_TASKS"
echo "Maximum concurrent tasks: $MAX_CONCURRENT"

if [[ "$N_TASKS" -eq 0 ]]; then
  rm -f -- "$TASK_LIST"
  echo "All selected driver-mutation output bundles are complete; qsub was not called."
  exit 0
fi

echo "Manifest list: $TASK_LIST"

qsub \
  -N "merge_drv" \
  -t "1-${N_TASKS}" \
  -tc "$MAX_CONCURRENT" \
  -o "$LOG_DIR" \
  -e "$LOG_DIR" \
  "$ARRAY_JOB_SCRIPT" \
  "$R_SCRIPT" \
  "$TASK_LIST" \
  "$PROJECT_DIR"
