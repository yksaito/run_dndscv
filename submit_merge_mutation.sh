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
  MAX_CONCURRENT   Maximum simultaneously running tasks in each array (default: 20)

The normal and MSIHexcluded analyses are submitted as two independent array jobs.
Completed outputs are excluded from each array input list.
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

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_DIR=$(pwd)
R_SCRIPT="${SCRIPT_DIR}/src/R/${R_SCRIPT_NAME}"
ARRAY_JOB_SCRIPT="${SCRIPT_DIR}/src/merge_mutation.sh"
MANIFEST_DIR="${PROJECT_DIR}/output/manifest"
OUTPUT_DIR="${PROJECT_DIR}/output/merged_mutation"
LIST_DIR="${PROJECT_DIR}/output/manifest_lists"
LOG_DIR="${PROJECT_DIR}/log/merge_mutation"
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
    exit 1
  fi
done

mkdir -p "$LIST_DIR"

TIMESTAMP=$(date '+%Y%m%d_%H%M%S')
CATEGORY_LABEL=$(IFS=-; echo "${CATEGORIES[*]}")
TARGET_CHRS=(
  chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12
  chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY
)

declare -A INPUT_LISTS
declare -A N_CANDIDATES
declare -A N_SKIPPED
declare -A N_TASKS

check_output_status() {
  local manifest_id="$1"

  if [[ "$MODE" == "dndscv" ]]; then
    local output_file="${OUTPUT_DIR}/${manifest_id}_merged_mutation_for_dndscv.txt"
    if [[ -s "$output_file" ]]; then
      return 0
    fi
    if [[ -e "$output_file" ]]; then
      echo "Error: Empty output exists: $output_file" >&2
      return 2
    fi
    return 1
  fi

  local output_dir="${OUTPUT_DIR}/${manifest_id}_by_chr"
  if [[ ! -e "$output_dir" ]]; then
    return 1
  fi
  if [[ ! -d "$output_dir" ]]; then
    echo "Error: Expected output directory but found another file type: $output_dir" >&2
    return 2
  fi

  local chr_name
  local output_file
  for chr_name in "${TARGET_CHRS[@]}"; do
    output_file="${output_dir}/${manifest_id}_merged_mutation_${chr_name}.txt"
    if [[ ! -s "$output_file" ]]; then
      echo "Error: Incomplete chromosome-split output: $output_file" >&2
      return 2
    fi
  done
  return 0
}

prepare_analysis_list() {
  local analysis_type="$1"
  local list_label
  local manifest_list
  local candidate_list

  if [[ "$analysis_type" == "normal" ]]; then
    list_label="normal"
  else
    list_label="MSIHexcluded"
  fi

  manifest_list="${LIST_DIR}/${MODE}_${CATEGORY_LABEL}_${list_label}_${TIMESTAMP}.txt"
  candidate_list="${manifest_list}.candidates"
  : > "$candidate_list"

  for category in "${CATEGORIES[@]}"; do
    if [[ "$analysis_type" == "normal" ]]; then
      find "$MANIFEST_DIR" -maxdepth 1 -type f \
        -name "${category}_*_manifest.txt" \
        ! -name "*_MSIHexcluded_manifest.txt" -print
    else
      find "$MANIFEST_DIR" -maxdepth 1 -type f \
        -name "${category}_*_MSIHexcluded_manifest.txt" -print
    fi
  done | sort -u > "$candidate_list"

  local n_candidates
  n_candidates=$(wc -l < "$candidate_list" | tr -d ' ')
  if [[ "$n_candidates" -eq 0 ]]; then
    echo "Error: No ${list_label} manifest files matched the requested categories." >&2
    rm -f "$candidate_list" "$manifest_list"
    return 1
  fi

  : > "$manifest_list"
  local n_skipped=0
  local manifest_file
  local manifest_id
  local status

  while IFS= read -r manifest_file; do
    manifest_id=$(basename "$manifest_file")
    manifest_id=${manifest_id%_manifest.txt}

    if check_output_status "$manifest_id"; then
      echo "Skip completed ${list_label} cohort: $manifest_id"
      n_skipped=$((n_skipped + 1))
    else
      status=$?
      if [[ "$status" -eq 2 ]]; then
        rm -f "$candidate_list" "$manifest_list"
        return 1
      fi
      printf '%s\n' "$manifest_file" >> "$manifest_list"
    fi
  done < "$candidate_list"

  rm -f "$candidate_list"

  local n_tasks
  n_tasks=$(wc -l < "$manifest_list" | tr -d ' ')
  if [[ "$n_tasks" -eq 0 ]]; then
    rm -f "$manifest_list"
    manifest_list=""
  fi

  INPUT_LISTS[$analysis_type]="$manifest_list"
  N_CANDIDATES[$analysis_type]="$n_candidates"
  N_SKIPPED[$analysis_type]="$n_skipped"
  N_TASKS[$analysis_type]="$n_tasks"
}

# Validate both analysis sets before submitting either job.
prepare_analysis_list normal
prepare_analysis_list msihexcluded

submit_analysis() {
  local analysis_type="$1"
  local list_label
  local job_suffix

  if [[ "$analysis_type" == "normal" ]]; then
    list_label="normal"
    job_suffix="norm"
  else
    list_label="MSIHexcluded"
    job_suffix="msi"
  fi

  echo "========================================"
  echo "Mode                    : $MODE"
  echo "Analysis type           : $list_label"
  echo "Categories              : ${CATEGORIES[*]}"
  echo "Matched cohorts         : ${N_CANDIDATES[$analysis_type]}"
  echo "Skipped completed       : ${N_SKIPPED[$analysis_type]}"
  echo "Number of array tasks   : ${N_TASKS[$analysis_type]}"
  echo "Maximum concurrent tasks: $MAX_CONCURRENT"

  if [[ "${N_TASKS[$analysis_type]}" -eq 0 ]]; then
    echo "All ${list_label} cohorts already have complete outputs; qsub was not called."
    return 0
  fi

  local job_prefix
  if [[ "$MODE" == "merge" ]]; then
    job_prefix="merge_mut"
  else
    job_prefix="prep_dnds"
  fi

  echo "Input list              : ${INPUT_LISTS[$analysis_type]}"
  qsub \
    -N "${job_prefix}_${job_suffix}" \
    -t "1-${N_TASKS[$analysis_type]}" \
    -tc "$MAX_CONCURRENT" \
    -o "$LOG_DIR" \
    -e "$LOG_DIR" \
    "$ARRAY_JOB_SCRIPT" \
    "$R_SCRIPT" \
    "${INPUT_LISTS[$analysis_type]}" \
    "$PROJECT_DIR"
}

submit_analysis normal
submit_analysis msihexcluded
