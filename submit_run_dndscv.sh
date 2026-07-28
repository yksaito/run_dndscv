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
  MAX_CONCURRENT   Maximum simultaneously running tasks in each array (default: 2)
  REFDB            dndscv reference database passed to run_dndscv.R (default: hg38)

The normal and MSIHexcluded analyses are submitted as two independent array jobs.
Completed outputs are excluded from each array input list.
USAGE
}

if [[ $# -gt 0 && ( "$1" == "-h" || "$1" == "--help" ) ]]; then
  usage
  exit 0
fi

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_DIR=$(pwd)
R_SCRIPT="${SCRIPT_DIR}/src/R/run_dndscv.R"
ARRAY_JOB_SCRIPT="${SCRIPT_DIR}/src/run_dndscv.sh"
INPUT_DIR="${PROJECT_DIR}/output/merged_mutation"
OUTPUT_DIR="${PROJECT_DIR}/output/dndscv"
LIST_DIR="${PROJECT_DIR}/output/dndscv_lists"
LOG_DIR="${PROJECT_DIR}/log/dndscv"
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

mkdir -p "$LIST_DIR" "$LOG_DIR"

TIMESTAMP=$(date '+%Y%m%d_%H%M%S')
CATEGORY_LABEL=$(IFS=-; echo "${CATEGORIES[*]}")

declare -A INPUT_LISTS
declare -A N_CANDIDATES
declare -A N_SKIPPED
declare -A N_TASKS

prepare_analysis_list() {
  local analysis_type="$1"
  local list_label
  local input_list
  local candidate_list

  if [[ "$analysis_type" == "normal" ]]; then
    list_label="normal"
  else
    list_label="MSIHexcluded"
  fi

  input_list="${LIST_DIR}/dndscv_${CATEGORY_LABEL}_${list_label}_${TIMESTAMP}.txt"
  candidate_list="${input_list}.candidates"
  : > "$candidate_list"

  for category in "${CATEGORIES[@]}"; do
    if [[ "$analysis_type" == "normal" ]]; then
      find "$INPUT_DIR" -maxdepth 1 -type f \
        -name "${category}_*_merged_mutation_for_dndscv.txt" \
        ! -name "*_MSIHexcluded_merged_mutation_for_dndscv.txt" -print
    else
      find "$INPUT_DIR" -maxdepth 1 -type f \
        -name "${category}_*_MSIHexcluded_merged_mutation_for_dndscv.txt" -print
    fi
  done | sort -u > "$candidate_list"

  local n_candidates
  n_candidates=$(wc -l < "$candidate_list" | tr -d ' ')
  if [[ "$n_candidates" -eq 0 ]]; then
    echo "Error: No ${list_label} merged mutation files matched the requested categories." >&2
    rm -f "$candidate_list" "$input_list"
    return 1
  fi

  : > "$input_list"
  local n_skipped=0
  local input_file
  local input_basename
  local analysis_id
  local sel_cv_output
  local all_output

  while IFS= read -r input_file; do
    input_basename=$(basename "$input_file")
    analysis_id=${input_basename%_merged_mutation_for_dndscv.txt}
    sel_cv_output="${OUTPUT_DIR}/sel_cv_${analysis_id}.csv"
    all_output="${OUTPUT_DIR}/dndsCV_${analysis_id}.xlsx"

    if [[ -s "$sel_cv_output" && -s "$all_output" ]]; then
      echo "Skip completed ${list_label} cohort: $analysis_id"
      n_skipped=$((n_skipped + 1))
    elif [[ -e "$sel_cv_output" || -e "$all_output" ]]; then
      echo "Error: Incomplete dNdScv output for ${list_label} cohort: $analysis_id" >&2
      echo "Expected two non-empty files:" >&2
      echo "  $sel_cv_output" >&2
      echo "  $all_output" >&2
      rm -f "$candidate_list" "$input_list"
      return 1
    else
      printf '%s\n' "$input_file" >> "$input_list"
    fi
  done < "$candidate_list"

  rm -f "$candidate_list"

  local n_tasks
  n_tasks=$(wc -l < "$input_list" | tr -d ' ')
  if [[ "$n_tasks" -eq 0 ]]; then
    rm -f "$input_list"
    input_list=""
  fi

  INPUT_LISTS[$analysis_type]="$input_list"
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
  echo "Analysis type           : $list_label"
  echo "Categories              : ${CATEGORIES[*]}"
  echo "Reference database      : $REFDB"
  echo "Matched cohorts         : ${N_CANDIDATES[$analysis_type]}"
  echo "Skipped completed       : ${N_SKIPPED[$analysis_type]}"
  echo "Number of array tasks   : ${N_TASKS[$analysis_type]}"
  echo "Maximum concurrent tasks: $MAX_CONCURRENT"

  if [[ "${N_TASKS[$analysis_type]}" -eq 0 ]]; then
    echo "All ${list_label} cohorts already have complete dNdScv outputs; qsub was not called."
    return 0
  fi

  echo "Input list              : ${INPUT_LISTS[$analysis_type]}"
  qsub \
    -N "run_dnds_${job_suffix}" \
    -t "1-${N_TASKS[$analysis_type]}" \
    -tc "$MAX_CONCURRENT" \
    -o "$LOG_DIR" \
    -e "$LOG_DIR" \
    "$ARRAY_JOB_SCRIPT" \
    "$R_SCRIPT" \
    "${INPUT_LISTS[$analysis_type]}" \
    "$PROJECT_DIR" \
    "$REFDB"
}

submit_analysis normal
submit_analysis msihexcluded
