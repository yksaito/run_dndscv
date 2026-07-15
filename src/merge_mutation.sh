#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -e log/merge_mutation/
#$ -o log/merge_mutation/
#$ -l s_vmem=32G

set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 <R_script> <manifest_list> <project_dir>" >&2
  exit 1
fi

R_SCRIPT="$1"
MANIFEST_LIST="$2"
PROJECT_DIR="$3"
TASK_ID="${SGE_TASK_ID:-1}"

cd "$PROJECT_DIR"

MANIFEST_FILE=$(sed -n "${TASK_ID}p" "$MANIFEST_LIST")

if [[ -z "$MANIFEST_FILE" ]]; then
  echo "Error: No manifest was found for SGE_TASK_ID=${TASK_ID}" >&2
  exit 1
fi

if [[ ! -f "$MANIFEST_FILE" ]]; then
  echo "Error: Manifest file not found: $MANIFEST_FILE" >&2
  exit 1
fi

if [[ ! -f "$R_SCRIPT" ]]; then
  echo "Error: R script not found: $R_SCRIPT" >&2
  exit 1
fi

echo "SGE_TASK_ID: $TASK_ID"
echo "Host: $(hostname)"
echo "R script: $R_SCRIPT"
echo "Manifest: $MANIFEST_FILE"
echo "Started: $(date '+%Y-%m-%d %H:%M:%S')"

Rscript "$R_SCRIPT" "$MANIFEST_FILE"

echo "Finished: $(date '+%Y-%m-%d %H:%M:%S')"
