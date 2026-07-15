#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -e log/dndscv/
#$ -o log/dndscv/
#$ -l s_vmem=32G
#$ -pe def_slot 2

set -euo pipefail

module use /usr/local/package/modulefiles/
module load R/4.4.0

if [[ $# -ne 4 ]]; then
  echo "Usage: $0 <R_script> <input_list> <project_dir> <refdb>" >&2
  exit 1
fi

R_SCRIPT="$1"
INPUT_LIST="$2"
PROJECT_DIR="$3"
REFDB="$4"
TASK_ID="${SGE_TASK_ID:-1}"

# 1 task内で意図せず多数のCPUスレッドを使用しないようにする
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

cd "$PROJECT_DIR"

INPUT_FILE=$(sed -n "${TASK_ID}p" "$INPUT_LIST")

if [[ -z "$INPUT_FILE" ]]; then
  echo "Error: No input file was found for SGE_TASK_ID=${TASK_ID}" >&2
  exit 1
fi

if [[ ! -f "$INPUT_FILE" ]]; then
  echo "Error: Input file not found: $INPUT_FILE" >&2
  exit 1
fi

if [[ ! -f "$R_SCRIPT" ]]; then
  echo "Error: R script not found: $R_SCRIPT" >&2
  exit 1
fi

echo "========================================"
echo "SGE_TASK_ID: $TASK_ID"
echo "Host       : $(hostname)"
echo "R script   : $R_SCRIPT"
echo "Input file : $INPUT_FILE"
echo "Reference  : $REFDB"
echo "Started    : $(date '+%Y-%m-%d %H:%M:%S')"
echo "========================================"

Rscript "$R_SCRIPT" "$INPUT_FILE" "$REFDB"

echo "Finished   : $(date '+%Y-%m-%d %H:%M:%S')"
