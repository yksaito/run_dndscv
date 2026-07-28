#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l s_vmem=32G

set -euo pipefail

module use /usr/local/package/modulefiles/
module load R/4.4.0

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <R_script> <project_dir>" >&2
  exit 1
fi

R_SCRIPT="$1"
PROJECT_DIR="$2"

cd "$PROJECT_DIR"

if [[ ! -f "$R_SCRIPT" ]]; then
  echo "Error: R script not found: $R_SCRIPT" >&2
  exit 1
fi

echo "Host     : $(hostname)"
echo "R script : $R_SCRIPT"
echo "Started  : $(date '+%Y-%m-%d %H:%M:%S')"

Rscript "$R_SCRIPT"

echo "Finished : $(date '+%Y-%m-%d %H:%M:%S')"
