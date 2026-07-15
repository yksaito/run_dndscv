#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -e log/manifest/
#$ -o log/manifest/
#$ -l s_vmem=8G

module use /usr/local/package/modulefiles/
module load R/4.4.0

Rscript src/R/make_manifest.R
