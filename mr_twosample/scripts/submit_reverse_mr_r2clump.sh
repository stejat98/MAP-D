#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-00:45
#SBATCH --mem=16G
#SBATCH -c 2
#SBATCH -J rev_r2clump
#SBATCH -o /n/groups/patel/sivateja/regenie_pipeline/slurm/rev_r2clump_%A_%a.out
#SBATCH -e /n/groups/patel/sivateja/regenie_pipeline/slurm/rev_r2clump_%A_%a.err
#SBATCH --array=1-32
#
# Sensitivity reverse-MR with r2<0.01 LD clumping.
# Array layout: tasks 1-16 = BMI chunks 1-16 ; tasks 17-32 = TRIG_HDL_RATIO chunks 1-16
#
set -euo pipefail

PIPE=/n/groups/patel/sivateja/regenie_pipeline

module load gcc/14.2.0 R/4.4.2 plink/1.90b7.7_20241022
export PLINK_BIN=$(which plink)
# OpenGWAS JWT (only needed if clumping falls back to the remote API; local plink clump does not use it)
if [ -z "${OPENGWAS_JWT:-}" ] && [ -f "$PIPE/.secrets/opengwas_jwt" ]; then
  OPENGWAS_JWT=$(tr -d '\r\n' < "$PIPE/.secrets/opengwas_jwt")
  export OPENGWAS_JWT
fi

N_CHUNKS=16
TASK=${SLURM_ARRAY_TASK_ID}

if [ "$TASK" -le "$N_CHUNKS" ]; then
  PHENO="BMI"
  CHUNK=$TASK
else
  PHENO="TRIG_HDL_RATIO"
  CHUNK=$((TASK - N_CHUNKS))
fi

echo "Task $TASK -> PHENO=$PHENO CHUNK=$CHUNK / $N_CHUNKS  (plink=$PLINK_BIN)"

Rscript /n/groups/patel/sivateja/regenie_pipeline/scripts/reverse_mr_r2clump_sensitivity.R \
  "$PHENO" "$CHUNK" "$N_CHUNKS"
