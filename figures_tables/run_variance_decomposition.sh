#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -t 0-06:00
#SBATCH -p medium
#SBATCH -o slurm/variance_decomp_%j.out
#SBATCH -e slurm/variance_decomp_%j.err
#SBATCH -J var_decomp

module load gcc/14.2.0 R/4.5.2

# Optional overrides (see variance_decomposition_methodology_notes.txt):
# export VARDECOMP_MAX_N=30000
# export VARDECOMP_MAX_PROT_MISS=0.30
# export VARDECOMP_NFOLDS=3

Rscript /path/to/project/regenie_pipeline/variance_decomposition_bmi.R
