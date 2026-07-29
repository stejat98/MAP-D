#!/bin/bash
#SBATCH --job-name=reverse_PEWAS
#SBATCH --partition=short
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH --output=/path/to/project/regenie_pipeline/slurm/reverse_pewas_%j.out
#SBATCH --error=/path/to/project/regenie_pipeline/slurm/reverse_pewas_%j.err

module load gcc/14.2.0
module load R/4.4.2

export R_LIBS_USER="/path/to/home/R/x86_64-pc-linux-gnu-library/4.4"

cd /path/to/project/regenie_pipeline

echo "Starting Reverse PEWAS: Phenotype -> Protein regressions"
echo "Job ID: $SLURM_JOB_ID"
date

Rscript scripts/reverse_observational_pewas.R

echo "Finished"
date
