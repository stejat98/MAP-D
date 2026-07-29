#!/bin/bash
#SBATCH --job-name=LEP_cis_trans_MR
#SBATCH --partition=short
#SBATCH --time=1:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=2
#SBATCH --output=/path/to/project/regenie_pipeline/slurm/LEP_cis_trans_%j.out
#SBATCH --error=/path/to/project/regenie_pipeline/slurm/LEP_cis_trans_%j.err

# Load R module
module load gcc/14.2.0
module load R/4.4.2

# Set library path
export R_LIBS_USER="/path/to/project/.R/library"

cd /path/to/project/regenie_pipeline

echo "Starting LEP cis vs trans MR analysis..."
date

Rscript scripts/twosampleMR_leptin_cis_trans.R

echo "Completed!"
date
