#!/bin/bash
#SBATCH --job-name=DECODE_LEP_BMI_MR
#SBATCH --partition=short
#SBATCH --time=1:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --output=/n/groups/patel/sivateja/regenie_pipeline/slurm/decode_lep_bmi_%j.out
#SBATCH --error=/n/groups/patel/sivateja/regenie_pipeline/slurm/decode_lep_bmi_%j.err

# Load R module
module load gcc/14.2.0
module load R/4.4.2

# Set library path
export R_LIBS_USER="/n/groups/patel/sivateja/.R/library"

cd /n/groups/patel/sivateja/regenie_pipeline

echo "Starting DECODE LEP cis-pQTL -> BMI MR analysis..."
date

Rscript scripts/mr_decode_lep_bmi.R

echo "Completed!"
date
