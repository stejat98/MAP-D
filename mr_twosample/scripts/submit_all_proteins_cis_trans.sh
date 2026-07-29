#!/bin/bash
#SBATCH --job-name=MR_cis_trans
#SBATCH --partition=short
#SBATCH --time=11:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=/n/groups/patel/sivateja/regenie_pipeline/slurm/MR_cis_trans_%j.out
#SBATCH --error=/n/groups/patel/sivateja/regenie_pipeline/slurm/MR_cis_trans_%j.err

# Load modules
module load gcc/14.2.0
module load R/4.4.2

# Set library path
export R_LIBS_USER="/n/groups/patel/sivateja/.R/library"

cd /n/groups/patel/sivateja/regenie_pipeline

echo "Starting ALL proteins cis vs trans MR analysis..."
echo "Date: $(date)"
echo "Hostname: $(hostname)"
echo ""

Rscript scripts/twosampleMR_all_proteins_cis_trans.R

echo ""
echo "Completed!"
echo "Date: $(date)"
