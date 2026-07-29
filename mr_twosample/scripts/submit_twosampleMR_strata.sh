#!/bin/bash
#SBATCH --job-name=twosampleMR_strata
#SBATCH --partition=short
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=/n/groups/patel/sivateja/regenie_pipeline/slurm/twosampleMR_strata_%A_%a.out
#SBATCH --error=/n/groups/patel/sivateja/regenie_pipeline/slurm/twosampleMR_strata_%A_%a.err
#SBATCH --array=1-4

# TwoSampleMR Analysis: Validated Proteins -> BMI/HbA1c (Prediabetes & Diabetes Strata)
# This script submits array jobs for all 4 MR analyses

# Load R module
module load gcc/9.2.0
module load R/4.2.1

# Set working directory
cd /n/groups/patel/sivateja/regenie_pipeline

# Define scripts for each array task
SCRIPTS=(
    "scripts/twosampleMR_proteins_bmi_prediabetes.R"
    "scripts/twosampleMR_proteins_bmi_diabetes.R"
    "scripts/twosampleMR_proteins_hba1c_prediabetes.R"
    "scripts/twosampleMR_proteins_hba1c_diabetes.R"
)

# Get the script for this array task (0-indexed)
SCRIPT_IDX=$((SLURM_ARRAY_TASK_ID - 1))
SCRIPT=${SCRIPTS[$SCRIPT_IDX]}

echo "=============================================="
echo "TwoSampleMR Analysis - Array Task ${SLURM_ARRAY_TASK_ID}"
echo "Script: ${SCRIPT}"
echo "Started: $(date)"
echo "=============================================="

# Run the R script
Rscript ${SCRIPT}

echo "=============================================="
echo "Completed: $(date)"
echo "=============================================="
