#!/bin/bash
#SBATCH -c 8
#SBATCH -t 2-00:00:00
#SBATCH --mem=50G
#SBATCH -p medium
#SBATCH -o /n/groups/patel/sivateja/regenie_pipeline/slurm/test_prediabetes_%j.out
#SBATCH -e /n/groups/patel/sivateja/regenie_pipeline/slurm/test_prediabetes_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sivateja_tangirala@hms.harvard.edu
#SBATCH --job-name=regenie_test_pred

# Test REGENIE for prediabetes stratum
/n/groups/patel/sivateja/regenie_pipeline/scripts/run_regenie_array.sh prediabetes proteins_only protein_x.1004
