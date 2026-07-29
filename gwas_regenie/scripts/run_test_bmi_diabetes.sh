#!/bin/bash
#SBATCH -c 8
#SBATCH -t 11:59:00
#SBATCH --mem=50G
#SBATCH -p short
#SBATCH -o /n/groups/patel/sivateja/regenie_pipeline/slurm/test_bmi_diabetes_%j.out
#SBATCH -e /n/groups/patel/sivateja/regenie_pipeline/slurm/test_bmi_diabetes_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sivateja_tangirala@hms.harvard.edu

# Conda environment setup (IMPORTANT: Set before any conda commands)
export HOME=/n/groups/patel/sivateja
export CONDA_NO_PLUGINS=true
export CONDA_PKGS_DIRS=/n/groups/patel/sivateja/.conda/pkgs
export CONDA_ENVS_PATH=/n/groups/patel/sivateja/.conda/envs

cd /n/groups/patel/sivateja
bash regenie_pipeline/scripts/test_bmi_diabetes.sh
