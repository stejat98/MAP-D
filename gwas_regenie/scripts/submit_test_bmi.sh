#!/bin/bash
#SBATCH -c 4
#SBATCH -t 2:00:00
#SBATCH --mem=30G
#SBATCH -p short
#SBATCH -o /path/to/project/regenie_pipeline/slurm/test_bmi_%j.out
#SBATCH -e /path/to/project/regenie_pipeline/slurm/test_bmi_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@example.com

# Quick test of REGENIE with BMI phenotype
# Tests chromosome 22 only with small block sizes

# Conda environment setup (IMPORTANT: Set before any conda commands)
export HOME=/path/to/project
export CONDA_NO_PLUGINS=true
export CONDA_PKGS_DIRS=/path/to/project/.conda/pkgs
export CONDA_ENVS_PATH=/path/to/project/.conda/envs

cd /path/to/project
bash regenie_pipeline/scripts/test_bmi_quick_subset.sh
