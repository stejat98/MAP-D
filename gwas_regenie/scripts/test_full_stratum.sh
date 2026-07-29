#!/bin/bash
#SBATCH -c 8
#SBATCH -t 2-00:00:00
#SBATCH --mem=50G
#SBATCH -p medium
#SBATCH -o /path/to/project/regenie_pipeline/slurm/test_full_%j.out
#SBATCH -e /path/to/project/regenie_pipeline/slurm/test_full_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@example.com
#SBATCH --job-name=regenie_test_full

# Test REGENIE for full stratum
/path/to/project/regenie_pipeline/scripts/run_regenie_array.sh full proteins_only protein_x.1004
