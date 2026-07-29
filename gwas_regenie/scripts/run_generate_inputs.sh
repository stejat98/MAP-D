#!/bin/bash
#SBATCH -c 4
#SBATCH -t 4:00:00
#SBATCH --mem=100G
#SBATCH -p short
#SBATCH -o /path/to/project/regenie_pipeline/slurm/generate_inputs_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@example.com

# Run generate_inputs.R on a compute node with more memory

module load gcc/14.2.0 R/4.4.2

cd /path/to/project

Rscript regenie_pipeline/scripts/generate_inputs.R 2>&1 | tee regenie_pipeline/generate_inputs.log

echo "Script completed. Check log file: regenie_pipeline/generate_inputs.log"
