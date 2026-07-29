#!/bin/bash
# Submit REGENIE GWAS array jobs
# Generates phenotype list and submits SLURM array jobs

set -e

SCRIPT_DIR="/path/to/project/regenie_pipeline/scripts"
INPUT_BASE="/path/to/project/regenie_pipeline/inputs"
PHENO_LIST_FILE="${SCRIPT_DIR}/pheno_list.txt"

echo "=========================================="
echo "REGENIE GWAS Job Submission"
echo "=========================================="

# Generate phenotype list
echo "Generating phenotype list..."
rm -f "$PHENO_LIST_FILE"

STRATA=("full")
SAMPLE_TYPES=("proteins_only" "hallmarks_heldout")

for STRATUM in "${STRATA[@]}"; do
    for SAMPLE_TYPE in "${SAMPLE_TYPES[@]}"; do
        INPUT_DIR="${INPUT_BASE}/${STRATUM}/${SAMPLE_TYPE}"
        
        if [ ! -d "$INPUT_DIR" ]; then
            echo "WARNING: Input directory not found: $INPUT_DIR"
            continue
        fi
        
        # Find all phenotype directories
        for PHENO_DIR in "${INPUT_DIR}"/*; do
            if [ -d "$PHENO_DIR" ] && [ -f "${PHENO_DIR}/pheno.txt" ] && [ -f "${PHENO_DIR}/covar.txt" ]; then
                PHENO_NAME=$(basename "$PHENO_DIR")
                echo "${STRATUM}|${SAMPLE_TYPE}|${PHENO_NAME}" >> "$PHENO_LIST_FILE"
            fi
        done
    done
done

TOTAL_JOBS=$(wc -l < "$PHENO_LIST_FILE")
echo "Found $TOTAL_JOBS phenotypes to process"

if [ "$TOTAL_JOBS" -eq 0 ]; then
    echo "ERROR: No phenotypes found. Run generate_inputs.R first."
    exit 1
fi

# Show first few entries
echo ""
echo "First 10 entries:"
head -10 "$PHENO_LIST_FILE"
echo ""

# Ask for confirmation (unless --yes flag is provided)
if [ "$1" != "--yes" ]; then
    read -p "Submit $TOTAL_JOBS array jobs? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Cancelled."
        exit 0
    fi
fi

# Submit array job
echo "Submitting SLURM array job..."
sbatch \
    --array=1-${TOTAL_JOBS} \
    "${SCRIPT_DIR}/run_regenie_array.sh"

echo ""
echo "Job submitted!"
echo "Monitor with: squeue -u $USER"
echo "Check logs in: /path/to/project/regenie_pipeline/slurm/"
