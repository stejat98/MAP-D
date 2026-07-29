#!/bin/bash
# Resume REGENIE GWAS pipeline
# Checks completed jobs and submits remaining phenotypes

set -e

SCRIPT_DIR="/path/to/project/regenie_pipeline/scripts"
INPUT_BASE="/path/to/project/regenie_pipeline/inputs"
OUTPUT_ROOT="/path/to/scratch/regenie_test"
PHENO_LIST_FILE="${SCRIPT_DIR}/pheno_list.txt"
PHENO_LIST_REMAINING="${SCRIPT_DIR}/pheno_list_remaining.txt"
SLURM_LOG_DIR="${SCRIPT_DIR}/../slurm"

echo "=========================================="
echo "REGENIE GWAS Resume Script"
echo "=========================================="

# Create slurm logs directory if it doesn't exist
mkdir -p "$SLURM_LOG_DIR"

# Generate full phenotype list if it doesn't exist
if [ ! -f "$PHENO_LIST_FILE" ]; then
    echo "Generating full phenotype list..."
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
    
    TOTAL_PHENOS=$(wc -l < "$PHENO_LIST_FILE")
    echo "Generated phenotype list with $TOTAL_PHENOS phenotypes"
else
    TOTAL_PHENOS=$(wc -l < "$PHENO_LIST_FILE")
    echo "Using existing phenotype list with $TOTAL_PHENOS phenotypes"
fi

# Check which phenotypes are completed
echo ""
echo "Checking completed jobs..."
echo "Output directory: $OUTPUT_ROOT"

rm -f "$PHENO_LIST_REMAINING"
COMPLETED=0
INCOMPLETE=0

while IFS='|' read -r STRATUM SAMPLE_TYPE PHENO_NAME; do
    # Check if step2 output exists (indicates completion)
    STEP2_OUTPUT="${OUTPUT_ROOT}/${PHENO_NAME}/${STRATUM}/step2/regenie_step2_${PHENO_NAME}.regenie.gz"
    
    if [ -f "$STEP2_OUTPUT" ]; then
        # Verify file is not empty
        if [ -s "$STEP2_OUTPUT" ]; then
            COMPLETED=$((COMPLETED + 1))
            echo "  ✓ Completed: ${STRATUM}|${SAMPLE_TYPE}|${PHENO_NAME}"
        else
            echo "  ✗ Incomplete (empty output): ${STRATUM}|${SAMPLE_TYPE}|${PHENO_NAME}"
            echo "${STRATUM}|${SAMPLE_TYPE}|${PHENO_NAME}" >> "$PHENO_LIST_REMAINING"
            INCOMPLETE=$((INCOMPLETE + 1))
        fi
    else
        # Check if step1 exists but step2 doesn't (partial completion)
        STEP1_PRED="${OUTPUT_ROOT}/${PHENO_NAME}/${STRATUM}/step1/regenie_step1_${PHENO_NAME}_pred.list"
        if [ -f "$STEP1_PRED" ]; then
            echo "  ⚠ Partial (step1 done, step2 missing): ${STRATUM}|${SAMPLE_TYPE}|${PHENO_NAME}"
        else
            echo "  ✗ Not started: ${STRATUM}|${SAMPLE_TYPE}|${PHENO_NAME}"
        fi
        echo "${STRATUM}|${SAMPLE_TYPE}|${PHENO_NAME}" >> "$PHENO_LIST_REMAINING"
        INCOMPLETE=$((INCOMPLETE + 1))
    fi
done < "$PHENO_LIST_FILE"

echo ""
echo "Summary:"
echo "  Completed: $COMPLETED"
echo "  Incomplete/Remaining: $INCOMPLETE"
echo "  Total: $TOTAL_PHENOS"

if [ ! -f "$PHENO_LIST_REMAINING" ] || [ ! -s "$PHENO_LIST_REMAINING" ]; then
    echo ""
    echo "All phenotypes are complete! Nothing to resume."
    exit 0
fi

REMAINING_COUNT=$(wc -l < "$PHENO_LIST_REMAINING")
echo ""
echo "Remaining phenotypes to process: $REMAINING_COUNT"
echo ""
echo "First 10 remaining entries:"
head -10 "$PHENO_LIST_REMAINING"
echo ""

# Ask for confirmation (unless --yes flag is provided)
if [ "$1" != "--yes" ]; then
    read -p "Submit $REMAINING_COUNT array jobs? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Cancelled."
        exit 0
    fi
fi

# Create temporary phenotype list for array job
PHENO_LIST_TEMP="${SCRIPT_DIR}/pheno_list_remaining_${SLURM_JOB_ID:-$$}.txt"
cp "$PHENO_LIST_REMAINING" "$PHENO_LIST_TEMP"

# Submit array job with remaining phenotypes
echo "Submitting SLURM array job for remaining phenotypes..."
echo "Using phenotype list: $PHENO_LIST_TEMP"

JOB_ID=$(sbatch \
    --array=1-${REMAINING_COUNT} \
    --export=PHENO_LIST_FILE="${PHENO_LIST_TEMP}" \
    "${SCRIPT_DIR}/run_regenie_array.sh" | grep -oP '\d+')

echo ""
echo "Job submitted with ID: $JOB_ID"
echo "Monitor with: squeue -j $JOB_ID"
echo "Check logs in: $SLURM_LOG_DIR"
echo ""
echo "To check progress, run this script again:"
echo "  bash ${SCRIPT_DIR}/resume_regenie.sh"
