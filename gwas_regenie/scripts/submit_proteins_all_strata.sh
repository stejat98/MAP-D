#!/bin/bash
# Submit REGENIE jobs for all proteins across all strata
# Uses same covariate approach as hallmarks (should be faster)

SCRIPT_DIR="/n/groups/patel/sivateja/regenie_pipeline/scripts"
REGENIE_SCRIPT="${SCRIPT_DIR}/run_regenie_array.sh"

# Get list of all proteins from the inputs directory
PROTEINS_DIR="/n/groups/patel/sivateja/regenie_pipeline/inputs/full/proteins_only"
STRATA=("diabetes" "prediabetes" "full")

echo "=========================================="
echo "Submitting REGENIE jobs for proteins"
echo "=========================================="

# Get all protein names
PROTEINS=($(ls "$PROTEINS_DIR" | sort))

echo "Found ${#PROTEINS[@]} proteins"
echo "Strata: ${STRATA[@]}"
echo ""

# Create pheno list file
PHENO_LIST="${SCRIPT_DIR}/pheno_list_proteins.txt"
> "$PHENO_LIST"

JOB_COUNT=0
for STRATUM in "${STRATA[@]}"; do
    for PROTEIN in "${PROTEINS[@]}"; do
        # Check if input files exist
        INPUT_PATH="/n/groups/patel/sivateja/regenie_pipeline/inputs/${STRATUM}/proteins_only/${PROTEIN}"
        if [ -f "${INPUT_PATH}/pheno.txt" ] && [ -f "${INPUT_PATH}/covar.txt" ]; then
            echo "${STRATUM}|proteins_only|${PROTEIN}" >> "$PHENO_LIST"
            ((JOB_COUNT++))
        fi
    done
done

echo "Total jobs to submit: $JOB_COUNT"
echo "Pheno list file: $PHENO_LIST"
echo ""

# Estimate time: Proteins should be similar to hallmarks
# Step 1: ~8-10 hours, Step 2: ~22-44 hours (all chromosomes)
# Total: ~1.25-2.25 days per job

# Submit diabetes and prediabetes first (smaller, faster)
DIAB_PREDIAB_COUNT=$(grep -E "^diabetes|^prediabetes" "$PHENO_LIST" | wc -l)
if [ $DIAB_PREDIAB_COUNT -gt 0 ]; then
    echo "Submitting diabetes and prediabetes protein jobs (${DIAB_PREDIAB_COUNT} jobs)..."
    sbatch \
      --array=1-${DIAB_PREDIAB_COUNT} \
      --time=2-00:00:00 \
      --cpus-per-task=8 \
      --mem=50G \
      --partition=medium \
      --job-name=regenie_proteins \
      --output=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_proteins_%A_%a.out \
      --error=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_proteins_%A_%a.err \
      --export=PHENO_LIST_FILE="$PHENO_LIST" \
      "$REGENIE_SCRIPT"
    
    if [ $? -eq 0 ]; then
        echo "  ✓ Diabetes and prediabetes protein jobs submitted"
    else
        echo "  ✗ Failed to submit diabetes/prediabetes jobs"
    fi
fi

# Submit full stratum jobs
FULL_START=$(grep -n "^full" "$PHENO_LIST" | head -1 | cut -d: -f1)
FULL_END=$(wc -l < "$PHENO_LIST")
FULL_COUNT=$((FULL_END - FULL_START + 1))

if [ $FULL_COUNT -gt 0 ]; then
    echo ""
    echo "Submitting full stratum protein jobs (${FULL_COUNT} jobs)..."
    sbatch \
      --array=${FULL_START}-${FULL_END} \
      --time=3-00:00:00 \
      --cpus-per-task=8 \
      --mem=50G \
      --partition=medium \
      --job-name=regenie_proteins_full \
      --output=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_proteins_full_%A_%a.out \
      --error=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_proteins_full_%A_%a.err \
      --export=PHENO_LIST_FILE="$PHENO_LIST" \
      "$REGENIE_SCRIPT"
    
    if [ $? -eq 0 ]; then
        echo "  ✓ Full stratum protein jobs submitted"
    else
        echo "  ✗ Failed to submit full stratum jobs"
    fi
fi

echo ""
echo "=========================================="
echo "Job submission complete!"
echo "=========================================="
echo "Total jobs: $JOB_COUNT"
echo "  - Diabetes/Prediabetes: $DIAB_PREDIAB_COUNT"
echo "  - Full: $FULL_COUNT"
echo ""
echo "Check status with: squeue -u st320"
echo "Monitor logs in: /n/groups/patel/sivateja/regenie_pipeline/slurm/"
echo ""
