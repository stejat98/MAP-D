#!/bin/bash
# Submit REGENIE jobs for BMI and HbA1c hallmarks for diabetes and prediabetes only
# Updated with longer time limits to prevent timeout (need ~50+ hours for all 22 chromosomes)

# Hallmarks (subset for testing)
HALLMARKS=("BMI" "HbA1c")

# Strata (submit in order: diabetes, prediabetes - NOT full)
STRATA=("diabetes" "prediabetes")

SCRIPT_DIR="/n/groups/patel/sivateja/regenie_pipeline/scripts"
REGENIE_SCRIPT="${SCRIPT_DIR}/run_regenie_array.sh"

echo "=========================================="
echo "Submitting REGENIE jobs for BMI and HbA1c"
echo "=========================================="
echo "Hallmarks: ${HALLMARKS[@]}"
echo "Strata: ${STRATA[@]}"
echo ""

# Create pheno list file for array job submission
PHENO_LIST="${SCRIPT_DIR}/pheno_list_hallmarks_bmi_hba1c.txt"
> "$PHENO_LIST"  # Clear/create file

JOB_COUNT=0
for STRATUM in "${STRATA[@]}"; do
    for PHENO in "${HALLMARKS[@]}"; do
        # Check if input files exist
        INPUT_PATH="/n/groups/patel/sivateja/regenie_pipeline/inputs/${STRATUM}/hallmarks_heldout/${PHENO}"
        if [ -f "${INPUT_PATH}/pheno.txt" ] && [ -f "${INPUT_PATH}/covar.txt" ]; then
            echo "${STRATUM}|hallmarks_heldout|${PHENO}" >> "$PHENO_LIST"
            ((JOB_COUNT++))
            echo "  Added: ${STRATUM} - ${PHENO}"
        else
            echo "  WARNING: Files not found for ${STRATUM} - ${PHENO}, skipping..."
        fi
    done
done

echo ""
echo "Total jobs to submit: $JOB_COUNT"
echo "Pheno list file: $PHENO_LIST"
echo ""

# Time estimates based on test job 24864971 (diabetes, BMI)
# Test: Step 1 ~57min, Step 2 ~11min (chr22 only) → Production ~3-4h (all chr)
# Scaling: Prediabetes ~5.3x, Full ~49.5x

echo "=========================================="
echo "Submitting jobs (optimized for fairshare)..."
echo "=========================================="

# Submit diabetes first (shortest, highest priority partition)
DIAB_COUNT=$(grep -c "^diabetes" "$PHENO_LIST" 2>/dev/null || echo "0")
if [ $DIAB_COUNT -gt 0 ]; then
    DIAB_START=$(grep -n "^diabetes" "$PHENO_LIST" | head -1 | cut -d: -f1)
    DIAB_END=$(grep -n "^diabetes" "$PHENO_LIST" | tail -1 | cut -d: -f1)
    
    echo "Submitting diabetes jobs (${DIAB_COUNT} jobs)..."
    echo "  Partition: medium (needed for longer time limit)"
    echo "  Time: 3-00:00:00 (3 days - estimated ~50+ hours needed for all 22 chromosomes)"
    echo "  Note: Previous jobs timed out with 12 hour limit (only completed 2-3/22 chromosomes)"
    echo "  Strategy: Using medium partition for longer time limit"
    
    sbatch \
      --array=${DIAB_START}-${DIAB_END} \
      --time=3-00:00:00 \
      --cpus-per-task=8 \
      --mem=50G \
      --partition=medium \
      --job-name=regenie_hallmarks_diab \
      --output=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_hallmarks_diab_%A_%a.out \
      --error=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_hallmarks_diab_%A_%a.err \
      --export=PHENO_LIST_FILE="$PHENO_LIST" \
      "$REGENIE_SCRIPT"
    
    if [ $? -eq 0 ]; then
        echo "  ✓ Diabetes jobs submitted to short partition"
    else
        echo "  ✗ Failed to submit diabetes jobs"
    fi
else
    echo "  No diabetes jobs to submit"
fi

# Submit prediabetes next (medium priority)
PREDIAB_COUNT=$(grep -c "^prediabetes" "$PHENO_LIST" 2>/dev/null || echo "0")
if [ $PREDIAB_COUNT -gt 0 ]; then
    PREDIAB_START=$(grep -n "^prediabetes" "$PHENO_LIST" | head -1 | cut -d: -f1)
    PREDIAB_END=$(grep -n "^prediabetes" "$PHENO_LIST" | tail -1 | cut -d: -f1)
    
    echo ""
    echo "Submitting prediabetes jobs (${PREDIAB_COUNT} jobs)..."
    echo "  Partition: medium (PriorityJobFactor=6)"
    echo "  Sample size: ~45K individuals (~5.3x diabetes)"
    echo "  Time: 5-00:00:00 (5 days - estimated ~110+ hours needed for all 22 chromosomes)"
    echo "  Note: Previous jobs timed out with 24 hour limit (only completed 2/22 chromosomes)"
    
    sbatch \
      --array=${PREDIAB_START}-${PREDIAB_END} \
      --time=5-00:00:00 \
      --cpus-per-task=8 \
      --mem=50G \
      --partition=medium \
      --job-name=regenie_hallmarks_prediab \
      --output=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_hallmarks_prediab_%A_%a.out \
      --error=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_hallmarks_prediab_%A_%a.err \
      --export=PHENO_LIST_FILE="$PHENO_LIST" \
      "$REGENIE_SCRIPT"
    
    if [ $? -eq 0 ]; then
        echo "  ✓ Prediabetes jobs submitted"
    else
        echo "  ✗ Failed to submit prediabetes jobs"
    fi
else
    echo ""
    echo "  No prediabetes jobs to submit"
fi

# Skip full stratum - only submitting diabetes and prediabetes
echo ""
echo "Skipping full stratum jobs (only submitting diabetes and prediabetes)"

echo ""
echo "=========================================="
echo "Job submission complete!"
echo "=========================================="
echo ""
echo "Summary:"
echo "  ✓ Diabetes: ${DIAB_COUNT} jobs (medium partition, ~50+ hours, 3-00:00:00 limit)"
echo "  ✓ Prediabetes: ${PREDIAB_COUNT} jobs (medium partition, ~110+ hours, 5-00:00:00 limit)"
echo "  ⚠ Full stratum: SKIPPED (only submitting BMI and HbA1c for diabetes and prediabetes)"
echo ""
echo "Monitor jobs:"
echo "  squeue -u st320"
echo "  squeue -u st320 -p medium   # Diabetes and prediabetes jobs"
echo ""
echo "Logs: /n/groups/patel/sivateja/regenie_pipeline/slurm/"
echo ""
