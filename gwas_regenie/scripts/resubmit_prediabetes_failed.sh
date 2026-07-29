#!/bin/bash
# Resubmit failed prediabetes jobs
# - Tasks 15 and 37: Chromosome 15 failures (NaN error at block 49)
# - Task 24: Chromosome 2 timeout (needs more time)

SCRIPT_DIR="/n/groups/patel/sivateja/regenie_pipeline/scripts"
REGENIE_SCRIPT="${SCRIPT_DIR}/run_regenie_array_per_chr.sh"
PHENO_LIST="${SCRIPT_DIR}/pheno_list_prediabetes_per_chr.txt"

echo "=========================================="
echo "Resubmitting Failed Prediabetes Jobs"
echo "=========================================="
echo ""
echo "Failed jobs to resubmit:"
echo "  Task 15: BMI, Chromosome 15 (NaN error at block 49)"
echo "  Task 24: HbA1c, Chromosome 2 (TIMEOUT after 6 hours)"
echo "  Task 37: HbA1c, Chromosome 15 (NaN error at block 49)"
echo ""

# Check pheno list exists
if [ ! -f "$PHENO_LIST" ]; then
    echo "ERROR: Pheno list file not found: $PHENO_LIST"
    exit 1
fi

# For timeout job (task 24), use longer time limit
# For NaN error jobs (tasks 15, 37), use same time limit but may need investigation
TIMEOUT_TASK=24
CHR15_TASKS="15,37"

echo "=========================================="
echo "Submitting jobs with updated parameters..."
echo "=========================================="
echo ""

# Submit timeout job with increased time limit (12 hours instead of 6)
echo "Submitting task 24 (HbA1c chr2) with 12-hour time limit..."
sbatch \
  --array=24 \
  --time=12:00:00 \
  --cpus-per-task=8 \
  --mem=50G \
  --partition=short \
  --job-name=regenie_prediabetes_perchr_retry \
  --output=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_prediabetes_perchr_retry_%A_%a.out \
  --error=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_prediabetes_perchr_retry_%A_%a.err \
  --export=PHENO_LIST_FILE="$PHENO_LIST" \
  "$REGENIE_SCRIPT"

TIMEOUT_EXIT=$?

# Submit chromosome 15 jobs (tasks 15 and 37) - these failed with NaN error
# Note: The NaN error suggests a data quality issue. We'll resubmit to see if it's transient,
# but may need to investigate block 49 variants if it fails again.
echo ""
echo "Submitting tasks 15,37 (chr15) - investigating NaN error..."
sbatch \
  --array=15,37 \
  --time=6:00:00 \
  --cpus-per-task=8 \
  --mem=50G \
  --partition=short \
  --job-name=regenie_prediabetes_perchr_retry \
  --output=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_prediabetes_perchr_retry_%A_%a.out \
  --error=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_prediabetes_perchr_retry_%A_%a.err \
  --export=PHENO_LIST_FILE="$PHENO_LIST" \
  "$REGENIE_SCRIPT"

CHR15_EXIT=$?

if [ $TIMEOUT_EXIT -eq 0 ] && [ $CHR15_EXIT -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "✓ All jobs submitted successfully"
    echo "=========================================="
    echo ""
    echo "Summary:"
    echo "  ✓ Task 24 (HbA1c chr2): Resubmitted with 12-hour time limit"
    echo "  ✓ Tasks 15,37 (chr15): Resubmitted (investigating NaN error)"
    echo ""
    echo "Monitor jobs:"
    echo "  squeue -u \$USER -p short | grep regenie_prediabetes_perchr_retry"
    echo ""
    echo "Note on chromosome 15 NaN error:"
    echo "  - Error occurred at block 49/6920 on chromosome 15"
    echo "  - This suggests a numerical issue (possibly zero variance variant)"
    echo "  - If it fails again, may need to investigate variants in block 49"
    echo "  - Block 49 corresponds to approximately variants 19,201-19,600 (49 * 400 block size)"
    echo ""
else
    echo ""
    echo "✗ Some jobs failed to submit"
    echo "  Timeout job exit code: $TIMEOUT_EXIT"
    echo "  Chr15 jobs exit code: $CHR15_EXIT"
    exit 1
fi
