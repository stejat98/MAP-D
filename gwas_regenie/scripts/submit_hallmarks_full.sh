#!/bin/bash
# Submit REGENIE jobs for the hallmark traits in the held-out (non-Olink) sample.
# Full cohort only.

HALLMARKS=("BMI" "HbA1c" "TRIG_HDL_RATIO")

PIPE="${REGENIE_PIPELINE_ROOT:-/path/to/project/regenie_pipeline}"
SCRIPT_DIR="${PIPE}/scripts"
REGENIE_SCRIPT="${SCRIPT_DIR}/run_regenie_array.sh"

echo "=========================================="
echo "Submitting REGENIE jobs for hallmarks"
echo "=========================================="
echo "Hallmarks: ${HALLMARKS[@]}"
echo ""

# Build the pheno list consumed by the array job
PHENO_LIST="${SCRIPT_DIR}/pheno_list_hallmarks.txt"
> "$PHENO_LIST"

JOB_COUNT=0
for PHENO in "${HALLMARKS[@]}"; do
    INPUT_PATH="${PIPE}/inputs/full/hallmarks_heldout/${PHENO}"
    if [ -f "${INPUT_PATH}/pheno.txt" ] && [ -f "${INPUT_PATH}/covar.txt" ]; then
        echo "full|hallmarks_heldout|${PHENO}" >> "$PHENO_LIST"
        ((JOB_COUNT++))
        echo "  Added: ${PHENO}"
    else
        echo "  WARNING: Files not found for ${PHENO}, skipping..."
    fi
done

echo ""
echo "Total jobs to submit: $JOB_COUNT"
echo "Pheno list file: $PHENO_LIST"
echo ""

if [ "$JOB_COUNT" -eq 0 ]; then
    echo "Nothing to submit. Run generate_inputs.R first."
    exit 1
fi

# The held-out sample is large (~400K), so these run on the long partition.
sbatch \
  --array=1-${JOB_COUNT} \
  --time=10-00:00:00 \
  --cpus-per-task=8 \
  --mem=50G \
  --partition=long \
  --job-name=regenie_hallmarks_full \
  --output=${PIPE}/slurm/regenie_hallmarks_full_%A_%a.out \
  --error=${PIPE}/slurm/regenie_hallmarks_full_%A_%a.err \
  --export=PHENO_LIST_FILE="$PHENO_LIST" \
  "$REGENIE_SCRIPT"

if [ $? -eq 0 ]; then
    echo "  ✓ Hallmark jobs submitted"
else
    echo "  ✗ Failed to submit hallmark jobs"
    exit 1
fi

echo ""
echo "Monitor jobs:  squeue -u $USER -p long"
echo "Logs:          ${PIPE}/slurm/"
echo ""
