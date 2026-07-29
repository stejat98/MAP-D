#!/bin/bash
# Resubmit chromosome 15 prediabetes jobs with MAC 10
# This is a targeted fix for the NaN error issue

SCRIPT_DIR="/n/groups/patel/sivateja/regenie_pipeline/scripts"
REGENIE_SCRIPT="${SCRIPT_DIR}/run_regenie_array_per_chr.sh"
PHENO_LIST="${SCRIPT_DIR}/pheno_list_prediabetes_per_chr.txt"

echo "=========================================="
echo "Resubmitting Chromosome 15 Prediabetes Jobs"
echo "WITH MAC 10 (targeted fix for NaN error)"
echo "=========================================="
echo ""
echo "Jobs to resubmit:"
echo "  Task 15: BMI, Chromosome 15"
echo "  Task 37: HbA1c, Chromosome 15"
echo ""
echo "MAC Threshold: 10 (increased from default 5 for chr15 prediabetes)"
echo "  - Prevents NaN errors from low-MAC variants"
echo "  - Targeted fix for problematic chr15 prediabetes jobs"
echo "  - All other analyses continue using MAC 5 (consistent)"
echo ""

# Check pheno list exists
if [ ! -f "$PHENO_LIST" ]; then
    echo "ERROR: Pheno list file not found: $PHENO_LIST"
    exit 1
fi

# Verify script has been updated
if ! grep -q "TARGET_CHR.*15.*prediabetes" "$REGENIE_SCRIPT"; then
    echo "WARNING: Script may not have chr15 prediabetes MAC 10 logic"
    echo "Please ensure run_regenie_array_per_chr.sh sets MAC 10 for chr15 prediabetes"
fi

echo "=========================================="
echo "Submitting jobs..."
echo "=========================================="
echo ""

sbatch \
  --array=15,37 \
  --time=6:00:00 \
  --cpus-per-task=8 \
  --mem=50G \
  --partition=short \
  --job-name=regenie_prediabetes_perchr_chr15_mac10 \
  --output=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_prediabetes_perchr_chr15_mac10_%A_%a.out \
  --error=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_prediabetes_perchr_chr15_mac10_%A_%a.err \
  --export=PHENO_LIST_FILE="$PHENO_LIST" \
  "$REGENIE_SCRIPT"

EXIT_CODE=$?

if [ $EXIT_CODE -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "✓ Jobs submitted successfully"
    echo "=========================================="
    echo ""
    echo "Summary:"
    echo "  ✓ Tasks 15,37 (chr15 prediabetes): Resubmitted with MAC 10"
    echo "  ✓ MAC 10 is targeted fix for NaN error issue"
    echo "  ✓ All other analyses continue using MAC 5 (consistent)"
    echo ""
    echo "Monitor jobs:"
    echo "  squeue -u \$USER -p short | grep regenie_prediabetes_perchr_chr15_mac10"
    echo ""
    echo "Note:"
    echo "  - MAC 10 is used only for chr15 prediabetes (exception)"
    echo "  - All other jobs continue with MAC 5 (consistent with earlier analyses)"
    echo "  - This should prevent the NaN error while maintaining consistency"
    echo ""
else
    echo ""
    echo "✗ Failed to submit jobs (exit code: $EXIT_CODE)"
    exit 1
fi
