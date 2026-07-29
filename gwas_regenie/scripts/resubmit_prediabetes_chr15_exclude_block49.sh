#!/bin/bash
# Resubmit chromosome 15 prediabetes jobs with block 49 exclusion
# MAC filtering (5, 10, 100) all failed - issue is specific to block 49 variants

SCRIPT_DIR="/n/groups/patel/sivateja/regenie_pipeline/scripts"
REGENIE_SCRIPT="${SCRIPT_DIR}/run_regenie_array_per_chr.sh"
PHENO_LIST="${SCRIPT_DIR}/pheno_list_prediabetes_per_chr.txt"

echo "=========================================="
echo "Resubmitting Chromosome 15 Prediabetes Jobs"
echo "WITH BLOCK 49 EXCLUSION (MAC filtering failed)"
echo "=========================================="
echo ""
echo "Jobs to resubmit:"
echo "  Task 15: BMI, Chromosome 15"
echo "  Task 37: HbA1c, Chromosome 15"
echo ""
echo "Solution: Exclude block 49 variants (400 variants)"
echo "  - MAC 5, 10, 100 all failed at block 49"
echo "  - MAC filtering is NOT the issue"
echo "  - Block 49 has specific problematic variant(s)"
echo "  - Excluding 400 variants (0.014% of chr15)"
echo "  - Diabetes/full strata passed without exclusion"
echo ""

# Check pheno list exists
if [ ! -f "$PHENO_LIST" ]; then
    echo "ERROR: Pheno list file not found: $PHENO_LIST"
    exit 1
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
  --job-name=regenie_prediabetes_perchr_chr15_exclude49 \
  --output=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_prediabetes_perchr_chr15_exclude49_%A_%a.out \
  --error=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_prediabetes_perchr_chr15_exclude49_%A_%a.err \
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
    echo "  ✓ Tasks 15,37 (chr15 prediabetes): Resubmitted with block 49 exclusion"
    echo "  ✓ MAC filtering (5, 10, 100) all failed - not the issue"
    echo "  ✓ Excluding 400 variants from block 49 (targeted fix)"
    echo "  ✓ All other analyses continue using MAC 5 (consistent)"
    echo ""
    echo "Monitor jobs:"
    echo "  squeue -u \$USER -p short | grep regenie_prediabetes_perchr_chr15_exclude49"
    echo ""
    echo "Note:"
    echo "  - Block 49 exclusion is necessary (MAC filtering didn't work)"
    echo "  - Only 400 variants excluded (0.014% of chr15)"
    echo "  - Minimal impact on analysis"
    echo ""
else
    echo ""
    echo "✗ Failed to submit jobs (exit code: $EXIT_CODE)"
    exit 1
fi
