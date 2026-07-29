#!/bin/bash
# Resubmit chromosome 15 prediabetes jobs with block 49 exclusion fix
# This version uses the modified script that excludes problematic block 49 variants

SCRIPT_DIR="/n/groups/patel/sivateja/regenie_pipeline/scripts"
REGENIE_SCRIPT="${SCRIPT_DIR}/run_regenie_array_per_chr.sh"
PHENO_LIST="${SCRIPT_DIR}/pheno_list_prediabetes_per_chr.txt"

echo "=========================================="
echo "Resubmitting Chromosome 15 Prediabetes Jobs"
echo "WITH BLOCK 49 EXCLUSION FIX"
echo "=========================================="
echo ""
echo "Jobs to resubmit:"
echo "  Task 15: BMI, Chromosome 15 (will exclude block 49 variants)"
echo "  Task 37: HbA1c, Chromosome 15 (will exclude block 49 variants)"
echo ""
echo "Fix: Excluding 400 variants in block 49 (variants 19,201-19,600)"
echo "     that cause NaN errors in prediabetes stratum"
echo ""

# Check pheno list exists
if [ ! -f "$PHENO_LIST" ]; then
    echo "ERROR: Pheno list file not found: $PHENO_LIST"
    exit 1
fi

# Check that script has been modified
if ! grep -q "chr15_block49_exclude" "$REGENIE_SCRIPT"; then
    echo "ERROR: Script has not been modified to exclude block 49 variants"
    echo "Please ensure run_regenie_array_per_chr.sh includes the block 49 exclusion logic"
    exit 1
fi

echo "=========================================="
echo "Submitting jobs with block 49 exclusion..."
echo "=========================================="
echo ""

sbatch \
  --array=15,37 \
  --time=6:00:00 \
  --cpus-per-task=8 \
  --mem=50G \
  --partition=short \
  --job-name=regenie_prediabetes_perchr_chr15fix \
  --output=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_prediabetes_perchr_chr15fix_%A_%a.out \
  --error=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_prediabetes_perchr_chr15fix_%A_%a.err \
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
    echo "  ✓ Tasks 15,37 (chr15): Resubmitted with block 49 exclusion"
    echo "  ✓ Block 49 variants (19,201-19,600) will be excluded"
    echo "  ✓ This should prevent the NaN error"
    echo ""
    echo "Monitor jobs:"
    echo "  squeue -u \$USER -p short | grep regenie_prediabetes_perchr_chr15fix"
    echo ""
    echo "Note:"
    echo "  - 400 variants from block 49 will be excluded from analysis"
    echo "  - These variants cause NaN errors specifically in prediabetes stratum"
    echo "  - Diabetes and full strata completed successfully without exclusion"
    echo ""
else
    echo ""
    echo "✗ Failed to submit jobs (exit code: $EXIT_CODE)"
    exit 1
fi
