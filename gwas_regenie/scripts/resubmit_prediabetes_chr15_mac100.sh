#!/bin/bash
# Resubmit chromosome 15 prediabetes jobs with MAC 100
# MAC 100 matches Step 1 prefiltering threshold (consistent with pipeline)

SCRIPT_DIR="/n/groups/patel/sivateja/regenie_pipeline/scripts"
REGENIE_SCRIPT="${SCRIPT_DIR}/run_regenie_array_per_chr.sh"
PHENO_LIST="${SCRIPT_DIR}/pheno_list_prediabetes_per_chr.txt"

echo "=========================================="
echo "Resubmitting Chromosome 15 Prediabetes Jobs"
echo "WITH MAC 100 (matches Step 1 threshold)"
echo "=========================================="
echo ""
echo "Jobs to resubmit:"
echo "  Task 15: BMI, Chromosome 15"
echo "  Task 37: HbA1c, Chromosome 15"
echo ""
echo "MAC Threshold: 100"
echo "  - Matches Step 1 prefiltering threshold (consistent with pipeline)"
echo "  - Standard UK Biobank practice for prefiltering"
echo "  - For prediabetes (37K samples): MAC 100 ≈ MAF 0.27%"
echo "  - More stringent than typical Step 2 (MAC 5-20), but principled"
echo "  - Should prevent NaN errors from low-MAC variants"
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
  --job-name=regenie_prediabetes_perchr_chr15_mac100 \
  --output=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_prediabetes_perchr_chr15_mac100_%A_%a.out \
  --error=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_prediabetes_perchr_chr15_mac100_%A_%a.err \
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
    echo "  ✓ Tasks 15,37 (chr15 prediabetes): Resubmitted with MAC 100"
    echo "  ✓ MAC 100 matches Step 1 prefiltering (consistent)"
    echo "  ✓ All other analyses continue using MAC 5 (consistent)"
    echo ""
    echo "Monitor jobs:"
    echo "  squeue -u \$USER -p short | grep regenie_prediabetes_perchr_chr15_mac100"
    echo ""
    echo "Note:"
    echo "  - MAC 100 is stringent but matches Step 1 threshold"
    echo "  - If this still fails, may need to exclude block 49 variants"
    echo "  - MAC 100 excludes ~1-2% of variants (more than MAC 20)"
    echo ""
else
    echo ""
    echo "✗ Failed to submit jobs (exit code: $EXIT_CODE)"
    exit 1
fi
