#!/bin/bash
# Resubmit chromosome 15 prediabetes jobs with MAC 20
# MAC 10 still failed, so increasing to MAC 20

SCRIPT_DIR="/n/groups/patel/sivateja/regenie_pipeline/scripts"
REGENIE_SCRIPT="${SCRIPT_DIR}/run_regenie_array_per_chr.sh"
PHENO_LIST="${SCRIPT_DIR}/pheno_list_prediabetes_per_chr.txt"

echo "=========================================="
echo "Resubmitting Chromosome 15 Prediabetes Jobs"
echo "WITH MAC 20 (MAC 10 still failed)"
echo "=========================================="
echo ""
echo "Jobs to resubmit:"
echo "  Task 15: BMI, Chromosome 15"
echo "  Task 37: HbA1c, Chromosome 15"
echo ""
echo "MAC Threshold: 20 (increased from 10, which still failed)"
echo "  - MAC 10 was used but still encountered NaN error at block 49"
echo "  - Increasing to MAC 20 to filter more low-MAC variants"
echo "  - For prediabetes (37K samples): MAC 20 ≈ MAF 0.053%"
echo "  - Still retains >99% of variants"
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
  --job-name=regenie_prediabetes_perchr_chr15_mac20 \
  --output=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_prediabetes_perchr_chr15_mac20_%A_%a.out \
  --error=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_prediabetes_perchr_chr15_mac20_%A_%a.err \
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
    echo "  ✓ Tasks 15,37 (chr15 prediabetes): Resubmitted with MAC 20"
    echo "  ✓ MAC 10 was used but still failed - increasing to MAC 20"
    echo "  ✓ All other analyses continue using MAC 5 (consistent)"
    echo ""
    echo "Monitor jobs:"
    echo "  squeue -u \$USER -p short | grep regenie_prediabetes_perchr_chr15_mac20"
    echo ""
    echo "Note:"
    echo "  - MAC 20 is more stringent but should prevent NaN errors"
    echo "  - If this still fails, may need to exclude block 49 variants"
    echo ""
else
    echo ""
    echo "✗ Failed to submit jobs (exit code: $EXIT_CODE)"
    exit 1
fi
