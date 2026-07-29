#!/bin/bash
# Resubmit chromosome 15 prediabetes jobs
# Option 1: Keep MAC 5 (consistent with earlier analyses) - may still fail
# Option 2: Use MAC 10 (more stringent, fixes issue) - breaks consistency

SCRIPT_DIR="/n/groups/patel/sivateja/regenie_pipeline/scripts"
REGENIE_SCRIPT="${SCRIPT_DIR}/run_regenie_array_per_chr.sh"
PHENO_LIST="${SCRIPT_DIR}/pheno_list_prediabetes_per_chr.txt"

echo "=========================================="
echo "Resubmitting Chromosome 15 Prediabetes Jobs"
echo "=========================================="
echo ""
echo "MAC Threshold Options:"
echo "  1. MAC 5 (default, consistent with earlier analyses)"
echo "  2. MAC 10 (more stringent, fixes NaN error)"
echo ""
read -p "Choose option (1 or 2, default 2): " MAC_OPTION
MAC_OPTION=${MAC_OPTION:-2}

if [ "$MAC_OPTION" = "1" ]; then
    MAC_THRESHOLD=5
    echo "Using MAC 5 (consistent with earlier analyses)"
    echo "WARNING: May still encounter NaN errors"
elif [ "$MAC_OPTION" = "2" ]; then
    MAC_THRESHOLD=10
    echo "Using MAC 10 (more stringent, should fix NaN error)"
else
    echo "Invalid option, using MAC 10"
    MAC_THRESHOLD=10
fi

# Temporarily modify script to use chosen MAC threshold
TEMP_SCRIPT="/tmp/run_regenie_array_per_chr_${MAC_THRESHOLD}.sh"
cp "$REGENIE_SCRIPT" "$TEMP_SCRIPT"
sed -i "s/MIN_MAC_THRESHOLD=.*/MIN_MAC_THRESHOLD=${MAC_THRESHOLD}/" "$TEMP_SCRIPT"

echo ""
echo "Jobs to resubmit:"
echo "  Task 15: BMI, Chromosome 15"
echo "  Task 37: HbA1c, Chromosome 15"
echo ""

# Check pheno list exists
if [ ! -f "$PHENO_LIST" ]; then
    echo "ERROR: Pheno list file not found: $PHENO_LIST"
    exit 1
fi

echo "=========================================="
echo "Submitting jobs with MAC ${MAC_THRESHOLD}..."
echo "=========================================="
echo ""

sbatch \
  --array=15,37 \
  --time=6:00:00 \
  --cpus-per-task=8 \
  --mem=50G \
  --partition=short \
  --job-name=regenie_prediabetes_perchr_mac${MAC_THRESHOLD} \
  --output=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_prediabetes_perchr_mac${MAC_THRESHOLD}_%A_%a.out \
  --error=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_prediabetes_perchr_mac${MAC_THRESHOLD}_%A_%a.err \
  --export=PHENO_LIST_FILE="$PHENO_LIST" \
  "$TEMP_SCRIPT"

EXIT_CODE=$?

if [ $EXIT_CODE -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "✓ Jobs submitted successfully"
    echo "=========================================="
    echo ""
    echo "Summary:"
    echo "  ✓ Tasks 15,37 (chr15): Resubmitted with --minMAC ${MAC_THRESHOLD}"
    echo ""
    if [ "$MAC_THRESHOLD" = "5" ]; then
        echo "Note: Using MAC 5 (consistent with earlier analyses)"
        echo "      If NaN error persists, consider using MAC 10"
    else
        echo "Note: Using MAC 10 (more stringent than earlier analyses)"
        echo "      This should prevent NaN errors but excludes more variants"
    fi
    echo ""
    echo "Monitor jobs:"
    echo "  squeue -u \$USER -p short | grep regenie_prediabetes_perchr_mac${MAC_THRESHOLD}"
    echo ""
else
    echo ""
    echo "✗ Failed to submit jobs (exit code: $EXIT_CODE)"
    rm -f "$TEMP_SCRIPT"
    exit 1
fi

rm -f "$TEMP_SCRIPT"
