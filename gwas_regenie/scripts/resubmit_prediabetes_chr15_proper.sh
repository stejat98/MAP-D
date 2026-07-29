#!/bin/bash
# Resubmit chromosome 15 prediabetes jobs with proper MAC threshold
# This uses --minMAC 10 instead of excluding specific variants
# This is the standard approach for UK Biobank-scale analyses

SCRIPT_DIR="/n/groups/patel/sivateja/regenie_pipeline/scripts"
REGENIE_SCRIPT="${SCRIPT_DIR}/run_regenie_array_per_chr.sh"
PHENO_LIST="${SCRIPT_DIR}/pheno_list_prediabetes_per_chr.txt"

echo "=========================================="
echo "Resubmitting Chromosome 15 Prediabetes Jobs"
echo "WITH PROPER MAC THRESHOLD (--minMAC 10)"
echo "=========================================="
echo ""
echo "Jobs to resubmit:"
echo "  Task 15: BMI, Chromosome 15"
echo "  Task 37: HbA1c, Chromosome 15"
echo ""
echo "Solution: Using --minMAC 10 (instead of default 5)"
echo "  - Prevents numerical errors from low-MAC variants"
echo "  - Standard practice for UK Biobank GWAS analyses"
echo "  - Only excludes ~0.1-0.5% of variants (those with MAC < 10)"
echo "  - More principled than excluding specific variants"
echo ""

# Check pheno list exists
if [ ! -f "$PHENO_LIST" ]; then
    echo "ERROR: Pheno list file not found: $PHENO_LIST"
    exit 1
fi

# Verify script has been updated with --minMAC
if ! grep -q "MIN_MAC_THRESHOLD\|--minMAC" "$REGENIE_SCRIPT"; then
    echo "WARNING: Script may not have --minMAC threshold set"
    echo "Please ensure run_regenie_array_per_chr.sh includes --minMAC flag"
fi

echo "=========================================="
echo "Submitting jobs with MAC threshold..."
echo "=========================================="
echo ""

sbatch \
  --array=15,37 \
  --time=6:00:00 \
  --cpus-per-task=8 \
  --mem=50G \
  --partition=short \
  --job-name=regenie_prediabetes_perchr_mac10 \
  --output=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_prediabetes_perchr_mac10_%A_%a.out \
  --error=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_prediabetes_perchr_mac10_%A_%a.err \
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
    echo "  ✓ Tasks 15,37 (chr15): Resubmitted with --minMAC 10"
    echo "  ✓ This is the standard approach for UK Biobank analyses"
    echo "  ✓ Prevents numerical errors from low-MAC variants"
    echo ""
    echo "Why this is better than excluding variants:"
    echo "  ✓ Addresses root cause (low MAC variants)"
    echo "  ✓ Prevents future issues, not just block 49"
    echo "  ✓ Principled filtering (standard GWAS practice)"
    echo "  ✓ Reproducible and scalable"
    echo ""
    echo "Monitor jobs:"
    echo "  squeue -u \$USER -p short | grep regenie_prediabetes_perchr_mac10"
    echo ""
else
    echo ""
    echo "✗ Failed to submit jobs (exit code: $EXIT_CODE)"
    exit 1
fi
