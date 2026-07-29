#!/bin/bash
# Submit REGENIE jobs for full strata split by chromosome (to avoid timeouts)
# Each chromosome gets its own job (~13 hours each, well within 2-day limit)
# For remaining hallmarks: HDL, LDL, TRIG_HDL_RATIO, systolic_BP, diastolic_BP

SCRIPT_DIR="/n/groups/patel/sivateja/regenie_pipeline/scripts"
REGENIE_SCRIPT="${SCRIPT_DIR}/run_regenie_array_per_chr.sh"
PHENO_LIST="${SCRIPT_DIR}/pheno_list_full_per_chr_remaining.txt"

echo "=========================================="
echo "Submitting REGENIE jobs for full strata (split by chromosome)"
echo "=========================================="
echo "Phenotypes: HDL, LDL, TRIG_HDL_RATIO, systolic_BP, diastolic_BP"
echo "Chromosomes: 1-22 (110 jobs total: 5 phenotypes × 22 chromosomes)"
echo ""

# Check pheno list exists
if [ ! -f "$PHENO_LIST" ]; then
    echo "ERROR: Pheno list file not found: $PHENO_LIST"
    exit 1
fi

TOTAL_JOBS=$(wc -l < "$PHENO_LIST")
echo "Total jobs to submit: $TOTAL_JOBS"
echo "Pheno list file: $PHENO_LIST"
echo ""

echo "=========================================="
echo "Submitting jobs..."
echo "=========================================="
echo "  Partition: medium (2-day limit sufficient for ~13 hour jobs)"
echo "  Time limit: 2-00:00:00 (48 hours - safe margin for ~13 hour jobs)"
echo "  Strategy: Split by chromosome to avoid timeout issues"
echo ""

sbatch \
  --array=1-${TOTAL_JOBS} \
  --time=2-00:00:00 \
  --cpus-per-task=8 \
  --mem=50G \
  --partition=medium \
  --job-name=regenie_full_perchr_remaining \
  --output=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_full_perchr_remaining_%A_%a.out \
  --error=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_full_perchr_remaining_%A_%a.err \
  --export=PHENO_LIST_FILE="$PHENO_LIST" \
  "$REGENIE_SCRIPT"

if [ $? -eq 0 ]; then
    echo "  ✓ Full stratum jobs submitted (split by chromosome)"
    echo ""
    echo "Summary:"
    echo "  ✓ Total jobs: ${TOTAL_JOBS} (5 phenotypes × 22 chromosomes)"
    echo "  ✓ Each job processes one chromosome (~13 hours)"
    echo "  ✓ Time limit: 2 days (plenty of margin)"
    echo "  ✓ Jobs can run in parallel (faster completion)"
    echo ""
    echo "Monitor jobs:"
    echo "  squeue -u st320 -p medium | grep regenie_full_perchr_remaining"
    echo ""
    echo "Logs: /n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_full_perchr_remaining_*.out"
else
    echo "  ✗ Failed to submit jobs"
    exit 1
fi
