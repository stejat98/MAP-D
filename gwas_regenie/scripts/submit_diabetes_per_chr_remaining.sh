#!/bin/bash
# Submit REGENIE jobs for diabetes strata split by chromosome (to avoid timeouts)
# Each chromosome gets its own job (~2.67 hours each, well within 12-hour limit)
# Uses SHORT queue for higher priority and faster start
# For remaining hallmarks: HDL, LDL, TRIG_HDL_RATIO, systolic_BP, diastolic_BP

SCRIPT_DIR="/n/groups/patel/sivateja/regenie_pipeline/scripts"
REGENIE_SCRIPT="${SCRIPT_DIR}/run_regenie_array_per_chr.sh"
PHENO_LIST="${SCRIPT_DIR}/pheno_list_diabetes_per_chr_remaining.txt"

echo "=========================================="
echo "Submitting REGENIE jobs for diabetes strata (split by chromosome)"
echo "=========================================="
echo "Phenotypes: HDL, LDL, TRIG_HDL_RATIO, systolic_BP, diastolic_BP"
echo "Chromosomes: 1-22 (110 jobs total: 5 phenotypes × 22 chromosomes)"
echo "Partition: SHORT (higher priority, faster start)"
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
echo "  Partition: short (PriorityJobFactor=12, highest priority)"
echo "  Time limit: 6:00:00 (6 hours - safe margin for ~2.67 hour jobs)"
echo "  Strategy: Split by chromosome to avoid timeout issues"
echo "  Benefits: Higher priority, faster start, better fairshare"
echo ""

sbatch \
  --array=1-${TOTAL_JOBS} \
  --time=6:00:00 \
  --cpus-per-task=8 \
  --mem=50G \
  --partition=short \
  --job-name=regenie_diabetes_perchr_remaining \
  --output=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_diabetes_perchr_remaining_%A_%a.out \
  --error=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_diabetes_perchr_remaining_%A_%a.err \
  --export=PHENO_LIST_FILE="$PHENO_LIST" \
  "$REGENIE_SCRIPT"

if [ $? -eq 0 ]; then
    echo "  ✓ Diabetes jobs submitted (split by chromosome, short queue)"
    echo ""
    echo "Summary:"
    echo "  ✓ Total jobs: ${TOTAL_JOBS} (5 phenotypes × 22 chromosomes)"
    echo "  ✓ Each job processes one chromosome (~2.67 hours)"
    echo "  ✓ Time limit: 6 hours (plenty of margin)"
    echo "  ✓ Partition: short (higher priority, faster start)"
    echo "  ✓ Jobs can run in parallel (faster completion)"
    echo ""
    echo "Monitor jobs:"
    echo "  squeue -u st320 -p short | grep regenie_diabetes_perchr_remaining"
    echo ""
    echo "Logs: /n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_diabetes_perchr_remaining_*.out"
else
    echo "  ✗ Failed to submit jobs"
    exit 1
fi
