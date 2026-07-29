#!/bin/bash
# Submit REGENIE jobs for all hallmarks across all strata
# Submits diabetes and prediabetes first (smaller, faster), then full

# Hallmarks (held-out sample GWAS)
HALLMARKS=("BMI" "HbA1c" "HDL" "LDL" "TRIG_HDL_RATIO" "systolic_BP" "diastolic_BP" "ALST")

# Strata (submit in order: diabetes, prediabetes, then full)
STRATA=("diabetes" "prediabetes" "full")

SCRIPT_DIR="/n/groups/patel/sivateja/regenie_pipeline/scripts"
REGENIE_SCRIPT="${SCRIPT_DIR}/run_regenie_array.sh"

echo "=========================================="
echo "Submitting REGENIE jobs for hallmarks"
echo "=========================================="
echo "Hallmarks: ${HALLMARKS[@]}"
echo "Strata: ${STRATA[@]}"
echo ""

# Create pheno list file for array job submission
PHENO_LIST="${SCRIPT_DIR}/pheno_list_hallmarks.txt"
> "$PHENO_LIST"  # Clear/create file

JOB_COUNT=0
for STRATUM in "${STRATA[@]}"; do
    for PHENO in "${HALLMARKS[@]}"; do
        # Check if input files exist
        INPUT_PATH="/n/groups/patel/sivateja/regenie_pipeline/inputs/${STRATUM}/hallmarks_heldout/${PHENO}"
        if [ -f "${INPUT_PATH}/pheno.txt" ] && [ -f "${INPUT_PATH}/covar.txt" ]; then
            echo "${STRATUM}|hallmarks_heldout|${PHENO}" >> "$PHENO_LIST"
            ((JOB_COUNT++))
            echo "  Added: ${STRATUM} - ${PHENO}"
        else
            echo "  WARNING: Files not found for ${STRATUM} - ${PHENO}, skipping..."
        fi
    done
done

echo ""
echo "Total jobs to submit: $JOB_COUNT"
echo "Pheno list file: $PHENO_LIST"
echo ""

# Estimate time based on actual test run (diabetes stratum, BMI, job 24864971)
# Test results (Dec 29, 2025): 
#   - Step 1: 3444s (~57.4 min) for chr22, bsize 500, 4 threads, ~8,435 individuals
#   - Step 2: 667s (~11.1 min) for chr22 only, bsize 200, 4 threads
#   - Total wall time: 2:21:09 (2.35 hours)
# Production config: bsize 1000, 8 threads (2x block size, 2x threads = ~4x faster per block)
# Step 1: Processes chr22 only (same as test, but faster with production settings)
# Step 2: Processes all 22 chromosomes (22x more work than test)
#
# Production time estimates (accounting for optimizations):
#   Step 1: 57.4 min / 4 (speedup) = ~14.4 min
#   Step 2: 11.1 min * 22 chromosomes / 2 (thread speedup) = ~122 min (~2 hours)
#   Total diabetes: ~2.5 hours + overhead = ~3-4 hours (use 4:00:00 with buffer)
#
# Sample sizes and scaling (REGENIE scales roughly linearly with sample size):
#   Diabetes:    ~8,435 individuals (baseline)
#   Prediabetes: ~45,000 individuals (~5.3x diabetes)
#   Full:        ~418,000 individuals (~49.5x diabetes)
#
# Revised estimates accounting for size differences and test results:
#   Diabetes (~8,435): ~3-4 hours → use 4:00:00 (fits in short partition, highest priority)
#   Prediabetes (~45K): ~3-4h * 5.3 = ~16-21h → use 1-00:00:00 (1 day with buffer)
#   Full (~418K): ~3-4h * 49.5 = ~149-198h (~6.2-8.2 days) → use 10-00:00:00 (10 days with buffer)
#
# Time limits include buffer for retries, I/O variability, SNP exclusions, and multi-chromosome overhead
# Fairshare optimization: Submit shortest jobs first (diabetes) to build fairshare score

# Fairshare optimization strategy:
# 1. Submit shortest jobs first (diabetes) to build fairshare score
# 2. Use highest priority partition when possible (short > medium > long)
# 3. Use shortest time limit that fits the job
# 4. Submit in separate batches by stratum for better scheduling
#
# Partition priorities: short (PriorityJobFactor=12) > medium (6) > long (4)
# - Diabetes (~3-4h): Use short partition (fits in 12h max, highest priority)
# - Prediabetes (~16-21h): Use medium partition (fits in 5d max, better than long)
# - Full (~6-8d): Use long partition (needs >5d max of medium)

echo "=========================================="
echo "Submitting jobs (optimized for fairshare)..."
echo "=========================================="

# Submit diabetes first (shortest, highest priority partition)
DIAB_COUNT=$(grep -c "^diabetes" "$PHENO_LIST" 2>/dev/null || echo "0")
if [ $DIAB_COUNT -gt 0 ]; then
    DIAB_START=$(grep -n "^diabetes" "$PHENO_LIST" | head -1 | cut -d: -f1)
    DIAB_END=$(grep -n "^diabetes" "$PHENO_LIST" | tail -1 | cut -d: -f1)
    
    echo "Submitting diabetes jobs (${DIAB_COUNT} jobs)..."
    echo "  Partition: short (PriorityJobFactor=12, highest priority)"
    echo "  Time: 12:00:00 (updated based on actual runtime - worst case ~9.4h, using max short partition limit with buffer)"
    echo "  Strategy: Submit first to build fairshare score"
    echo "  Actual runtime analysis (Dec 29, 2025): Slowest job needs ~9.4h, using 12h limit (max for short partition)"
    
    sbatch \
      --array=${DIAB_START}-${DIAB_END} \
      --time=12:00:00 \
      --cpus-per-task=8 \
      --mem=50G \
      --partition=short \
      --job-name=regenie_hallmarks_diab \
      --output=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_hallmarks_diab_%A_%a.out \
      --error=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_hallmarks_diab_%A_%a.err \
      --export=PHENO_LIST_FILE="$PHENO_LIST" \
      "$REGENIE_SCRIPT"
    
    if [ $? -eq 0 ]; then
        echo "  ✓ Diabetes jobs submitted to short partition"
    else
        echo "  ✗ Failed to submit diabetes jobs"
    fi
else
    echo "  No diabetes jobs to submit"
fi

# Submit prediabetes next (medium priority)
PREDIAB_COUNT=$(grep -c "^prediabetes" "$PHENO_LIST" 2>/dev/null || echo "0")
if [ $PREDIAB_COUNT -gt 0 ]; then
    PREDIAB_START=$(grep -n "^prediabetes" "$PHENO_LIST" | head -1 | cut -d: -f1)
    PREDIAB_END=$(grep -n "^prediabetes" "$PHENO_LIST" | tail -1 | cut -d: -f1)
    
    echo ""
    echo "Submitting prediabetes jobs (${PREDIAB_COUNT} jobs)..."
    echo "  Partition: medium (PriorityJobFactor=6)"
    echo "  Sample size: ~45K individuals (~5.3x diabetes)"
    echo "  Time: 1-00:00:00 (estimated ~16-21 hours, 1 day with buffer for size scaling)"
    echo "  Scaling: diabetes ~3-4h * 5.3x = ~16-21h"
    
    sbatch \
      --array=${PREDIAB_START}-${PREDIAB_END} \
      --time=1-00:00:00 \
      --cpus-per-task=8 \
      --mem=50G \
      --partition=medium \
      --job-name=regenie_hallmarks_prediab \
      --output=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_hallmarks_prediab_%A_%a.out \
      --error=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_hallmarks_prediab_%A_%a.err \
      --export=PHENO_LIST_FILE="$PHENO_LIST" \
      "$REGENIE_SCRIPT"
    
    if [ $? -eq 0 ]; then
        echo "  ✓ Prediabetes jobs submitted"
    else
        echo "  ✗ Failed to submit prediabetes jobs"
    fi
else
    echo ""
    echo "  No prediabetes jobs to submit"
fi

# Submit full stratum last (longest, lowest priority partition)
FULL_START=$(grep -n "^full" "$PHENO_LIST" | head -1 | cut -d: -f1)
FULL_END=$(wc -l < "$PHENO_LIST")
FULL_COUNT=$((FULL_END - FULL_START + 1))

if [ $FULL_COUNT -gt 0 ]; then
    echo ""
    echo "Submitting full stratum jobs (${FULL_COUNT} jobs)..."
    echo "  Partition: long (PriorityJobFactor=4, needed for >5 day jobs)"
    echo "  Sample size: ~418K individuals (~49.5x diabetes)"
    echo "  Time: 10-00:00:00 (estimated ~6-8 days, 10 days with buffer for size scaling)"
    echo "  Scaling: diabetes ~3-4h * 49.5x = ~149-198h (~6.2-8.2 days)"
    echo "  Strategy: Submit after shorter jobs to improve fairshare"
    
    sbatch \
      --array=${FULL_START}-${FULL_END} \
      --time=10-00:00:00 \
      --cpus-per-task=8 \
      --mem=50G \
      --partition=long \
      --job-name=regenie_hallmarks_full \
      --output=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_hallmarks_full_%A_%a.out \
      --error=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_hallmarks_full_%A_%a.err \
      --export=PHENO_LIST_FILE="$PHENO_LIST" \
      "$REGENIE_SCRIPT"
    
    if [ $? -eq 0 ]; then
        echo "  ✓ Full stratum jobs submitted"
    else
        echo "  ✗ Failed to submit full stratum jobs"
    fi
else
    echo ""
    echo "  No full stratum jobs to submit"
fi

echo ""
echo "=========================================="
echo "Job submission complete!"
echo "=========================================="
echo ""
echo "Fairshare optimization summary:"
echo "  ✓ Diabetes: short partition (PriorityJobFactor=12, ~9-10h actual runtime, 12:00:00 limit)"
echo "              Sample: ~8,435 individuals - Submit FIRST to build fairshare"
echo "              Updated Dec 29, 2025: Actual runtime analysis shows ~9.4h worst case, using 12h (max short partition)"
echo "  ✓ Prediabetes: medium partition (PriorityJobFactor=6, ~16-21h, 1-00:00:00 limit)"
echo "                 Sample: ~45K individuals (~5.3x diabetes)"
echo "  ✓ Full: long partition (PriorityJobFactor=4, ~6-8d, 10-00:00:00 limit)"
echo "          Sample: ~418K individuals (~49.5x diabetes) - Submit LAST"
echo ""
echo "Strategy:"
echo "  - Shortest jobs submitted first to build fairshare score"
echo "  - Highest priority partition used when possible"
echo "  - Shortest time limits that fit each job"
echo ""
echo "Monitor jobs:"
echo "  squeue -u st320"
echo "  squeue -u st320 -p short    # Diabetes jobs"
echo "  squeue -u st320 -p medium   # Prediabetes jobs"
echo "  squeue -u st320 -p long     # Full stratum jobs"
echo ""
echo "Logs: /n/groups/patel/sivateja/regenie_pipeline/slurm/"
echo ""
