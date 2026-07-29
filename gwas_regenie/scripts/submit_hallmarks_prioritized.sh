#!/bin/bash
# Submit REGENIE jobs prioritizing BMI and HbA1c first
# Then submit remaining hallmarks
# Submits diabetes first (shortest), then prediabetes, then full

# Priority hallmarks (submit first)
PRIORITY_HALLMARKS=("BMI" "HbA1c")

# Remaining hallmarks
OTHER_HALLMARKS=("HDL" "LDL" "TRIG_HDL_RATIO" "systolic_BP" "diastolic_BP")

# All hallmarks
ALL_HALLMARKS=("${PRIORITY_HALLMARKS[@]}" "${OTHER_HALLMARKS[@]}")

# Strata (submit in order: diabetes, prediabetes, then full)
STRATA=("diabetes" "prediabetes" "full")

SCRIPT_DIR="/n/groups/patel/sivateja/regenie_pipeline/scripts"
REGENIE_SCRIPT="${SCRIPT_DIR}/run_regenie_array.sh"

echo "=========================================="
echo "Submitting REGENIE jobs for hallmarks"
echo "PRIORITY: BMI and HbA1c first"
echo "=========================================="
echo "Priority hallmarks: ${PRIORITY_HALLMARKS[@]}"
echo "Other hallmarks: ${OTHER_HALLMARKS[@]}"
echo "Strata: ${STRATA[@]}"
echo ""

# Create pheno list files
PRIORITY_LIST="${SCRIPT_DIR}/pheno_list_priority.txt"
OTHER_LIST="${SCRIPT_DIR}/pheno_list_other.txt"
> "$PRIORITY_LIST"  # Clear/create file
> "$OTHER_LIST"      # Clear/create file

PRIORITY_COUNT=0
OTHER_COUNT=0

# Build priority list (BMI and HbA1c)
for STRATUM in "${STRATA[@]}"; do
    for PHENO in "${PRIORITY_HALLMARKS[@]}"; do
        INPUT_PATH="/n/groups/patel/sivateja/regenie_pipeline/inputs/${STRATUM}/hallmarks_heldout/${PHENO}"
        if [ -f "${INPUT_PATH}/pheno.txt" ] && [ -f "${INPUT_PATH}/covar.txt" ]; then
            echo "${STRATUM}|hallmarks_heldout|${PHENO}" >> "$PRIORITY_LIST"
            ((PRIORITY_COUNT++))
            echo "  [PRIORITY] Added: ${STRATUM} - ${PHENO}"
        fi
    done
done

# Build other hallmarks list
for STRATUM in "${STRATA[@]}"; do
    for PHENO in "${OTHER_HALLMARKS[@]}"; do
        INPUT_PATH="/n/groups/patel/sivateja/regenie_pipeline/inputs/${STRATUM}/hallmarks_heldout/${PHENO}"
        if [ -f "${INPUT_PATH}/pheno.txt" ] && [ -f "${INPUT_PATH}/covar.txt" ]; then
            echo "${STRATUM}|hallmarks_heldout|${PHENO}" >> "$OTHER_LIST"
            ((OTHER_COUNT++))
            echo "  [OTHER] Added: ${STRATUM} - ${PHENO}"
        fi
    done
done

echo ""
echo "Priority jobs: $PRIORITY_COUNT"
echo "Other jobs: $OTHER_COUNT"
echo "Total jobs: $((PRIORITY_COUNT + OTHER_COUNT))"
echo ""

# Function to submit jobs for a given list
submit_stratum_jobs() {
    local LIST_FILE="$1"
    local STRATUM_NAME="$2"
    local TIME_LIMIT="$3"
    local PARTITION="$4"
    local JOB_NAME_PREFIX="$5"
    
    local COUNT=$(wc -l < "$LIST_FILE" 2>/dev/null || echo "0")
    if [ $COUNT -eq 0 ]; then
        echo "  No ${STRATUM_NAME} jobs in this list"
        return
    fi
    
    local START=1
    local END=$COUNT
    
    echo "Submitting ${STRATUM_NAME} jobs (${COUNT} jobs from ${LIST_FILE})..."
    echo "  Partition: $PARTITION"
    echo "  Time: $TIME_LIMIT"
    
    sbatch \
      --array=${START}-${END} \
      --time=$TIME_LIMIT \
      --cpus-per-task=8 \
      --mem=50G \
      --partition=$PARTITION \
      --job-name=${JOB_NAME_PREFIX}_${STRATUM_NAME} \
      --output=/n/groups/patel/sivateja/regenie_pipeline/slurm/${JOB_NAME_PREFIX}_${STRATUM_NAME}_%A_%a.out \
      --error=/n/groups/patel/sivateja/regenie_pipeline/slurm/${JOB_NAME_PREFIX}_${STRATUM_NAME}_%A_%a.err \
      --export=PHENO_LIST_FILE="$LIST_FILE" \
      "$REGENIE_SCRIPT"
    
    if [ $? -eq 0 ]; then
        echo "  ✓ ${STRATUM_NAME} jobs submitted"
    else
        echo "  ✗ Failed to submit ${STRATUM_NAME} jobs"
    fi
}

echo "=========================================="
echo "PHASE 1: Priority Hallmarks (BMI, HbA1c)"
echo "=========================================="

# Submit priority hallmarks - diabetes first
DIAB_PRIORITY=$(grep -c "^diabetes" "$PRIORITY_LIST" 2>/dev/null || echo "0")
if [ $DIAB_PRIORITY -gt 0 ]; then
    DIAB_START=$(grep -n "^diabetes" "$PRIORITY_LIST" | head -1 | cut -d: -f1)
    DIAB_END=$(grep -n "^diabetes" "$PRIORITY_LIST" | tail -1 | cut -d: -f1)
    
    echo "Submitting PRIORITY diabetes jobs (${DIAB_PRIORITY} jobs)..."
    echo "  Partition: short (highest priority)"
    echo "  Time: 11:59:00"
    
    sbatch \
      --array=${DIAB_START}-${DIAB_END} \
      --time=11:59:00 \
      --cpus-per-task=8 \
      --mem=50G \
      --partition=short \
      --job-name=regenie_priority_diab \
      --output=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_priority_diab_%A_%a.out \
      --error=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_priority_diab_%A_%a.err \
      --export=PHENO_LIST_FILE="$PRIORITY_LIST" \
      "$REGENIE_SCRIPT"
    
    if [ $? -eq 0 ]; then
        echo "  ✓ Priority diabetes jobs submitted"
    fi
fi

# Submit priority hallmarks - prediabetes
PREDIAB_PRIORITY=$(grep -c "^prediabetes" "$PRIORITY_LIST" 2>/dev/null || echo "0")
if [ $PREDIAB_PRIORITY -gt 0 ]; then
    PREDIAB_START=$(grep -n "^prediabetes" "$PRIORITY_LIST" | head -1 | cut -d: -f1)
    PREDIAB_END=$(grep -n "^prediabetes" "$PRIORITY_LIST" | tail -1 | cut -d: -f1)
    
    echo ""
    echo "Submitting PRIORITY prediabetes jobs (${PREDIAB_PRIORITY} jobs)..."
    echo "  Partition: medium"
    echo "  Time: 3-00:00:00"
    
    sbatch \
      --array=${PREDIAB_START}-${PREDIAB_END} \
      --time=3-00:00:00 \
      --cpus-per-task=8 \
      --mem=50G \
      --partition=medium \
      --job-name=regenie_priority_prediab \
      --output=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_priority_prediab_%A_%a.out \
      --error=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_priority_prediab_%A_%a.err \
      --export=PHENO_LIST_FILE="$PRIORITY_LIST" \
      "$REGENIE_SCRIPT"
    
    if [ $? -eq 0 ]; then
        echo "  ✓ Priority prediabetes jobs submitted"
    fi
fi

# Submit priority hallmarks - full
FULL_PRIORITY_START=$(grep -n "^full" "$PRIORITY_LIST" | head -1 | cut -d: -f1 2>/dev/null)
FULL_PRIORITY_END=$(grep -n "^full" "$PRIORITY_LIST" | tail -1 | cut -d: -f1 2>/dev/null)
FULL_PRIORITY_COUNT=$((FULL_PRIORITY_END - FULL_PRIORITY_START + 1))

if [ ! -z "$FULL_PRIORITY_START" ] && [ $FULL_PRIORITY_COUNT -gt 0 ]; then
    echo ""
    echo "Submitting PRIORITY full stratum jobs (${FULL_PRIORITY_COUNT} jobs)..."
    echo "  Partition: long"
    echo "  Time: 20-00:00:00"
    
    sbatch \
      --array=${FULL_PRIORITY_START}-${FULL_PRIORITY_END} \
      --time=20-00:00:00 \
      --cpus-per-task=8 \
      --mem=50G \
      --partition=long \
      --job-name=regenie_priority_full \
      --output=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_priority_full_%A_%a.out \
      --error=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_priority_full_%A_%a.err \
      --export=PHENO_LIST_FILE="$PRIORITY_LIST" \
      "$REGENIE_SCRIPT"
    
    if [ $? -eq 0 ]; then
        echo "  ✓ Priority full stratum jobs submitted"
    fi
fi

echo ""
echo "=========================================="
echo "PHASE 2: Other Hallmarks"
echo "=========================================="
echo "Waiting 30 seconds before submitting other hallmarks..."
sleep 30

# Submit other hallmarks - diabetes
DIAB_OTHER=$(grep -c "^diabetes" "$OTHER_LIST" 2>/dev/null || echo "0")
if [ $DIAB_OTHER -gt 0 ]; then
    DIAB_START=$(grep -n "^diabetes" "$OTHER_LIST" | head -1 | cut -d: -f1)
    DIAB_END=$(grep -n "^diabetes" "$OTHER_LIST" | tail -1 | cut -d: -f1)
    
    echo "Submitting other diabetes jobs (${DIAB_OTHER} jobs)..."
    echo "  Partition: short"
    echo "  Time: 11:59:00"
    
    sbatch \
      --array=${DIAB_START}-${DIAB_END} \
      --time=11:59:00 \
      --cpus-per-task=8 \
      --mem=50G \
      --partition=short \
      --job-name=regenie_other_diab \
      --output=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_other_diab_%A_%a.out \
      --error=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_other_diab_%A_%a.err \
      --export=PHENO_LIST_FILE="$OTHER_LIST" \
      "$REGENIE_SCRIPT"
    
    if [ $? -eq 0 ]; then
        echo "  ✓ Other diabetes jobs submitted"
    fi
fi

# Submit other hallmarks - prediabetes
PREDIAB_OTHER=$(grep -c "^prediabetes" "$OTHER_LIST" 2>/dev/null || echo "0")
if [ $PREDIAB_OTHER -gt 0 ]; then
    PREDIAB_START=$(grep -n "^prediabetes" "$OTHER_LIST" | head -1 | cut -d: -f1)
    PREDIAB_END=$(grep -n "^prediabetes" "$OTHER_LIST" | tail -1 | cut -d: -f1)
    
    echo ""
    echo "Submitting other prediabetes jobs (${PREDIAB_OTHER} jobs)..."
    echo "  Partition: medium"
    echo "  Time: 3-00:00:00"
    
    sbatch \
      --array=${PREDIAB_START}-${PREDIAB_END} \
      --time=3-00:00:00 \
      --cpus-per-task=8 \
      --mem=50G \
      --partition=medium \
      --job-name=regenie_other_prediab \
      --output=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_other_prediab_%A_%a.out \
      --error=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_other_prediab_%A_%a.err \
      --export=PHENO_LIST_FILE="$OTHER_LIST" \
      "$REGENIE_SCRIPT"
    
    if [ $? -eq 0 ]; then
        echo "  ✓ Other prediabetes jobs submitted"
    fi
fi

# Submit other hallmarks - full
FULL_OTHER_START=$(grep -n "^full" "$OTHER_LIST" | head -1 | cut -d: -f1 2>/dev/null)
FULL_OTHER_END=$(grep -n "^full" "$OTHER_LIST" | tail -1 | cut -d: -f1 2>/dev/null)
FULL_OTHER_COUNT=$((FULL_OTHER_END - FULL_OTHER_START + 1))

if [ ! -z "$FULL_OTHER_START" ] && [ $FULL_OTHER_COUNT -gt 0 ]; then
    echo ""
    echo "Submitting other full stratum jobs (${FULL_OTHER_COUNT} jobs)..."
    echo "  Partition: long"
    echo "  Time: 20-00:00:00"
    
    sbatch \
      --array=${FULL_OTHER_START}-${FULL_OTHER_END} \
      --time=20-00:00:00 \
      --cpus-per-task=8 \
      --mem=50G \
      --partition=long \
      --job-name=regenie_other_full \
      --output=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_other_full_%A_%a.out \
      --error=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_other_full_%A_%a.err \
      --export=PHENO_LIST_FILE="$OTHER_LIST" \
      "$REGENIE_SCRIPT"
    
    if [ $? -eq 0 ]; then
        echo "  ✓ Other full stratum jobs submitted"
    fi
fi

echo ""
echo "=========================================="
echo "Job Submission Complete!"
echo "=========================================="
echo ""
echo "Submission Summary:"
echo "  Priority (BMI, HbA1c):"
echo "    - Diabetes: ${DIAB_PRIORITY} jobs (short partition)"
echo "    - Prediabetes: ${PREDIAB_PRIORITY} jobs (medium partition)"
echo "    - Full: ${FULL_PRIORITY_COUNT} jobs (long partition)"
echo ""
echo "  Other hallmarks:"
echo "    - Diabetes: ${DIAB_OTHER} jobs (short partition)"
echo "    - Prediabetes: ${PREDIAB_OTHER} jobs (medium partition)"
echo "    - Full: ${FULL_OTHER_COUNT} jobs (long partition)"
echo ""
echo "Monitor jobs:"
echo "  squeue -u st320 | grep regenie_priority"
echo "  squeue -u st320 | grep regenie_other"
echo ""
