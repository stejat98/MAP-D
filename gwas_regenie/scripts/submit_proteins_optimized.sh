#!/bin/bash
# Submit REGENIE jobs for proteins with fairshare optimization
# Submits in smaller batches to manage fairshare score

SCRIPT_DIR="/n/groups/patel/sivateja/regenie_pipeline/scripts"
REGENIE_SCRIPT="${SCRIPT_DIR}/run_regenie_array.sh"

# Get list of all proteins
PROTEINS_DIR="/n/groups/patel/sivateja/regenie_pipeline/inputs/full/proteins_only"
PROTEINS=($(ls "$PROTEINS_DIR" | sort))
TOTAL_PROTEINS=${#PROTEINS[@]}

STRATA=("diabetes" "prediabetes" "full")

echo "=========================================="
echo "Submitting REGENIE jobs for proteins (optimized)"
echo "=========================================="
echo "Total proteins: $TOTAL_PROTEINS"
echo "Strata: ${STRATA[@]}"
echo ""

# Create pheno list file
PHENO_LIST="${SCRIPT_DIR}/pheno_list_proteins.txt"
> "$PHENO_LIST"

JOB_COUNT=0
for STRATUM in "${STRATA[@]}"; do
    for PROTEIN in "${PROTEINS[@]}"; do
        INPUT_PATH="/n/groups/patel/sivateja/regenie_pipeline/inputs/${STRATUM}/proteins_only/${PROTEIN}"
        if [ -f "${INPUT_PATH}/pheno.txt" ] && [ -f "${INPUT_PATH}/covar.txt" ]; then
            echo "${STRATUM}|proteins_only|${PROTEIN}" >> "$PHENO_LIST"
            ((JOB_COUNT++))
        fi
    done
done

echo "Total jobs: $JOB_COUNT"
echo "Pheno list file: $PHENO_LIST"
echo ""

# Fairshare optimization: Submit in smaller batches
# Batch size: 50-100 jobs at a time to avoid overwhelming scheduler
BATCH_SIZE=50

# Calculate batches for diabetes/prediabetes
DIAB_PREDIAB_COUNT=$(grep -E "^diabetes|^prediabetes" "$PHENO_LIST" | wc -l)
DIAB_PREDIAB_BATCHES=$(( (DIAB_PREDIAB_COUNT + BATCH_SIZE - 1) / BATCH_SIZE ))

# Calculate batches for full stratum
FULL_START=$(grep -n "^full" "$PHENO_LIST" | head -1 | cut -d: -f1)
FULL_END=$(wc -l < "$PHENO_LIST")
FULL_COUNT=$((FULL_END - FULL_START + 1))
FULL_BATCHES=$(( (FULL_COUNT + BATCH_SIZE - 1) / BATCH_SIZE ))

echo "Batch strategy:"
echo "  Diabetes/Prediabetes: $DIAB_PREDIAB_COUNT jobs in $DIAB_PREDIAB_BATCHES batches"
echo "  Full: $FULL_COUNT jobs in $FULL_BATCHES batches"
echo "  Batch size: $BATCH_SIZE jobs"
echo ""

# Submit diabetes/prediabetes in batches
if [ $DIAB_PREDIAB_COUNT -gt 0 ]; then
    echo "=========================================="
    echo "Submitting diabetes/prediabetes batches..."
    echo "=========================================="
    
    PREV_JOB=""
    for ((BATCH=1; BATCH<=DIAB_PREDIAB_BATCHES; BATCH++)); do
        BATCH_START=$(( (BATCH - 1) * BATCH_SIZE + 1 ))
        BATCH_END=$(( BATCH * BATCH_SIZE ))
        if [ $BATCH_END -gt $DIAB_PREDIAB_COUNT ]; then
            BATCH_END=$DIAB_PREDIAB_COUNT
        fi
        
        echo "  Batch $BATCH: jobs $BATCH_START-$BATCH_END"
        
        SBATCH_CMD="sbatch"
        if [ -n "$PREV_JOB" ]; then
            # Add dependency on previous batch to stagger submissions
            SBATCH_CMD="$SBATCH_CMD --dependency=afterok:$PREV_JOB"
        fi
        
        JOB_ID=$($SBATCH_CMD \
          --array=${BATCH_START}-${BATCH_END} \
          --time=1-12:00:00 \
          --cpus-per-task=8 \
          --mem=50G \
          --partition=medium \
          --job-name=regenie_prot_${BATCH} \
          --output=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_proteins_%A_%a.out \
          --error=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_proteins_%A_%a.err \
          --export=PHENO_LIST_FILE="$PHENO_LIST" \
          "$REGENIE_SCRIPT" 2>&1 | grep -oP 'Submitted batch job \K[0-9]+')
        
        if [ -n "$JOB_ID" ]; then
            echo "    ✓ Submitted batch job $JOB_ID"
            PREV_JOB="$JOB_ID"
            # Small delay between batches to help fairshare
            sleep 2
        else
            echo "    ✗ Failed to submit batch $BATCH"
        fi
    done
fi

# Submit full stratum in batches (after diabetes/prediabetes)
if [ $FULL_COUNT -gt 0 ]; then
    echo ""
    echo "=========================================="
    echo "Submitting full stratum batches..."
    echo "=========================================="
    
    PREV_JOB=""
    for ((BATCH=1; BATCH<=FULL_BATCHES; BATCH++)); do
        BATCH_START=$(( FULL_START + (BATCH - 1) * BATCH_SIZE ))
        BATCH_END=$(( FULL_START + BATCH * BATCH_SIZE - 1 ))
        if [ $BATCH_END -gt $FULL_END ]; then
            BATCH_END=$FULL_END
        fi
        
        echo "  Batch $BATCH: jobs $BATCH_START-$BATCH_END"
        
        SBATCH_CMD="sbatch"
        if [ -n "$PREV_JOB" ]; then
            SBATCH_CMD="$SBATCH_CMD --dependency=afterok:$PREV_JOB"
        fi
        
        JOB_ID=$($SBATCH_CMD \
          --array=${BATCH_START}-${BATCH_END} \
          --time=2-00:00:00 \
          --cpus-per-task=8 \
          --mem=50G \
          --partition=medium \
          --job-name=regenie_prot_full_${BATCH} \
          --output=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_proteins_full_%A_%a.out \
          --error=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_proteins_full_%A_%a.err \
          --export=PHENO_LIST_FILE="$PHENO_LIST" \
          "$REGENIE_SCRIPT" 2>&1 | grep -oP 'Submitted batch job \K[0-9]+')
        
        if [ -n "$JOB_ID" ]; then
            echo "    ✓ Submitted batch job $JOB_ID"
            PREV_JOB="$JOB_ID"
            sleep 2
        else
            echo "    ✗ Failed to submit batch $BATCH"
        fi
    done
fi

echo ""
echo "=========================================="
echo "Job submission complete!"
echo "=========================================="
echo "Total jobs: $JOB_COUNT"
echo "  - Diabetes/Prediabetes: $DIAB_PREDIAB_COUNT (in $DIAB_PREDIAB_BATCHES batches)"
echo "  - Full: $FULL_COUNT (in $FULL_BATCHES batches)"
echo ""
echo "Fairshare optimization:"
echo "  - Batched submissions (50 jobs/batch)"
echo "  - Staggered with dependencies"
echo "  - Shorter time limits for smaller samples"
echo ""
echo "Monitor with: squeue -u st320"
echo ""
