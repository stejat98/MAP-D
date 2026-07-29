#!/bin/bash
# Modified version of run_regenie_array_per_chr.sh to handle chromosome 15 block 49 issue
# Excludes problematic variants in block 49 that cause NaN errors in prediabetes stratum

# Copy the original script logic but add exclude file for chr15 block 49
SCRIPT_DIR="/n/groups/patel/sivateja/regenie_pipeline/scripts"
ORIGINAL_SCRIPT="${SCRIPT_DIR}/run_regenie_array_per_chr.sh"

# Check if this is chromosome 15 for prediabetes
PHENO_LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$PHENO_LIST_FILE")
STRATUM=$(echo "$PHENO_LINE" | cut -d'|' -f1)
TARGET_CHR=$(echo "$PHENO_LINE" | cut -d'|' -f4)

# Create exclude file for block 49 variants if this is chr15
if [ "$TARGET_CHR" = "15" ] && [ "$STRATUM" = "prediabetes" ]; then
    BLOCK_SIZE=400
    BLOCK_NUM=49
    START_VAR=$(( (BLOCK_NUM - 1) * BLOCK_SIZE + 1 ))
    END_VAR=$(( BLOCK_NUM * BLOCK_SIZE ))
    BIM_FILE="/n/no_backup2/patel/uk_biobank/ukb_genetics/22881/bgen_converted/ukb15.bim"
    EXCLUDE_FILE="/tmp/chr15_block49_exclude_${SLURM_JOB_ID}.txt"
    
    echo "=========================================="
    echo "CHROMOSOME 15 BLOCK 49 FIX"
    echo "=========================================="
    echo "Excluding block 49 variants to avoid NaN error"
    echo "Block 49: variants $START_VAR to $END_VAR"
    echo "Creating exclude file: $EXCLUDE_FILE"
    
    if [ -f "$BIM_FILE" ]; then
        sed -n "${START_VAR},${END_VAR}p" "$BIM_FILE" | awk '{print $2}' > "$EXCLUDE_FILE"
        EXCLUDE_COUNT=$(wc -l < "$EXCLUDE_FILE")
        echo "Excluding $EXCLUDE_COUNT variants from block 49"
        export CHR15_BLOCK49_EXCLUDE="$EXCLUDE_FILE"
    else
        echo "WARNING: BIM file not found, cannot create exclude file"
        export CHR15_BLOCK49_EXCLUDE=""
    fi
    echo ""
fi

# Source the original script
source "$ORIGINAL_SCRIPT"
