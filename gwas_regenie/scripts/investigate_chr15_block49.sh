#!/bin/bash
# Investigate chromosome 15 block 49 issue
# The NaN error suggests a numerical problem - possibly zero variance variants

echo "=========================================="
echo "Chromosome 15 Block 49 Investigation"
echo "=========================================="
echo ""
echo "Error: Chi Square parameter was -nan at block 49/6920"
echo "Block size: 400 variants per block"
echo "Block 49 contains approximately variants 19,201-19,600"
echo ""

BIM_FILE="/n/no_backup2/patel/uk_biobank/ukb_genetics/22881/bgen_converted/ukb15.bim"
BLOCK_SIZE=400
BLOCK_NUM=49

if [ ! -f "$BIM_FILE" ]; then
    echo "ERROR: BIM file not found: $BIM_FILE"
    exit 1
fi

TOTAL_VARIANTS=$(wc -l < "$BIM_FILE")
echo "Total variants in chromosome 15: $TOTAL_VARIANTS"
echo "Total blocks: $((TOTAL_VARIANTS / BLOCK_SIZE))"
echo ""

# Calculate variant range for block 49
START_VAR=$(( (BLOCK_NUM - 1) * BLOCK_SIZE + 1 ))
END_VAR=$(( BLOCK_NUM * BLOCK_SIZE ))

echo "Block 49 variant range: $START_VAR to $END_VAR"
echo ""
echo "Variants in block 49:"
sed -n "${START_VAR},${END_VAR}p" "$BIM_FILE" | head -20
echo "..."
sed -n "${START_VAR},${END_VAR}p" "$BIM_FILE" | tail -5
echo ""

# Check for potential issues
echo "=========================================="
echo "Checking for potential data quality issues..."
echo "=========================================="

# Check if there are any variants with unusual patterns
echo ""
echo "Checking variant positions in block 49:"
sed -n "${START_VAR},${END_VAR}p" "$BIM_FILE" | awk '{print $4}' | sort -n | head -5
echo "..."
sed -n "${START_VAR},${END_VAR}p" "$BIM_FILE" | awk '{print $4}' | sort -n | tail -5

echo ""
echo "=========================================="
echo "Recommendations:"
echo "=========================================="
echo "1. The NaN error typically occurs when:"
echo "   - A variant has zero variance (all samples have same genotype)"
echo "   - Division by zero in test statistic calculation"
echo "   - Numerical instability in chi-square calculation"
echo ""
echo "2. If resubmission fails again, consider:"
echo "   - Checking variant MAF in block 49"
echo "   - Verifying genotype data quality for these variants"
echo "   - Using --minMAC flag to filter low MAC variants"
echo "   - Checking if other strata had similar issues"
echo ""
echo "3. REGENIE uses --minMAC 5 by default, but block 49 may have"
echo "   variants that pass MAC filter but cause numerical issues"
echo ""
