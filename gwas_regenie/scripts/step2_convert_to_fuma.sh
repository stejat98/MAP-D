#!/bin/bash
#
# Step 2: Convert REGENIE files to FUMA format for hallmarks (full stratum)
#

set -e

RESULTS_DIR="/n/groups/patel/sivateja/regenie_pipeline/results"
SCRIPTS_DIR="/n/groups/patel/sivateja/regenie_pipeline/scripts"

PHENOTYPES=("HDL" "LDL" "systolic_BP" "diastolic_BP" "TRIG_HDL_RATIO")
STRATUM="full"

echo "=============================================="
echo "STEP 2: Convert to FUMA Format"
echo "=============================================="
echo "Phenotypes: ${PHENOTYPES[*]}"
echo "=============================================="
echo ""

total=${#PHENOTYPES[@]}
current=0

for PHENO in "${PHENOTYPES[@]}"; do
    current=$((current + 1))
    echo "[${current}/${total}] Converting: ${PHENO}"
    
    INPUT_FILE="${RESULTS_DIR}/${PHENO}_${STRATUM}_all_chr.regenie.gz"
    OUTPUT_FILE="${RESULTS_DIR}/${PHENO}_${STRATUM}_all_chr.fuma.gz"
    
    # Check input exists
    if [[ ! -f "$INPUT_FILE" ]]; then
        echo "  ERROR: Input file not found: $INPUT_FILE"
        echo "  Run step1_concatenate_hallmarks.sh first"
        continue
    fi
    
    echo "  Input: $(basename $INPUT_FILE)"
    
    # Convert using Python script
    python3 "${SCRIPTS_DIR}/convert_to_fuma_format.py" "$INPUT_FILE" --no-size-check
    
    # Report result
    if [[ -f "$OUTPUT_FILE" ]]; then
        file_size=$(ls -lh "$OUTPUT_FILE" | awk '{print $5}')
        variant_count=$(zcat "$OUTPUT_FILE" | tail -n +2 | wc -l)
        echo "  DONE: ${file_size}, ${variant_count} variants"
    else
        echo "  ERROR: Output file not created"
    fi
    echo ""
done

echo "=============================================="
echo "STEP 2 COMPLETE"
echo "=============================================="
echo "Created files:"
ls -lh "${RESULTS_DIR}"/*_full_all_chr.fuma.gz 2>/dev/null | grep -E "HDL|LDL|systolic|diastolic|TRIG" || echo "  (none found)"
echo ""
echo "Next: Run step3_filter_p0.1.sh"
