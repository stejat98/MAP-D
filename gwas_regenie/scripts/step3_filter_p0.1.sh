#!/bin/bash
#
# Step 3: Filter FUMA files to P < 0.1 for hallmarks (full stratum)
#

set -e

RESULTS_DIR="/n/groups/patel/sivateja/regenie_pipeline/results"

PHENOTYPES=("HDL" "LDL" "systolic_BP" "diastolic_BP" "TRIG_HDL_RATIO")
STRATUM="full"

echo "=============================================="
echo "STEP 3: Filter to P < 0.1"
echo "=============================================="
echo "Phenotypes: ${PHENOTYPES[*]}"
echo "=============================================="
echo ""

total=${#PHENOTYPES[@]}
current=0

for PHENO in "${PHENOTYPES[@]}"; do
    current=$((current + 1))
    echo "[${current}/${total}] Filtering: ${PHENO}"
    
    INPUT_FILE="${RESULTS_DIR}/${PHENO}_${STRATUM}_all_chr.fuma.gz"
    OUTPUT_FILE="${RESULTS_DIR}/${PHENO}_${STRATUM}_all_chr.p0.1.fuma.gz"
    
    # Check input exists
    if [[ ! -f "$INPUT_FILE" ]]; then
        echo "  ERROR: Input file not found: $INPUT_FILE"
        echo "  Run step2_convert_to_fuma.sh first"
        continue
    fi
    
    input_size=$(ls -lh "$INPUT_FILE" | awk '{print $5}')
    input_count=$(zcat "$INPUT_FILE" | tail -n +2 | wc -l)
    echo "  Input: $(basename $INPUT_FILE) (${input_size}, ${input_count} variants)"
    
    # Filter to P < 0.1 (column 6 is P-value)
    echo "  Filtering P < 0.1..."
    zcat "$INPUT_FILE" | awk 'NR==1 || $6 < 0.1' | gzip > "$OUTPUT_FILE"
    
    # Report result
    output_size=$(ls -lh "$OUTPUT_FILE" | awk '{print $5}')
    output_count=$(zcat "$OUTPUT_FILE" | tail -n +2 | wc -l)
    reduction=$(echo "scale=1; (1 - ($output_count / $input_count)) * 100" | bc)
    
    echo "  Output: $(basename $OUTPUT_FILE) (${output_size}, ${output_count} variants)"
    echo "  Reduction: ${reduction}%"
    echo ""
done

echo "=============================================="
echo "STEP 3 COMPLETE"
echo "=============================================="
echo ""
echo "All p0.1 FUMA files created:"
ls -lh "${RESULTS_DIR}"/*_full_all_chr.p0.1.fuma.gz 2>/dev/null | grep -E "BMI|HbA1c|TRIG_HDL_RATIO" || echo "  (none found)"
echo ""
echo "=============================================="
echo "ALL STEPS COMPLETE!"
echo "=============================================="
