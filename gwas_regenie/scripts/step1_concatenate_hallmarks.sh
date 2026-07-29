#!/bin/bash
#
# Step 1: Concatenate per-chromosome REGENIE files for hallmarks (full stratum)
#

set -e

SCRATCH_DIR="/n/scratch/users/s/st320/regenie"
RESULTS_DIR="/n/groups/patel/sivateja/regenie_pipeline/results"

PHENOTYPES=("HDL" "LDL" "systolic_BP" "diastolic_BP" "TRIG_HDL_RATIO")
STRATUM="full"

echo "=============================================="
echo "STEP 1: Concatenate Per-Chromosome Files"
echo "=============================================="
echo "Stratum: ${STRATUM}"
echo "Phenotypes: ${PHENOTYPES[*]}"
echo "Output directory: ${RESULTS_DIR}"
echo "=============================================="
echo ""

total=${#PHENOTYPES[@]}
current=0

for PHENO in "${PHENOTYPES[@]}"; do
    current=$((current + 1))
    echo "[${current}/${total}] Processing: ${PHENO}"
    
    INPUT_DIR="${SCRATCH_DIR}/${PHENO}/${STRATUM}/step2"
    OUTPUT_FILE="${RESULTS_DIR}/${PHENO}_${STRATUM}_all_chr.regenie.gz"
    
    # Check input exists
    if [[ ! -d "$INPUT_DIR" ]]; then
        echo "  ERROR: Input directory not found: $INPUT_DIR"
        continue
    fi
    
    # Get first file for header
    first_file="${INPUT_DIR}/regenie_step2_${PHENO}_chr1_${PHENO}.regenie"
    if [[ ! -f "$first_file" ]]; then
        echo "  ERROR: First chromosome file not found: $first_file"
        continue
    fi
    
    echo "  Reading from: ${INPUT_DIR}"
    echo "  Writing to: ${OUTPUT_FILE}"
    
    # Concatenate with progress
    (
        head -1 "$first_file"
        for chr in {1..22}; do
            chr_file="${INPUT_DIR}/regenie_step2_${PHENO}_chr${chr}_${PHENO}.regenie"
            if [[ -f "$chr_file" ]]; then
                tail -n +2 "$chr_file"
                echo "    chr${chr} done" >&2
            else
                echo "    WARNING: chr${chr} file missing" >&2
            fi
        done
    ) | gzip > "$OUTPUT_FILE"
    
    # Report result
    file_size=$(ls -lh "$OUTPUT_FILE" | awk '{print $5}')
    variant_count=$(zcat "$OUTPUT_FILE" | tail -n +2 | wc -l)
    echo "  DONE: ${file_size}, ${variant_count} variants"
    echo ""
done

echo "=============================================="
echo "STEP 1 COMPLETE"
echo "=============================================="
echo "Created files:"
ls -lh "${RESULTS_DIR}"/*_full_all_chr.regenie.gz 2>/dev/null | grep -E "HDL|LDL|systolic|diastolic|TRIG" || echo "  (none found)"
echo ""
echo "Next: Run step2_convert_to_fuma.sh"
