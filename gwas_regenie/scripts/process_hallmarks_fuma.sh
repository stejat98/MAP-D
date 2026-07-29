#!/bin/bash
#
# Process hallmark phenotypes: concatenate per-chr files, convert to FUMA, and filter p<0.1
# For full stratum only (as requested)
#

set -e

SCRATCH_DIR="/n/scratch/users/s/st320/regenie"
RESULTS_DIR="/n/groups/patel/sivateja/regenie_pipeline/results"
SCRIPTS_DIR="/n/groups/patel/sivateja/regenie_pipeline/scripts"

# Hallmarks to process (excluding BMI and HbA1c which are already done)
PHENOTYPES=("HDL" "LDL" "systolic_BP" "diastolic_BP" "TRIG_HDL_RATIO")
STRATUM="full"

echo "=============================================="
echo "Processing Hallmark Phenotypes for FUMA"
echo "Stratum: ${STRATUM}"
echo "=============================================="
echo ""

for PHENO in "${PHENOTYPES[@]}"; do
    echo "=============================================="
    echo "Processing: ${PHENO}"
    echo "=============================================="
    
    INPUT_DIR="${SCRATCH_DIR}/${PHENO}/${STRATUM}/step2"
    OUTPUT_BASE="${RESULTS_DIR}/${PHENO}_${STRATUM}_all_chr"
    
    # Check if input directory exists
    if [[ ! -d "$INPUT_DIR" ]]; then
        echo "  ERROR: Input directory not found: $INPUT_DIR"
        continue
    fi
    
    # Step 1: Concatenate per-chromosome files
    echo "  Step 1: Concatenating per-chromosome files..."
    REGENIE_FILE="${OUTPUT_BASE}.regenie.gz"
    
    # Get header from first file
    first_file=$(ls "${INPUT_DIR}"/*.regenie 2>/dev/null | head -1)
    if [[ -z "$first_file" ]]; then
        echo "  ERROR: No .regenie files found in $INPUT_DIR"
        continue
    fi
    
    # Concatenate: header from first file, then all data (skip headers)
    (
        head -1 "$first_file"
        for chr in {1..22}; do
            chr_file="${INPUT_DIR}/regenie_step2_${PHENO}_chr${chr}_${PHENO}.regenie"
            if [[ -f "$chr_file" ]]; then
                tail -n +2 "$chr_file"
            fi
        done
    ) | gzip > "$REGENIE_FILE"
    
    # Count variants
    variant_count=$(zcat "$REGENIE_FILE" | tail -n +2 | wc -l)
    file_size=$(ls -lh "$REGENIE_FILE" | awk '{print $5}')
    echo "    Created: $(basename $REGENIE_FILE) ($file_size, $variant_count variants)"
    
    # Step 2: Convert to FUMA format
    echo "  Step 2: Converting to FUMA format..."
    python3 "${SCRIPTS_DIR}/convert_to_fuma_format.py" "$REGENIE_FILE" --no-size-check
    
    FUMA_FILE="${OUTPUT_BASE}.fuma.gz"
    fuma_size=$(ls -lh "$FUMA_FILE" | awk '{print $5}')
    echo "    Created: $(basename $FUMA_FILE) ($fuma_size)"
    
    # Step 3: Filter to p < 0.1
    echo "  Step 3: Filtering to P < 0.1..."
    P01_FILE="${OUTPUT_BASE}.p0.1.fuma.gz"
    
    zcat "$FUMA_FILE" | awk 'NR==1 || $6 < 0.1' | gzip > "$P01_FILE"
    
    p01_size=$(ls -lh "$P01_FILE" | awk '{print $5}')
    p01_count=$(zcat "$P01_FILE" | tail -n +2 | wc -l)
    echo "    Created: $(basename $P01_FILE) ($p01_size, $p01_count variants)"
    
    echo ""
done

echo "=============================================="
echo "Summary"
echo "=============================================="
echo ""
echo "Files created in ${RESULTS_DIR}:"
ls -lh "${RESULTS_DIR}"/*_full_*.gz 2>/dev/null | grep -E "HDL|LDL|systolic|diastolic|TRIG"
echo ""
echo "Done!"
