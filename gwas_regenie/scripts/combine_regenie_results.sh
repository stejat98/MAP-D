#!/bin/bash
# Combine per-chromosome REGENIE results into single files
# Following UK Biobank GWAS best practices (Neale lab style)
# Usage: combine_regenie_results.sh [phenotype] [stratum]
#   If no arguments provided, combines all BMI and HbA1c results

set -e

# Default: combine the hallmark traits for the full cohort
PHENOTYPES="${1:-BMI HbA1c}"
STRATA="${2:-full}"

# Output directory (scratch for processing, then copy to permanent location)
SCRATCH_ROOT="/n/scratch/users/s/st320/regenie"
RESULTS_ROOT="/n/groups/patel/sivateja/regenie_pipeline/results"

# Create results directory if it doesn't exist
mkdir -p "$RESULTS_ROOT"

echo "=========================================="
echo "Combining REGENIE per-chromosome results"
echo "=========================================="
echo "Phenotypes: $PHENOTYPES"
echo "Strata: $STRATA"
echo ""

for PHENO in $PHENOTYPES; do
    for STRATUM in $STRATA; do
        STEP2_DIR="${SCRATCH_ROOT}/${PHENO}/${STRATUM}/step2"
        OUTPUT_FILE="${RESULTS_ROOT}/${PHENO}_${STRATUM}_all_chr.regenie.gz"
        
        # Check if step2 directory exists
        if [ ! -d "$STEP2_DIR" ]; then
            echo "⚠ WARNING: Directory not found: $STEP2_DIR"
            continue
        fi
        
        # Check if all chromosome files exist and verify naming consistency
        MISSING_CHR=""
        WRONG_NAME_CHR=""
        for CHR in {1..22}; do
            CHR_FILE="${STEP2_DIR}/regenie_step2_${PHENO}_chr${CHR}_${PHENO}.regenie"
            if [ ! -f "$CHR_FILE" ]; then
                MISSING_CHR="${MISSING_CHR} ${CHR}"
            else
                # Verify the file actually contains the correct phenotype name
                # Check that filename contains phenotype name (safety check)
                if [[ ! "$CHR_FILE" =~ ${PHENO} ]]; then
                    WRONG_NAME_CHR="${WRONG_NAME_CHR} ${CHR}"
                fi
            fi
        done
        
        if [ ! -z "$WRONG_NAME_CHR" ]; then
            echo "✗ ERROR: File naming mismatch for ${PHENO} ${STRATUM} (chromosomes:$WRONG_NAME_CHR)"
            echo "  Expected pattern: regenie_step2_${PHENO}_chr*_${PHENO}.regenie"
            continue
        fi
        
        if [ ! -z "$MISSING_CHR" ]; then
            echo "⚠ WARNING: Missing chromosome files for ${PHENO} ${STRATUM}:$MISSING_CHR"
            continue
        fi
        
        echo "Combining ${PHENO} ${STRATUM}..."
        echo "  Input: ${STEP2_DIR}/regenie_step2_${PHENO}_chr*_${PHENO}.regenie"
        echo "  Output: ${OUTPUT_FILE}"
        
        # Verify first file header contains expected phenotype name (safety check)
        FIRST_FILE="${STEP2_DIR}/regenie_step2_${PHENO}_chr1_${PHENO}.regenie"
        if ! grep -q "^CHROM" "$FIRST_FILE"; then
            echo "  ✗ ERROR: First file does not have expected header format"
            continue
        fi
        
        # Combine files: header from first file, then data from all files
        # Sort by chromosome and position (CHROM, GENPOS)
        # Use numeric sort for chromosome and position
        # Write to temp file first, then compress (more reliable for large files)
        TEMP_OUTPUT="${OUTPUT_FILE%.gz}"
        
        (
            # Get header from first chromosome file
            head -1 "$FIRST_FILE"
            
            # Combine all chromosome files (skip headers)
            # Process in order (chr1-22) to ensure proper sorting
            for CHR in {1..22}; do
                CHR_FILE="${STEP2_DIR}/regenie_step2_${PHENO}_chr${CHR}_${PHENO}.regenie"
                # Verify file exists and has data (safety check)
                if [ ! -f "$CHR_FILE" ]; then
                    echo "  ✗ ERROR: Missing chromosome ${CHR} file: $CHR_FILE" >&2
                    continue
                fi
                # Verify file contains the correct chromosome number in data
                # (First data line should have CHROM=$CHR)
                FIRST_DATA_LINE=$(tail -n +2 "$CHR_FILE" | head -1)
                FILE_CHROM=$(echo "$FIRST_DATA_LINE" | awk '{print $1}')
                if [ "$FILE_CHROM" != "$CHR" ]; then
                    echo "  ⚠ WARNING: Chromosome ${CHR} file contains CHROM=${FILE_CHROM} data" >&2
                fi
                tail -n +2 "$CHR_FILE"
            done | sort -k1,1n -k2,2n  # Sort by CHROM (numeric), then GENPOS (numeric)
        ) > "$TEMP_OUTPUT"
        
        # Compress the combined file
        gzip "$TEMP_OUTPUT"
        
        # Verify compression succeeded
        if [ ! -f "$OUTPUT_FILE" ]; then
            echo "  ✗ ERROR: Failed to create compressed output file"
            [ -f "$TEMP_OUTPUT" ] && rm -f "$TEMP_OUTPUT"
            continue
        fi
        
        # Verify output
        if [ -f "$OUTPUT_FILE" ] && [ -s "$OUTPUT_FILE" ]; then
            N_VARIANTS=$(zcat "$OUTPUT_FILE" | tail -n +2 | wc -l)
            FILE_SIZE=$(du -h "$OUTPUT_FILE" | cut -f1)
            echo "  ✓ Success: ${N_VARIANTS} variants, ${FILE_SIZE}"
        else
            echo "  ✗ ERROR: Output file not created or empty"
        fi
        echo ""
    done
done

echo "=========================================="
echo "Combination complete!"
echo "=========================================="
echo ""
echo "Summary of combined files:"
for PHENO in $PHENOTYPES; do
    for STRATUM in $STRATA; do
        OUTPUT_FILE="${RESULTS_ROOT}/${PHENO}_${STRATUM}_all_chr.regenie.gz"
        if [ -f "$OUTPUT_FILE" ]; then
            N_VARIANTS=$(zcat "$OUTPUT_FILE" | tail -n +2 | wc -l)
            FILE_SIZE=$(du -h "$OUTPUT_FILE" | cut -f1)
            echo "  ${PHENO}_${STRATUM}: ${N_VARIANTS} variants, ${FILE_SIZE}"
        fi
    done
done
echo ""
echo "Output directory: ${RESULTS_ROOT}"
echo ""
echo "File format:"
echo "  Columns: CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTRA"
echo "  Sorted by: CHROM (numeric), GENPOS (numeric)"
echo "  Compressed: gzip"
echo ""
echo "To view a file:"
echo "  zcat ${RESULTS_ROOT}/BMI_diabetes_all_chr.regenie.gz | head -20"
echo "  zcat ${RESULTS_ROOT}/BMI_diabetes_all_chr.regenie.gz | grep -w 'rs12345'"
