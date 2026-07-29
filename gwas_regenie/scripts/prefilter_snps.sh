#!/bin/bash
# Pre-filter SNPs by MAC to avoid low-variance errors in REGENIE step 1
# Creates a SNP list file that can be used with --extract in REGENIE

# Configuration
GENO_PATH="/n/no_backup2/patel/uk_biobank/ukb_genetics/22881"
FAM_PATH="${GENO_PATH}/bgen_converted"
CHR=22  # Step 1 uses chromosome 22
MIN_MAC=100  # Minimum minor allele count (recommended: 100 for UK Biobank)

# Output file
OUTPUT_DIR="/n/groups/patel/sivateja/regenie_pipeline/filtered_snps"
mkdir -p "$OUTPUT_DIR"
SNP_LIST="${OUTPUT_DIR}/chr${CHR}_mac${MIN_MAC}.snplist"

echo "=========================================="
echo "Pre-filtering SNPs by MAC for REGENIE Step 1"
echo "=========================================="
echo "Chromosome: $CHR"
echo "Minimum MAC: $MIN_MAC"
echo "Output: $SNP_LIST"
echo ""

# Use PLINK 1.9 to filter bed files by MAC
# Note: We use bed files since PLINK2 is not available, and REGENIE step 2 uses bed format
# The SNP IDs should match between BGEN and bed formats for the same chromosome
module load plink/1.90b7.7_20241022 2>/dev/null || {
    echo "ERROR: PLINK not available. Please load plink module."
    exit 1
}

# Use bed files for filtering (REGENIE step 2 uses bed format, so SNP IDs will match)
BED_FILE="${FAM_PATH}/ukb${CHR}.bed"
BIM_FILE="${FAM_PATH}/ukb${CHR}.bim"
FAM_FILE="${FAM_PATH}/ukb${CHR}.fam"

if [ ! -f "$BED_FILE" ] || [ ! -f "$BIM_FILE" ] || [ ! -f "$FAM_FILE" ]; then
    echo "ERROR: Bed/Bim/Fam files not found for chromosome $CHR"
    echo "  Expected: $BED_FILE"
    echo "  Expected: $BIM_FILE"
    echo "  Expected: $FAM_FILE"
    exit 1
fi

echo "Filtering SNPs from bed file using PLINK..."
echo "  Input: $BED_FILE"
echo "  Minimum MAC: $MIN_MAC"
echo ""
echo "Note: BED files were converted from BGEN files, so SNP IDs should match."
echo "  The filtered list will work with REGENIE --extract when reading BGEN files."
echo ""

plink \
  --bed "$BED_FILE" \
  --bim "$BIM_FILE" \
  --fam "$FAM_FILE" \
  --mac $MIN_MAC \
  --write-snplist \
  --out "${SNP_LIST%.snplist}" 2>&1 | tee "${OUTPUT_DIR}/plink_filter.log"

if [ ${PIPESTATUS[0]} -ne 0 ]; then
    echo "ERROR: PLINK2 filtering failed"
    exit 1
fi

# Check if output was created
if [ ! -f "$SNP_LIST" ]; then
    echo "ERROR: SNP list file not created: $SNP_LIST"
    exit 1
fi

SNP_COUNT=$(wc -l < "$SNP_LIST")
echo ""
echo "✓ Successfully filtered SNPs"
echo "  Total SNPs passing MAC >= $MIN_MAC: $SNP_COUNT"
echo "  Output file: $SNP_LIST"
echo ""
echo "This file can be used with --extract in REGENIE step 1"
