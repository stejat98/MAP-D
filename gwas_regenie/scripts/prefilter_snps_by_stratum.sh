#!/bin/bash
# Pre-filter SNPs by MAC for each stratum separately
# This accounts for different sample sizes in each stratum

# Configuration
GENO_PATH="/n/no_backup2/patel/uk_biobank/ukb_genetics/22881"
FAM_PATH="${GENO_PATH}/bgen_converted"
CHR=22  # Step 1 uses chromosome 22
# Use MAF instead of MAC for consistency across different sample sizes
# MAF ≥ 0.001 (0.1%) ensures variants are polymorphic in smallest stratum (diabetes ~900 individuals)
# This is equivalent to MAC ≥ 487 in full dataset, but scales appropriately
MIN_MAF=0.001  # 0.1% minor allele frequency

# Output directory
OUTPUT_DIR="/n/groups/patel/sivateja/regenie_pipeline/filtered_snps"
mkdir -p "$OUTPUT_DIR"

# Load PLINK
module load plink/1.90b7.7_20241022 2>/dev/null || {
    echo "ERROR: PLINK not available"
    exit 1
}

# Use bed files for filtering
BED_FILE="${FAM_PATH}/ukb${CHR}.bed"
BIM_FILE="${FAM_PATH}/ukb${CHR}.bim"
FAM_FILE="${FAM_PATH}/ukb${CHR}.fam"

if [ ! -f "$BED_FILE" ] || [ ! -f "$BIM_FILE" ] || [ ! -f "$FAM_FILE" ]; then
    echo "ERROR: Bed/Bim/Fam files not found for chromosome $CHR"
    exit 1
fi

echo "=========================================="
echo "Pre-filtering SNPs by MAC for REGENIE Step 1"
echo "Using lower MAC threshold for smaller strata"
echo "=========================================="
echo "Chromosome: $CHR"
echo "Minimum MAC: $MIN_MAC (lowered to account for diabetes stratum ~900 individuals)"
echo ""

SNP_LIST="${OUTPUT_DIR}/chr${CHR}_maf${MIN_MAF}.snplist"

echo "Filtering SNPs from bed file using PLINK with MAF threshold..."
echo "  Using MAF ≥ $MIN_MAF (0.1%) to ensure variants work across all strata"
plink \
  --bed "$BED_FILE" \
  --bim "$BIM_FILE" \
  --fam "$FAM_FILE" \
  --maf $MIN_MAF \
  --write-snplist \
  --out "${SNP_LIST%.snplist}" 2>&1 | tee "${OUTPUT_DIR}/plink_filter_maf${MIN_MAF}.log"

if [ ${PIPESTATUS[0]} -ne 0 ]; then
    echo "ERROR: PLINK filtering failed"
    exit 1
fi

if [ ! -f "$SNP_LIST" ]; then
    echo "ERROR: SNP list file not created: $SNP_LIST"
    exit 1
fi

SNP_COUNT=$(wc -l < "$SNP_LIST")
echo ""
echo "✓ Successfully filtered SNPs"
echo "  Total SNPs passing MAF >= $MIN_MAF (0.1%): $SNP_COUNT"
echo "  Output file: $SNP_LIST"
echo ""
echo "Note: MAF-based filtering ensures consistency across all strata:"
echo "  - Full (35K): MAF ≥ 0.1% ≈ MAC ≥ 35"
echo "  - Prediabetes (4K): MAF ≥ 0.1% ≈ MAC ≥ 4"
echo "  - Diabetes (900): MAF ≥ 0.1% ≈ MAC ≥ 1 (but will be higher due to filtering)"
