#!/bin/bash
# Test version of REGENIE that runs on subset of data
# Uses smaller block sizes and fewer threads for testing

set -e

STRATUM="$1"
SAMPLE_TYPE="$2"
PHENO_NAME="$3"
INPUT_DIR="$4"
OUTDIR_STEP1="$5"
OUTDIR_STEP2="$6"

echo "  Stratum: $STRATUM"
echo "  Sample type: $SAMPLE_TYPE"
echo "  Phenotype: $PHENO_NAME"
echo "  Input dir: $INPUT_DIR"

# Paths
CONDA_ENV="${CONDA_ENV:-/path/to/project/.conda/envs/regenie_env}"
# IMPORTANT: Use ONLY genetic data from /n/no_backup2/patel/uk_biobank/ukb_genetics/22881
# to ensure ID matching across all input files (phenotype, genotype, sample files)
# BGEN files: ${GENO_PATH}/ukb_imp_chr*_v3.bgen (for Step 1)
# BED/BIM/FAM files: ${FAM_PATH}/ukb*.{bed,bim,fam} (for Step 2, converted from BGEN)
GENO_PATH="/n/no_backup2/patel/uk_biobank/ukb_genetics/22881"
FAM_PATH="${GENO_PATH}/bgen_converted"

# Create output directories
TMPDIR_STEP1="${OUTDIR_STEP1}/tmpdir_${SLURM_JOB_ID:-$$}"
mkdir -p "${OUTDIR_STEP1}" "${OUTDIR_STEP2}" "${TMPDIR_STEP1}"

# Check input files
PHENO_FILE="${INPUT_DIR}/pheno.txt"
COVAR_FILE="${INPUT_DIR}/covar.txt"

if [ ! -f "$PHENO_FILE" ]; then
    echo "ERROR: Phenotype file not found: $PHENO_FILE"
    exit 1
fi

if [ ! -f "$COVAR_FILE" ]; then
    echo "ERROR: Covariate file not found: $COVAR_FILE"
    exit 1
fi

# Count individuals
N_INDIV=$(tail -n +2 "$PHENO_FILE" | wc -l)
echo "  Number of individuals: $N_INDIV"

# Conda setup (IMPORTANT: Set environment variables before activating conda)
export HOME=/path/to/project
export CONDA_NO_PLUGINS=true
export CONDA_PKGS_DIRS=/path/to/project/.conda/pkgs
export CONDA_ENVS_PATH=/path/to/project/.conda/envs

module load conda/miniforge3/24.11.3-0
eval "$(conda shell.bash hook)"
conda activate "${CONDA_ENV}"

# Create sample file from FAM file (required for BGEN format, ensures ID matching)
# The sample file must match the IDs in the BGEN file exactly
# Using FAM file from ${FAM_PATH} ensures consistency with BED files used in Step 2
SAMPLE_FILE="${TMPDIR_STEP1}/sample_file.sample"
if [ ! -f "${FAM_PATH}/ukb22.fam" ]; then
    echo "ERROR: FAM file not found: ${FAM_PATH}/ukb22.fam"
    exit 1
fi
awk 'BEGIN {print "ID_1 ID_2 missing sex"; print "0 0 0 D"} {print $1, $2, "0", $5}' "${FAM_PATH}/ukb22.fam" > "$SAMPLE_FILE"
echo "  Created sample file from ${FAM_PATH}/ukb22.fam: $SAMPLE_FILE"

# Use pre-filtered SNP list to avoid low-variance SNPs
FILTERED_SNP_LIST="/path/to/project/regenie_pipeline/filtered_snps/chr22_maf0.001.snplist"
if [ ! -f "$FILTERED_SNP_LIST" ]; then
    FILTERED_SNP_LIST="/path/to/project/regenie_pipeline/filtered_snps/chr22_mac50.snplist"
fi
if [ ! -f "$FILTERED_SNP_LIST" ]; then
    FILTERED_SNP_LIST="/path/to/project/regenie_pipeline/filtered_snps/chr22_mac100.snplist"
fi

# Create exclude file for problematic low-variance SNPs that slip through
# (can happen when MAC filter is applied to full dataset but analysis uses subset)
EXCLUDE_FILE="${TMPDIR_STEP1}/exclude_snps.txt"
touch "$EXCLUDE_FILE"

if [ ! -f "$FILTERED_SNP_LIST" ]; then
    echo "WARNING: Pre-filtered SNP list not found, proceeding without filtering"
    echo "  This may cause low-variance SNP errors. Run prefilter_snps.sh first."
    EXTRACT_FLAG=""
else
    echo "  Using pre-filtered SNP list: $FILTERED_SNP_LIST"
    SNP_COUNT=$(wc -l < "$FILTERED_SNP_LIST")
    echo "  SNPs in filtered list: $SNP_COUNT"
    EXTRACT_FLAG="--extract ${FILTERED_SNP_LIST}"
fi

# Pre-scan and exclude low-variance SNPs for this specific analysis subset
# This is more efficient than retrying REGENIE multiple times
# Only do this if we have a filtered SNP list to compare against
if [ ! -z "$EXTRACT_FLAG" ] && [ -f "$FILTERED_SNP_LIST" ]; then
    echo "  Pre-scanning for low-variance SNPs in analysis subset..."
    # Extract sample IDs from phenotype file (first column, skip header)
    SAMPLE_ID_FILE="${TMPDIR_STEP1}/analysis_samples.txt"
    tail -n +2 "$PHENO_FILE" | awk '{print $1, $1}' > "$SAMPLE_ID_FILE"
    N_ANALYSIS_SAMPLES=$(wc -l < "$SAMPLE_ID_FILE")
    echo "  Analysis subset contains $N_ANALYSIS_SAMPLES samples"
    
    # Use PLINK to check MAC on the analysis subset
    # This will identify SNPs that become monomorphic/low-variance in this specific subset
    PLINK_CHECK_OUT="${TMPDIR_STEP1}/plink_mac_check"
    CHR_BED="${FAM_PATH}/ukb22.bed"
    CHR_BIM="${FAM_PATH}/ukb22.bim"
    CHR_FAM="${FAM_PATH}/ukb22.fam"
    
    if [ -f "$CHR_BED" ] && [ -f "$CHR_BIM" ] && [ -f "$CHR_FAM" ]; then
    # Try to load PLINK module first, then check if it's available
    # Load PLINK module for pre-scanning (needed to identify low-variance SNPs upfront)
    module load plink/1.90b7.7_20241022 2>/dev/null || true
    
    # Check if PLINK is now available
    if command -v plink &> /dev/null; then
        echo "  Checking MAC on analysis subset using PLINK..."
        # Use MAC >= 100 (recommended for UK Biobank, matches prefilter_snps.sh)
        # This is conservative but ensures SNPs remain polymorphic after residualization
        # Even though the filtered list was created with MAF 0.001 or MAC 50/100 on full dataset,
        # when analyzing subsets, SNPs may have lower MAC
        # Using MAC >= 100 matches the standard UK Biobank REGENIE practice
        MIN_MAC=100
        if [ ! -z "$EXTRACT_FLAG" ]; then
            plink \
              --bed "$CHR_BED" \
              --bim "$CHR_BIM" \
              --fam "$CHR_FAM" \
              --keep "$SAMPLE_ID_FILE" \
              --extract "$FILTERED_SNP_LIST" \
              --mac $MIN_MAC \
              --write-snplist \
              --out "$PLINK_CHECK_OUT" > "${PLINK_CHECK_OUT}.log" 2>&1 || true
        else
            plink \
              --bed "$CHR_BED" \
              --bim "$CHR_BIM" \
              --fam "$CHR_FAM" \
              --keep "$SAMPLE_ID_FILE" \
              --mac $MIN_MAC \
              --write-snplist \
              --out "$PLINK_CHECK_OUT" > "${PLINK_CHECK_OUT}.log" 2>&1 || true
        fi
        
        if [ -f "${PLINK_CHECK_OUT}.snplist" ]; then
            # Find SNPs that are in the filtered list but NOT in the PLINK output (low MAC in subset)
            if [ ! -z "$EXTRACT_FLAG" ]; then
                LOW_VAR_SNPS="${TMPDIR_STEP1}/low_variance_snps.txt"
                # Get SNPs in filtered list but not in PLINK output (comm -23: lines only in first file)
                comm -23 <(sort "$FILTERED_SNP_LIST") <(sort "${PLINK_CHECK_OUT}.snplist") > "$LOW_VAR_SNPS" 2>/dev/null || true
                
                if [ -s "$LOW_VAR_SNPS" ]; then
                    LOW_VAR_COUNT=$(wc -l < "$LOW_VAR_SNPS")
                    echo "  Found $LOW_VAR_COUNT SNPs with low variance in analysis subset"
                    echo "  Adding to exclude list..."
                    cat "$LOW_VAR_SNPS" >> "$EXCLUDE_FILE"
                else
                    echo "  All SNPs in filtered list have sufficient variance in analysis subset"
                fi
            fi
        else
            echo "  WARNING: PLINK check did not produce output file, will rely on retry mechanism"
        fi
    else
        echo "  PLINK not available after module load, skipping pre-scan (will rely on retry mechanism)"
    fi
    else
        echo "  BED/BIM/FAM files not found, skipping pre-scan (will rely on retry mechanism)"
    fi
    
    if [ -s "$EXCLUDE_FILE" ]; then
        EXCLUDE_COUNT=$(wc -l < "$EXCLUDE_FILE")
        echo "  Total SNPs to exclude: $EXCLUDE_COUNT"
    fi
else
    echo "  Skipping pre-scan (no filtered SNP list available for comparison)"
fi

# --- Step 1: null model (test with smaller block size and single chromosome) ---
echo "  Running REGENIE step 1 (test mode - chromosome 22 only)..."
# Retry mechanism to handle SNPs that become monomorphic in analysis subset
# Note: With PLINK pre-scan, retries should be rare. Increased limit for safety.
MAX_RETRIES=20
RETRY_COUNT=0
STEP1_EXIT=1

while [ $STEP1_EXIT -ne 0 ] && [ $RETRY_COUNT -lt $MAX_RETRIES ]; do
    if [ $RETRY_COUNT -gt 0 ]; then
        echo "  Retry attempt $RETRY_COUNT (excluding problematic SNPs)"
    fi
    
    regenie \
      --step 1 \
      --bgen "${GENO_PATH}/ukb_imp_chr22_v3.bgen" \
      --sample "${SAMPLE_FILE}" \
      --phenoFile "${PHENO_FILE}" \
      --covarFile "${COVAR_FILE}" \
      ${EXTRACT_FLAG} \
      $(if [ -s "$EXCLUDE_FILE" ]; then echo "--exclude ${EXCLUDE_FILE}"; fi) \
      --bsize 500 \
      --threads 4 \
      --lowmem \
      --lowmem-prefix "${TMPDIR_STEP1}/regenie_tmp_preds" \
      --force-step1 \
      --out "${OUTDIR_STEP1}/regenie_step1_${PHENO_NAME}" 2>&1 | tee "${OUTDIR_STEP1}/regenie_step1_${PHENO_NAME}.log"
    
    STEP1_EXIT=${PIPESTATUS[0]}
    
    if [ $STEP1_EXIT -ne 0 ]; then
        # Check if failure is due to low variance SNP
        FAILED_SNP=$(grep "has low variance" "${OUTDIR_STEP1}/regenie_step1_${PHENO_NAME}.log" | tail -1 | sed -E 's/.*SNP ([^ ]+) has low variance.*/\1/')
        if [ ! -z "$FAILED_SNP" ] && ! grep -q "^${FAILED_SNP}$" "$EXCLUDE_FILE" 2>/dev/null; then
            echo "$FAILED_SNP" >> "$EXCLUDE_FILE"
            echo "  Adding problematic SNP to exclude list: $FAILED_SNP"
            RETRY_COUNT=$((RETRY_COUNT + 1))
            # Clean up partial output files before retry
            rm -f "${OUTDIR_STEP1}/regenie_step1_${PHENO_NAME}"*.loco "${OUTDIR_STEP1}/regenie_step1_${PHENO_NAME}"*.pred 2>/dev/null
        else
            # Not a low variance issue or already excluded, break retry loop
            break
        fi
    fi
done

# If we hit the retry limit but excluded SNPs, try one final time with all excluded SNPs
if [ $STEP1_EXIT -ne 0 ] && [ $RETRY_COUNT -ge $MAX_RETRIES ] && [ -s "$EXCLUDE_FILE" ]; then
    echo "  Reached retry limit ($MAX_RETRIES), attempting final run with all excluded SNPs..."
    regenie \
      --step 1 \
      --bgen "${GENO_PATH}/ukb_imp_chr22_v3.bgen" \
      --sample "${SAMPLE_FILE}" \
      --phenoFile "${PHENO_FILE}" \
      --covarFile "${COVAR_FILE}" \
      ${EXTRACT_FLAG} \
      $(if [ -s "$EXCLUDE_FILE" ]; then echo "--exclude ${EXCLUDE_FILE}"; fi) \
      --bsize 500 \
      --threads 4 \
      --lowmem \
      --lowmem-prefix "${TMPDIR_STEP1}/regenie_tmp_preds" \
      --force-step1 \
      --out "${OUTDIR_STEP1}/regenie_step1_${PHENO_NAME}" 2>&1 | tee "${OUTDIR_STEP1}/regenie_step1_${PHENO_NAME}.log"
    STEP1_EXIT=${PIPESTATUS[0]}
    
    # If final attempt failed with another low-variance SNP, allow one more retry
    if [ $STEP1_EXIT -ne 0 ]; then
        FINAL_FAILED_SNP=$(grep "has low variance" "${OUTDIR_STEP1}/regenie_step1_${PHENO_NAME}.log" | tail -1 | sed -E 's/.*SNP ([^ ]+) has low variance.*/\1/')
        if [ ! -z "$FINAL_FAILED_SNP" ] && ! grep -q "^${FINAL_FAILED_SNP}$" "$EXCLUDE_FILE" 2>/dev/null; then
            echo "$FINAL_FAILED_SNP" >> "$EXCLUDE_FILE"
            echo "  Final attempt encountered another low-variance SNP: $FINAL_FAILED_SNP"
            echo "  Adding to exclude list and attempting one more time..."
            RETRY_COUNT=$((RETRY_COUNT + 1))
            # Clean up partial output files before retry
            rm -f "${OUTDIR_STEP1}/regenie_step1_${PHENO_NAME}"*.loco "${OUTDIR_STEP1}/regenie_step1_${PHENO_NAME}"*.pred 2>/dev/null
            # One final attempt with all excluded SNPs
            regenie \
              --step 1 \
              --bgen "${GENO_PATH}/ukb_imp_chr22_v3.bgen" \
              --sample "${SAMPLE_FILE}" \
              --phenoFile "${PHENO_FILE}" \
              --covarFile "${COVAR_FILE}" \
              ${EXTRACT_FLAG} \
              $(if [ -s "$EXCLUDE_FILE" ]; then echo "--exclude ${EXCLUDE_FILE}"; fi) \
              --bsize 500 \
              --threads 4 \
              --lowmem \
              --lowmem-prefix "${TMPDIR_STEP1}/regenie_tmp_preds" \
              --force-step1 \
              --out "${OUTDIR_STEP1}/regenie_step1_${PHENO_NAME}" 2>&1 | tee "${OUTDIR_STEP1}/regenie_step1_${PHENO_NAME}.log"
            STEP1_EXIT=${PIPESTATUS[0]}
        fi
    fi
fi

if [ $STEP1_EXIT -ne 0 ]; then
    echo "ERROR: REGENIE step 1 failed with exit code $STEP1_EXIT after $RETRY_COUNT retries"
    if [ -s "$EXCLUDE_FILE" ]; then
        echo "Excluded SNPs:"
        cat "$EXCLUDE_FILE"
    fi
    exit $STEP1_EXIT
fi

if [ $RETRY_COUNT -gt 0 ]; then
    echo "  Step 1 completed after excluding $RETRY_COUNT problematic SNP(s)"
fi

# Check that prediction file was created
PRED_FILE="${OUTDIR_STEP1}/regenie_step1_${PHENO_NAME}_pred.list"
if [ ! -f "$PRED_FILE" ]; then
    echo "ERROR: Prediction file not created: $PRED_FILE"
    exit 1
fi

echo "  Step 1 complete. Prediction file: $PRED_FILE"

# --- Step 2: association testing (test with smaller block size and single chromosome) ---
echo "  Running REGENIE step 2 (test mode - chromosome 22 only)..."
# For step 2, use BED format (converted from BGEN) for consistency and speed
# Using files from ${FAM_PATH} ensures IDs match with Step 1 sample file
CHR_BED="${FAM_PATH}/ukb22.bed"
CHR_BIM="${FAM_PATH}/ukb22.bim"
CHR_FAM="${FAM_PATH}/ukb22.fam"

if [ ! -f "$CHR_BED" ] || [ ! -f "$CHR_BIM" ] || [ ! -f "$CHR_FAM" ]; then
    echo "ERROR: BED/BIM/FAM files not found for chromosome 22"
    echo "  Expected: $CHR_BED, $CHR_BIM, $CHR_FAM"
    exit 1
fi

regenie \
  --step 2 \
  --bed "${FAM_PATH}/ukb22" \
  --phenoFile "${PHENO_FILE}" \
  --covarFile "${COVAR_FILE}" \
  --pred "${PRED_FILE}" \
  ${EXTRACT_FLAG} \
  $(if [ -s "$EXCLUDE_FILE" ]; then echo "--exclude ${EXCLUDE_FILE}"; fi) \
  --bsize 200 \
  --threads 4 \
  --out "${OUTDIR_STEP2}/regenie_step2_${PHENO_NAME}" 2>&1 | tee "${OUTDIR_STEP2}/regenie_step2_${PHENO_NAME}.log"

STEP2_EXIT=${PIPESTATUS[0]}
if [ $STEP2_EXIT -ne 0 ]; then
    echo "ERROR: REGENIE step 2 failed with exit code $STEP2_EXIT"
    exit $STEP2_EXIT
fi

# Check output file
OUTPUT_FILE="${OUTDIR_STEP2}/regenie_step2_${PHENO_NAME}.regenie.gz"
if [ -f "$OUTPUT_FILE" ]; then
    echo "  Step 2 complete. Output file: $OUTPUT_FILE"
    echo "  File size: $(du -h "$OUTPUT_FILE" | cut -f1)"
    echo "  Number of variants: $(zcat "$OUTPUT_FILE" | wc -l)"
else
    echo "WARNING: Expected output file not found: $OUTPUT_FILE"
fi

echo "  Test complete for $PHENO_NAME"
