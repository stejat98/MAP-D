#!/bin/bash
#SBATCH -c 8
#SBATCH -t 20-00:00:00
#SBATCH --mem=50G
#SBATCH -p long
# Note: Time limit should be overridden by submit scripts based on stratum:
#   - Full cohort: 10-00:00:00 (long partition, ~400K individuals)
#   - Full: 20-00:00:00 (long partition, ~418K individuals)
# These limits account for sample size scaling and fairshare optimization
#SBATCH -o /n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_%A_%a.out
#SBATCH -e /n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sivateja_tangirala@hms.harvard.edu

# REGENIE GWAS array job
# Runs Step 1 (null model) and Step 2 (association testing) for one phenotype

# Parse command line arguments or use SLURM_ARRAY_TASK_ID
if [ -z "$1" ]; then
    # Get phenotype info from array task ID
    PHENO_LIST_FILE="${PHENO_LIST_FILE:-/n/groups/patel/sivateja/regenie_pipeline/scripts/pheno_list_hallmarks.txt}"
    if [ ! -f "$PHENO_LIST_FILE" ]; then
        echo "ERROR: Phenotype list file not found: $PHENO_LIST_FILE"
        exit 1
    fi
    PHENO_LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$PHENO_LIST_FILE")
    STRATUM=$(echo "$PHENO_LINE" | cut -d'|' -f1)
    SAMPLE_TYPE=$(echo "$PHENO_LINE" | cut -d'|' -f2)
    PHENO_NAME=$(echo "$PHENO_LINE" | cut -d'|' -f3)
else
    # Direct specification (for testing)
    STRATUM="$1"
    SAMPLE_TYPE="$2"
    PHENO_NAME="$3"
fi

echo "=========================================="
echo "REGENIE GWAS"
echo "Stratum: $STRATUM"
echo "Sample type: $SAMPLE_TYPE"
echo "Phenotype: $PHENO_NAME"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID:-N/A}"
echo "=========================================="

# Paths
CONDA_ENV="${CONDA_ENV:-/n/groups/patel/sivateja/.conda/envs/regenie_env}"
# IMPORTANT: Use ONLY genetic data from /n/no_backup2/patel/uk_biobank/ukb_genetics/22881
# to ensure ID matching across all input files (phenotype, genotype, sample files)
# BGEN files: ${GENO_PATH}/ukb_imp_chr*_v3.bgen (for Step 1)
# BED/BIM/FAM files: ${FAM_PATH}/ukb*.{bed,bim,fam} (for Step 2, converted from BGEN)
GENO_PATH="/n/no_backup2/patel/uk_biobank/ukb_genetics/22881"
FAM_PATH="${GENO_PATH}/bgen_converted"
# Optional override: separate input trees (e.g. egfr_pwas/inputs) without moving shared inputs
REGENIE_INPUT_ROOT="${REGENIE_INPUT_ROOT:-/n/groups/patel/sivateja/regenie_pipeline/inputs}"
INPUT_PATH="${REGENIE_INPUT_ROOT}/${STRATUM}/${SAMPLE_TYPE}/${PHENO_NAME}"
OUTPUT_ROOT="${REGENIE_OUTPUT_ROOT:-/n/scratch/users/s/st320/regenie}"

# Per-phenotype output dirs & prefixes (using scratch convention)
OUTDIR_STEP1="${OUTPUT_ROOT}/${PHENO_NAME}/${STRATUM}/step1"
OUTDIR_STEP2="${OUTPUT_ROOT}/${PHENO_NAME}/${STRATUM}/step2"
TMPDIR_STEP1="${OUTDIR_STEP1}/tmpdir_${SLURM_JOB_ID:-$$}"

# Create output directories
mkdir -p "${OUTDIR_STEP1}" "${OUTDIR_STEP2}" "${TMPDIR_STEP1}"

# Check input files exist
PHENO_FILE="${INPUT_PATH}/pheno.txt"
COVAR_FILE="${INPUT_PATH}/covar.txt"

if [ ! -f "$PHENO_FILE" ]; then
    echo "ERROR: Phenotype file not found: $PHENO_FILE"
    exit 1
fi

if [ ! -f "$COVAR_FILE" ]; then
    echo "ERROR: Covariate file not found: $COVAR_FILE"
    exit 1
fi

# Initialize module system (required for SLURM jobs)
if [ -f /etc/profile.d/modules.sh ]; then
    source /etc/profile.d/modules.sh
fi

# Conda setup (IMPORTANT: Set environment variables before activating conda)
export HOME=/n/groups/patel/sivateja
export CONDA_NO_PLUGINS=true
export CONDA_PKGS_DIRS=/n/groups/patel/sivateja/.conda/pkgs
export CONDA_ENVS_PATH=/n/groups/patel/sivateja/.conda/envs

module load conda/miniforge3/24.11.3-0
eval "$(conda shell.bash hook)"
conda activate "${CONDA_ENV}"

# Check if Step 1 is already completed (for resume functionality)
PRED_FILE="${OUTDIR_STEP1}/regenie_step1_${PHENO_NAME}_pred.list"
STEP1_COMPLETE=0

if [ -f "$PRED_FILE" ] && [ -s "$PRED_FILE" ]; then
    echo "[${PHENO_NAME}] Step 1 already completed, skipping..."
    echo "  Using existing prediction file: $PRED_FILE"
    STEP1_COMPLETE=1
fi

# Determine if binary trait (check if phenotype has only 0/1/NA values)
# For now, assume continuous unless specified
BINARY_FLAG=""
if [ -n "$BINARY_TRAIT" ] && [ "$BINARY_TRAIT" = "1" ]; then
    BINARY_FLAG="--bt"
fi

# --- Step 1: null model (skip if already completed) ---
if [ $STEP1_COMPLETE -eq 0 ]; then
    echo "[${PHENO_NAME}] Running REGENIE step 1"
    echo "  Input: ${PHENO_FILE}"
    echo "  Output: ${OUTDIR_STEP1}/regenie_step1_${PHENO_NAME}"

    # Create sample file from FAM file (required for BGEN format, ensures ID matching)
    # The sample file must match the IDs in the BGEN file exactly
    # Using FAM file from ${FAM_PATH} ensures consistency with BED files used in Step 2
    SAMPLE_FILE="${TMPDIR_STEP1}/sample_file.sample"
    awk 'BEGIN {print "ID_1 ID_2 missing sex"; print "0 0 0 D"} {print $1, $2, "0", $5}' "${FAM_PATH}/ukb22.fam" > "$SAMPLE_FILE"
    
    # Use pre-filtered SNP list to avoid low-variance SNPs (per REGENIE recommendations)
    # Use MAF-based filtering for consistency across different sample sizes
    # Priority: MAF0.001 > MAC50 > MAC100 (most conservative first)
    FILTERED_SNP_LIST="/n/groups/patel/sivateja/regenie_pipeline/filtered_snps/chr22_maf0.001.snplist"
    # Fallback to MAC-based filters if MAF doesn't exist
    if [ ! -f "$FILTERED_SNP_LIST" ]; then
        FILTERED_SNP_LIST="/n/groups/patel/sivateja/regenie_pipeline/filtered_snps/chr22_mac50.snplist"
    fi
    if [ ! -f "$FILTERED_SNP_LIST" ]; then
        FILTERED_SNP_LIST="/n/groups/patel/sivateja/regenie_pipeline/filtered_snps/chr22_mac100.snplist"
    fi
    
    # Create exclude file for problematic low-variance SNPs that slip through
    # (can happen when MAC filter is applied to full dataset but analysis uses subset)
    EXCLUDE_FILE="${TMPDIR_STEP1}/exclude_snps.txt"
    touch "$EXCLUDE_FILE"
    
    if [ ! -f "$FILTERED_SNP_LIST" ]; then
        echo "ERROR: Pre-filtered SNP list not found: $FILTERED_SNP_LIST"
        echo "  Run prefilter_snps.sh first to create filtered SNP list"
        echo "  Cannot continue without SNP filtering (will encounter low-variance errors)"
        exit 1
    else
        echo "[${PHENO_NAME}] ✓ Using pre-filtered SNP list: $FILTERED_SNP_LIST"
        SNP_COUNT=$(wc -l < "$FILTERED_SNP_LIST")
        echo "  SNPs in filtered list: $SNP_COUNT"
        echo "  Filter type: $(basename "$FILTERED_SNP_LIST" .snplist | sed 's/chr22_//')"
        EXTRACT_FLAG="--extract ${FILTERED_SNP_LIST}"
    fi
    
    # Pre-scan and exclude low-variance SNPs for this specific analysis subset
    # This is more efficient than retrying REGENIE multiple times
    # Only do this if we have a filtered SNP list to compare against
    if [ ! -z "$EXTRACT_FLAG" ] && [ -f "$FILTERED_SNP_LIST" ]; then
        echo "[${PHENO_NAME}] Pre-scanning for low-variance SNPs in analysis subset..."
        # Extract sample IDs from phenotype file (first column, skip header)
        SAMPLE_ID_FILE="${TMPDIR_STEP1}/analysis_samples.txt"
        tail -n +2 "$PHENO_FILE" | awk '{print $1, $1}' > "$SAMPLE_ID_FILE"
        N_ANALYSIS_SAMPLES=$(wc -l < "$SAMPLE_ID_FILE")
        echo "[${PHENO_NAME}] Analysis subset contains $N_ANALYSIS_SAMPLES samples"
        
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
            echo "[${PHENO_NAME}] Checking MAC on analysis subset using PLINK..."
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
                LOW_VAR_SNPS="${TMPDIR_STEP1}/low_variance_snps.txt"
                # Get SNPs in filtered list but not in PLINK output (comm -23: lines only in first file)
                comm -23 <(sort "$FILTERED_SNP_LIST") <(sort "${PLINK_CHECK_OUT}.snplist") > "$LOW_VAR_SNPS" 2>/dev/null || true
                
                if [ -s "$LOW_VAR_SNPS" ]; then
                    LOW_VAR_COUNT=$(wc -l < "$LOW_VAR_SNPS")
                    echo "[${PHENO_NAME}] Found $LOW_VAR_COUNT SNPs with low variance in analysis subset"
                    echo "[${PHENO_NAME}] Adding to exclude list..."
                    cat "$LOW_VAR_SNPS" >> "$EXCLUDE_FILE"
                else
                    echo "[${PHENO_NAME}] All SNPs in filtered list have sufficient variance in analysis subset"
                fi
            else
                echo "[${PHENO_NAME}] WARNING: PLINK check did not produce output file, will rely on retry mechanism"
            fi
        else
            echo "[${PHENO_NAME}] PLINK not available after module load, skipping pre-scan (will rely on retry mechanism)"
        fi
        else
            echo "[${PHENO_NAME}] BED/BIM/FAM files not found, skipping pre-scan (will rely on retry mechanism)"
        fi
        
        if [ -s "$EXCLUDE_FILE" ]; then
            EXCLUDE_COUNT=$(wc -l < "$EXCLUDE_FILE")
            echo "[${PHENO_NAME}] Total SNPs to exclude: $EXCLUDE_COUNT"
        fi
    else
        echo "[${PHENO_NAME}] Skipping pre-scan (no filtered SNP list available for comparison)"
    fi
    
    # Retry mechanism to handle SNPs that become monomorphic in analysis subset
    MAX_RETRIES=10
    RETRY_COUNT=0
    STEP1_EXIT=1
    
    while [ $STEP1_EXIT -ne 0 ] && [ $RETRY_COUNT -lt $MAX_RETRIES ]; do
        if [ $RETRY_COUNT -gt 0 ]; then
            echo "[${PHENO_NAME}] Retry attempt $RETRY_COUNT (excluding problematic SNPs)"
        fi
        
        regenie \
          --step 1 \
          --bgen "${GENO_PATH}/ukb_imp_chr22_v3.bgen" \
          --sample "${SAMPLE_FILE}" \
          --phenoFile "${PHENO_FILE}" \
          --covarFile "${COVAR_FILE}" \
          ${BINARY_FLAG} \
          ${EXTRACT_FLAG} \
          $(if [ -s "$EXCLUDE_FILE" ]; then echo "--exclude ${EXCLUDE_FILE}"; fi) \
          --bsize 1000 \
          --threads 8 \
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
                echo "[${PHENO_NAME}] Adding problematic SNP to exclude list: $FAILED_SNP"
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
        echo "[${PHENO_NAME}] Reached retry limit ($MAX_RETRIES), attempting final run with all excluded SNPs..."
        regenie \
          --step 1 \
          --bgen "${GENO_PATH}/ukb_imp_chr22_v3.bgen" \
          --sample "${SAMPLE_FILE}" \
          --phenoFile "${PHENO_FILE}" \
          --covarFile "${COVAR_FILE}" \
          ${BINARY_FLAG} \
          ${EXTRACT_FLAG} \
          $(if [ -s "$EXCLUDE_FILE" ]; then echo "--exclude ${EXCLUDE_FILE}"; fi) \
          --bsize 1000 \
          --threads 8 \
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
                echo "[${PHENO_NAME}] Final attempt encountered another low-variance SNP: $FINAL_FAILED_SNP"
                echo "[${PHENO_NAME}] Adding to exclude list and attempting one more time..."
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
                  ${BINARY_FLAG} \
                  ${EXTRACT_FLAG} \
                  $(if [ -s "$EXCLUDE_FILE" ]; then echo "--exclude ${EXCLUDE_FILE}"; fi) \
                  --bsize 1000 \
                  --threads 8 \
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
        echo "[${PHENO_NAME}] Step 1 completed after excluding $RETRY_COUNT problematic SNP(s)"
    fi

    # Check that prediction file was created
    if [ ! -f "$PRED_FILE" ]; then
        echo "ERROR: Prediction file not created: $PRED_FILE"
        exit 1
    fi
fi

# --- Step 2: association testing (all chromosomes) ---
echo "[${PHENO_NAME}] Running REGENIE step 2 for all chromosomes"
echo "  Input: ${PHENO_FILE}"
echo "  Pred: ${PRED_FILE}"
echo "  Output: ${OUTDIR_STEP2}/regenie_step2_${PHENO_NAME}"

STEP2_FLAGS=""
if [ -n "$BINARY_FLAG" ]; then
    # For binary traits, use Firth correction
    STEP2_FLAGS="--bt --firth --approx --pThresh 0.01"
fi

# Step 2: Loop through all chromosomes (1-22)
# Use BED format for step 2 (faster than BGEN)
# Using files from ${FAM_PATH} ensures IDs match with Step 1 sample file
# REGENIE_STEP2_RESUME=1: skip chromosomes that already have a non-empty .regenie / .regenie.gz
STEP2_ANY_FAIL=0
for CHR in {1..22}; do
    CHR_BED="${FAM_PATH}/ukb${CHR}.bed"
    CHR_BIM="${FAM_PATH}/ukb${CHR}.bim"
    CHR_FAM="${FAM_PATH}/ukb${CHR}.fam"
    
    if [ ! -f "$CHR_BED" ] || [ ! -f "$CHR_BIM" ] || [ ! -f "$CHR_FAM" ]; then
        echo "WARNING: Chromosome ${CHR} files not found, skipping..."
        continue
    fi

    S2_REG="${OUTDIR_STEP2}/regenie_step2_${PHENO_NAME}_chr${CHR}_${PHENO_NAME}.regenie"
    S2_REG_GZ="${S2_REG}.gz"
    if [ "${REGENIE_STEP2_RESUME:-0}" = "1" ] && { [ -s "$S2_REG" ] || [ -s "$S2_REG_GZ" ]; }; then
        echo ""
        echo "[${PHENO_NAME}] Chromosome ${CHR} output already present, skipping (REGENIE_STEP2_RESUME=1)"
        continue
    fi
    
    echo ""
    echo "[${PHENO_NAME}] Processing chromosome ${CHR}..."
    
    regenie \
      --step 2 \
      --bed "${FAM_PATH}/ukb${CHR}" \
      --phenoFile "${PHENO_FILE}" \
      --covarFile "${COVAR_FILE}" \
      --pred "${PRED_FILE}" \
      ${STEP2_FLAGS} \
      --bsize 400 \
      --threads 8 \
      --out "${OUTDIR_STEP2}/regenie_step2_${PHENO_NAME}_chr${CHR}" 2>&1 | tee "${OUTDIR_STEP2}/regenie_step2_${PHENO_NAME}_chr${CHR}.log"
    
    CHR_EXIT=${PIPESTATUS[0]}
    if [ $CHR_EXIT -ne 0 ]; then
        echo "ERROR: REGENIE step 2 failed for chromosome ${CHR} with exit code $CHR_EXIT"
        STEP2_ANY_FAIL=1
        # Continue with other chromosomes instead of exiting
    else
        echo "[${PHENO_NAME}] Chromosome ${CHR} complete"
    fi
done

# Combine chromosome results (optional - REGENIE outputs separate files per chromosome)
echo ""
echo "[${PHENO_NAME}] Step 2 processing complete for all chromosomes"
echo "  Output files: ${OUTDIR_STEP2}/regenie_step2_${PHENO_NAME}_chr*.regenie.gz"

# Do not use PIPESTATUS here: after the loop the last command was often a successful
# echo, not the regenie pipeline — that mis-reported step 2 and could yield bogus exits.
if [ $STEP2_ANY_FAIL -ne 0 ]; then
    echo "ERROR: REGENIE step 2 had one or more chromosome failures (see log lines above)"
    exit 1
fi

echo "[${PHENO_NAME}] REGENIE GWAS complete"
echo "  Step 1 output: ${OUTDIR_STEP1}"
echo "  Step 2 output: ${OUTDIR_STEP2}"
