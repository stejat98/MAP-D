#!/bin/bash
# MINIMAL REGENIE test - runs in ~2-3 minutes
# Uses only 1000 individuals and ~1000 SNPs for quick testing

set -e

echo "=========================================="
echo "REGENIE Minimal Test - BMI (Tiny Subset)"
echo "=========================================="

# Configuration
STRATUM="full"
SAMPLE_TYPE="hallmarks_heldout"
PHENO_NAME="BMI"
CHR_TEST="22"

# Paths
INPUT_PATH="/n/groups/patel/sivateja/regenie_pipeline/inputs/${STRATUM}/${SAMPLE_TYPE}/${PHENO_NAME}"
OUTPUT_ROOT="/n/scratch/users/s/st320/regenie_test"
CONDA_ENV="/n/groups/patel/sivateja/.conda/envs/regenie_env"
GENO_PATH="/n/no_backup2/patel/uk_biobank/ukb_genetics/22881"
FAM_FILE="${GENO_PATH}/bgen_converted/ukb${CHR_TEST}.fam"

# Output directories
OUTDIR_STEP1="${OUTPUT_ROOT}/${PHENO_NAME}/${STRATUM}/step1_minimal"
OUTDIR_STEP2="${OUTPUT_ROOT}/${PHENO_NAME}/${STRATUM}/step2_minimal"
TMPDIR_STEP1="${OUTDIR_STEP1}/tmpdir_${SLURM_JOB_ID:-$$}"

# Create output directories
mkdir -p "${OUTDIR_STEP1}" "${OUTDIR_STEP2}" "${TMPDIR_STEP1}"

# Check input files
PHENO_FILE="${INPUT_PATH}/pheno.txt"
COVAR_FILE="${INPUT_PATH}/covar.txt"
BGEN_FILE="${GENO_PATH}/ukb_imp_chr${CHR_TEST}_v3.bgen"

if [ ! -f "$PHENO_FILE" ] || [ ! -f "$COVAR_FILE" ] || [ ! -f "$BGEN_FILE" ] || [ ! -f "$FAM_FILE" ]; then
    echo "ERROR: Required input files not found"
    exit 1
fi

# Create MINIMAL filtered files (5000 individuals for better SNP coverage)
echo "Creating minimal test files (5000 individuals)..."
MIN_PHENO="${TMPDIR_STEP1}/pheno_minimal.txt"
MIN_COVAR="${TMPDIR_STEP1}/covar_minimal.txt"

# Get first 5000 IIDs from fam file
FIRST_5000_IIDS="${TMPDIR_STEP1}/first_5000_iids.txt"
head -5001 "$FAM_FILE" | tail -n +2 | awk '{print $2}' | sort -n > "$FIRST_5000_IIDS"

# Filter phenotype file
head -1 "$PHENO_FILE" > "$MIN_PHENO"
awk -v iids_file="$FIRST_5000_IIDS" '
    BEGIN {
        while ((getline line < iids_file) > 0) {
            valid[line] = 1
        }
        close(iids_file)
    }
    NR > 1 && $2 in valid {
        print
        if (++count >= 5000) exit
    }
' "$PHENO_FILE" >> "$MIN_PHENO"

# Filter covariate file
head -1 "$COVAR_FILE" > "$MIN_COVAR"
awk -v iids_file="$FIRST_5000_IIDS" '
    BEGIN {
        while ((getline line < iids_file) > 0) {
            valid[line] = 1
        }
        close(iids_file)
    }
    NR > 1 && $2 in valid {
        print
        if (++count >= 5000) exit
    }
' "$COVAR_FILE" >> "$MIN_COVAR"

N_FILTERED=$(tail -n +2 "$MIN_PHENO" | wc -l)
echo "  Created minimal files with $N_FILTERED individuals"

# Create sample file (must match BGEN file exactly - use full file)
SAMPLE_FILE="${TMPDIR_STEP1}/ukb${CHR_TEST}.sample"
awk 'BEGIN {print "ID_1 ID_2 missing sex"; print "0 0 0 D"} {print $1, $2, "0", $5}' "$FAM_FILE" > "$SAMPLE_FILE"
echo "  Created sample file (full, matches BGEN)"

# Conda setup
export HOME=/n/groups/patel/sivateja
export CONDA_NO_PLUGINS=true
export CONDA_PKGS_DIRS=/n/groups/patel/sivateja/.conda/pkgs
export CONDA_ENVS_PATH=/n/groups/patel/sivateja/.conda/envs

module load conda/miniforge3/24.11.3-0
eval "$(conda shell.bash hook)"
conda activate "${CONDA_ENV}"

if ! command -v regenie &> /dev/null; then
    echo "ERROR: regenie command not found"
    exit 1
fi

echo ""
echo "REGENIE version: $(regenie --version)"
echo ""

# --- Step 1: null model (minimal settings) ---
echo "=========================================="
echo "Step 1: Running null model (minimal test)"
echo "=========================================="
echo "  Individuals: $N_FILTERED"
echo "  Block size: 100 (small for speed)"
echo "  Min MAC: 5 (filters low-variance SNPs)"
echo ""

# Create exclude file for problematic low-variance SNPs (empty for now, REGENIE will handle)
EXCLUDE_FILE="${TMPDIR_STEP1}/exclude_snps.txt"
touch "$EXCLUDE_FILE"

regenie \
  --step 1 \
  --bgen "${BGEN_FILE}" \
  --sample "${SAMPLE_FILE}" \
  --phenoFile "${MIN_PHENO}" \
  --covarFile "${MIN_COVAR}" \
  --bsize 200 \
  --exclude "${EXCLUDE_FILE}" \
  --threads 2 \
  --lowmem \
  --lowmem-prefix "${TMPDIR_STEP1}/regenie_tmp_preds" \
  --force-step1 \
  --out "${OUTDIR_STEP1}/regenie_step1_${PHENO_NAME}" 2>&1 | tee "${OUTDIR_STEP1}/regenie_step1_${PHENO_NAME}.log" || {
    # If it fails on low variance, extract the SNP and retry
    FAILED_SNP=$(grep "has low variance" "${OUTDIR_STEP1}/regenie_step1_${PHENO_NAME}.log" | tail -1 | awk '{print $5}')
    if [ ! -z "$FAILED_SNP" ]; then
        echo "$FAILED_SNP" >> "$EXCLUDE_FILE"
        echo "Retrying with excluded SNP: $FAILED_SNP"
        regenie \
          --step 1 \
          --bgen "${BGEN_FILE}" \
          --sample "${SAMPLE_FILE}" \
          --phenoFile "${MIN_PHENO}" \
          --covarFile "${MIN_COVAR}" \
          --bsize 200 \
          --exclude "${EXCLUDE_FILE}" \
          --threads 2 \
          --lowmem \
          --lowmem-prefix "${TMPDIR_STEP1}/regenie_tmp_preds" \
          --force-step1 \
          --out "${OUTDIR_STEP1}/regenie_step1_${PHENO_NAME}" 2>&1 | tee -a "${OUTDIR_STEP1}/regenie_step1_${PHENO_NAME}.log"
    else
        exit 1
    fi
}

STEP1_EXIT=${PIPESTATUS[0]}
if [ $STEP1_EXIT -ne 0 ]; then
    echo ""
    echo "ERROR: REGENIE step 1 failed with exit code $STEP1_EXIT"
    tail -20 "${OUTDIR_STEP1}/regenie_step1_${PHENO_NAME}.log"
    exit $STEP1_EXIT
fi

# Check prediction file
PRED_FILE="${OUTDIR_STEP1}/regenie_step1_${PHENO_NAME}_pred.list"
if [ ! -f "$PRED_FILE" ]; then
    echo "ERROR: Prediction file not created"
    exit 1
fi

echo ""
echo "✓ Step 1 complete!"
echo "  Prediction file: $PRED_FILE"
echo ""

# --- Step 2: association testing (minimal) ---
echo "=========================================="
echo "Step 2: Running association testing (minimal)"
echo "=========================================="
echo "  Using same minimal files"
echo ""

regenie \
  --step 2 \
  --bgen "${BGEN_FILE}" \
  --sample "${SAMPLE_FILE}" \
  --phenoFile "${MIN_PHENO}" \
  --covarFile "${MIN_COVAR}" \
  --pred "${PRED_FILE}" \
  --bsize 50 \
  --minMAC 1 \
  --threads 2 \
  --out "${OUTDIR_STEP2}/regenie_step2_${PHENO_NAME}" 2>&1 | tee "${OUTDIR_STEP2}/regenie_step2_${PHENO_NAME}.log"

STEP2_EXIT=${PIPESTATUS[0]}
if [ $STEP2_EXIT -ne 0 ]; then
    echo ""
    echo "ERROR: REGENIE step 2 failed with exit code $STEP2_EXIT"
    tail -20 "${OUTDIR_STEP2}/regenie_step2_${PHENO_NAME}.log"
    exit $STEP2_EXIT
fi

# Check output
OUTPUT_FILE="${OUTDIR_STEP2}/regenie_step2_${PHENO_NAME}.regenie.gz"
if [ -f "$OUTPUT_FILE" ]; then
    echo ""
    echo "✓ Step 2 complete!"
    echo "  Output file: $OUTPUT_FILE"
    N_VARIANTS=$(zcat "$OUTPUT_FILE" | tail -n +2 | wc -l)
    echo "  Number of variants: $N_VARIANTS"
    echo ""
    echo "=========================================="
    echo "✓✓✓ Minimal test completed successfully! ✓✓✓"
    echo "=========================================="
    echo ""
    echo "First 5 results:"
    zcat "$OUTPUT_FILE" | head -6
else
    echo ""
    echo "WARNING: Expected output file not found: $OUTPUT_FILE"
    exit 1
fi
