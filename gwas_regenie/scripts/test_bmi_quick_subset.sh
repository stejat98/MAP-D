#!/bin/bash
# Quick test of REGENIE pipeline with BMI phenotype - SMALL SUBSET TEST
# Uses bed/bim/fam files from /n/no_backup2/patel/uk_biobank/ukb_genetics/22881/bgen_converted
# Tests with chromosome 22 only and small block sizes for fast execution

set -e

echo "=========================================="
echo "REGENIE Quick Test - BMI (Small Subset)"
echo "=========================================="

# Configuration
STRATUM="full"
SAMPLE_TYPE="hallmarks_heldout"
PHENO_NAME="BMI"
CHR_TEST="22"  # Use chromosome 22 (smallest, fastest for testing)

# Paths
SCRIPT_DIR="/n/groups/patel/sivateja/regenie_pipeline/scripts"
INPUT_PATH="/n/groups/patel/sivateja/regenie_pipeline/inputs/${STRATUM}/${SAMPLE_TYPE}/${PHENO_NAME}"
OUTPUT_ROOT="/n/scratch/users/s/st320/regenie_test"
CONDA_ENV="/n/groups/patel/sivateja/.conda/envs/regenie_env"
GENO_PATH="/n/no_backup2/patel/uk_biobank/ukb_genetics/22881"
FAM_FILE="${GENO_PATH}/bgen_converted/ukb${CHR_TEST}.fam"

# Output directories
OUTDIR_STEP1="${OUTPUT_ROOT}/${PHENO_NAME}/${STRATUM}/step1_test_subset"
OUTDIR_STEP2="${OUTPUT_ROOT}/${PHENO_NAME}/${STRATUM}/step2_test_subset"
TMPDIR_STEP1="${OUTDIR_STEP1}/tmpdir_${SLURM_JOB_ID:-$$}"

# Create output directories
mkdir -p "${OUTDIR_STEP1}" "${OUTDIR_STEP2}" "${TMPDIR_STEP1}"

# Create sample file from fam file (REGENIE requires ID_1 ID_2 header)
SAMPLE_FILE="${TMPDIR_STEP1}/ukb${CHR_TEST}.sample"

# Check input files
PHENO_FILE="${INPUT_PATH}/pheno.txt"
COVAR_FILE="${INPUT_PATH}/covar.txt"
BGEN_FILE="${GENO_PATH}/ukb_imp_chr${CHR_TEST}_v3.bgen"

if [ ! -f "$PHENO_FILE" ]; then
    echo "ERROR: Phenotype file not found: $PHENO_FILE"
    exit 1
fi

if [ ! -f "$COVAR_FILE" ]; then
    echo "ERROR: Covariate file not found: $COVAR_FILE"
    exit 1
fi

if [ ! -f "$BGEN_FILE" ]; then
    echo "ERROR: BGEN file not found: $BGEN_FILE"
    exit 1
fi

if [ ! -f "$FAM_FILE" ]; then
    echo "ERROR: FAM file not found: $FAM_FILE"
    exit 1
fi

# Create sample file from fam file (REGENIE requires specific format)
echo "Creating sample file from FAM file..."
awk 'BEGIN {print "ID_1 ID_2 missing sex"; print "0 0 0 D"} {print $1, $2, "0", $5}' "$FAM_FILE" > "$SAMPLE_FILE"
echo "  Created: $SAMPLE_FILE ($(wc -l < "$SAMPLE_FILE") lines)"

# Count individuals
N_INDIV=$(tail -n +2 "$PHENO_FILE" | wc -l)
echo "Phenotype: $PHENO_NAME"
echo "Stratum: $STRATUM"
echo "Number of individuals: $N_INDIV"
echo "Testing chromosome: $CHR_TEST"
echo "Genotype file: ${BGEN_FILE}"
echo ""

# Conda setup (IMPORTANT: Set environment variables before activating conda)
export HOME=/n/groups/patel/sivateja
export CONDA_NO_PLUGINS=true
export CONDA_PKGS_DIRS=/n/groups/patel/sivateja/.conda/pkgs
export CONDA_ENVS_PATH=/n/groups/patel/sivateja/.conda/envs

module load conda/miniforge3/24.11.3-0
eval "$(conda shell.bash hook)"
conda activate "${CONDA_ENV}"

# Verify regenie is available
if ! command -v regenie &> /dev/null; then
    echo "ERROR: regenie command not found. Check conda environment."
    exit 1
fi

echo "REGENIE version:"
regenie --version
echo ""

# --- Step 1: null model (test with chromosome 22 only, small block size) ---
echo "=========================================="
echo "Step 1: Running null model (chromosome $CHR_TEST only)"
echo "=========================================="
echo "  Input: ${PHENO_FILE}"
echo "  Genotype: ${BED_FILE}"
echo "  Output: ${OUTDIR_STEP1}/regenie_step1_${PHENO_NAME}"
echo ""

regenie \
  --step 1 \
  --bgen "${BGEN_FILE}" \
  --sample "${SAMPLE_FILE}" \
  --phenoFile "${PHENO_FILE}" \
  --covarFile "${COVAR_FILE}" \
  --bsize 500 \
  --threads 4 \
  --lowmem \
  --lowmem-prefix "${TMPDIR_STEP1}/regenie_tmp_preds" \
  --force-step1 \
  --out "${OUTDIR_STEP1}/regenie_step1_${PHENO_NAME}" 2>&1 | tee "${OUTDIR_STEP1}/regenie_step1_${PHENO_NAME}.log"

STEP1_EXIT=${PIPESTATUS[0]}
if [ $STEP1_EXIT -ne 0 ]; then
    echo ""
    echo "ERROR: REGENIE step 1 failed with exit code $STEP1_EXIT"
    echo "Check log file: ${OUTDIR_STEP1}/regenie_step1_${PHENO_NAME}.log"
    exit $STEP1_EXIT
fi

# Check that prediction file was created
PRED_FILE="${OUTDIR_STEP1}/regenie_step1_${PHENO_NAME}_pred.list"
if [ ! -f "$PRED_FILE" ]; then
    echo "ERROR: Prediction file not created: $PRED_FILE"
    exit 1
fi

echo ""
echo "✓ Step 1 complete!"
echo "  Prediction file: $PRED_FILE"
echo "  File size: $(du -h "$PRED_FILE" | cut -f1)"
echo ""

# --- Step 2: association testing (test with chromosome 22 only, small block size) ---
echo "=========================================="
echo "Step 2: Running association testing (chromosome $CHR_TEST only)"
echo "=========================================="
echo "  Input: ${PHENO_FILE}"
echo "  Genotype: ${BED_FILE}"
echo "  Pred: ${PRED_FILE}"
echo "  Output: ${OUTDIR_STEP2}/regenie_step2_${PHENO_NAME}"
echo ""

regenie \
  --step 2 \
  --bgen "${BGEN_FILE}" \
  --sample "${SAMPLE_FILE}" \
  --phenoFile "${PHENO_FILE}" \
  --covarFile "${COVAR_FILE}" \
  --pred "${PRED_FILE}" \
  --bsize 200 \
  --threads 4 \
  --out "${OUTDIR_STEP2}/regenie_step2_${PHENO_NAME}" 2>&1 | tee "${OUTDIR_STEP2}/regenie_step2_${PHENO_NAME}.log"

STEP2_EXIT=${PIPESTATUS[0]}
if [ $STEP2_EXIT -ne 0 ]; then
    echo ""
    echo "ERROR: REGENIE step 2 failed with exit code $STEP2_EXIT"
    echo "Check log file: ${OUTDIR_STEP2}/regenie_step2_${PHENO_NAME}.log"
    exit $STEP2_EXIT
fi

# Check output file
OUTPUT_FILE="${OUTDIR_STEP2}/regenie_step2_${PHENO_NAME}.regenie.gz"
if [ -f "$OUTPUT_FILE" ]; then
    echo ""
    echo "✓ Step 2 complete!"
    echo "  Output file: $OUTPUT_FILE"
    echo "  File size: $(du -h "$OUTPUT_FILE" | cut -f1)"
    N_VARIANTS=$(zcat "$OUTPUT_FILE" | tail -n +2 | wc -l)
    echo "  Number of variants: $N_VARIANTS"
    echo ""
    echo "=========================================="
    echo "Test completed successfully!"
    echo "=========================================="
    echo "Outputs:"
    echo "  Step 1: ${OUTDIR_STEP1}"
    echo "  Step 2: ${OUTDIR_STEP2}"
    echo ""
    echo "To view results:"
    echo "  zcat ${OUTPUT_FILE} | head -20"
else
    echo ""
    echo "WARNING: Expected output file not found: $OUTPUT_FILE"
    exit 1
fi
