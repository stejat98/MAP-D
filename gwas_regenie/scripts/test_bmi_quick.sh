#!/bin/bash
# Quick test of REGENIE pipeline with BMI phenotype
# Tests with chromosome 22 only and small block sizes for fast execution

set -e

echo "=========================================="
echo "REGENIE Quick Test - BMI"
echo "=========================================="

# Configuration
STRATUM="full"
SAMPLE_TYPE="hallmarks_heldout"
PHENO_NAME="BMI"
CHR_TEST="22"  # Use chromosome 22 (smallest, fastest for testing)

# Paths
SCRIPT_DIR="/path/to/project/regenie_pipeline/scripts"
INPUT_PATH="/path/to/project/regenie_pipeline/inputs/${STRATUM}/${SAMPLE_TYPE}/${PHENO_NAME}"
OUTPUT_ROOT="/path/to/scratch/regenie_test"
CONDA_ENV="/path/to/project/.conda/envs/regenie_env"
GWAS_PATH="/path/to/shared_data/UKB/gwas"

# Output directories
OUTDIR_STEP1="${OUTPUT_ROOT}/${PHENO_NAME}/${STRATUM}/step1_test"
OUTDIR_STEP2="${OUTPUT_ROOT}/${PHENO_NAME}/${STRATUM}/step2_test"
TMPDIR_STEP1="${OUTDIR_STEP1}/tmpdir_${SLURM_JOB_ID:-$$}"

# Create output directories
mkdir -p "${OUTDIR_STEP1}" "${OUTDIR_STEP2}" "${TMPDIR_STEP1}"

# Check input files
PHENO_FILE="${INPUT_PATH}/pheno.txt"
COVAR_FILE="${INPUT_PATH}/covar.txt"
PSAM_FILE="${GWAS_PATH}/ukb_nonimputed_snps.psam"

if [ ! -f "$PHENO_FILE" ]; then
    echo "ERROR: Phenotype file not found: $PHENO_FILE"
    exit 1
fi

if [ ! -f "$COVAR_FILE" ]; then
    echo "ERROR: Covariate file not found: $COVAR_FILE"
    exit 1
fi

if [ ! -f "$PSAM_FILE" ]; then
    echo "ERROR: PSAM file not found: $PSAM_FILE"
    exit 1
fi

# Filter phenotype/covariate files to match genotype IIDs
# REGENIE requires matching IIDs - if there's no overlap, it will fail
echo "=========================================="
echo "Checking IID overlap between files"
echo "=========================================="

# Extract IIDs from genotype file
GENO_IIDS="${TMPDIR_STEP1}/geno_iids.txt"
awk -F'\t' 'NR>1 && $2>0 {print $2}' "$PSAM_FILE" | sort -n > "$GENO_IIDS"
N_GENO=$(wc -l < "$GENO_IIDS")
echo "  Genotype file has $N_GENO individuals"

# Extract IIDs from phenotype file  
PHENO_IIDS="${TMPDIR_STEP1}/pheno_iids.txt"
tail -n +2 "$PHENO_FILE" | awk '{print $2}' | sort -n > "$PHENO_IIDS"
N_PHENO=$(wc -l < "$PHENO_IIDS")
echo "  Phenotype file has $N_PHENO individuals"

# Check overlap
OVERLAP=$(comm -12 "$PHENO_IIDS" "$GENO_IIDS" | wc -l)
echo "  Overlap: $OVERLAP individuals"

if [ $OVERLAP -eq 0 ]; then
    echo ""
    echo "ERROR: Zero overlap between phenotype and genotype IIDs!"
    echo "  This means the phenotype files were generated from a different dataset."
    echo "  You need to regenerate the phenotype files to include only individuals"
    echo "  present in the genotype files, or use the correct genotype files."
    echo ""
    echo "  To fix: Modify generate_inputs.R to filter main_df to only include"
    echo "  eids that exist in the genotype psam files."
    exit 1
fi

echo "  Using $OVERLAP individuals for analysis"
echo ""

# Count individuals
N_INDIV=$(tail -n +2 "$PHENO_FILE" | wc -l)
echo "Phenotype: $PHENO_NAME"
echo "Stratum: $STRATUM"
echo "Number of individuals (after filtering): $N_INDIV"
echo "Testing chromosome: $CHR_TEST"
echo ""

# Conda setup (IMPORTANT: Set environment variables before activating conda)
export HOME=/path/to/project
export CONDA_NO_PLUGINS=true
export CONDA_PKGS_DIRS=/path/to/project/.conda/pkgs
export CONDA_ENVS_PATH=/path/to/project/.conda/envs

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
echo "  Output: ${OUTDIR_STEP1}/regenie_step1_${PHENO_NAME}"
echo ""

regenie \
  --step 1 \
  --pgen "${GWAS_PATH}/ukb_nonimputed_snps" \
  --phenoFile "${PHENO_FILE}" \
  --covarFile "${COVAR_FILE}" \
  --chr ${CHR_TEST} \
  --bsize 500 \
  --threads 4 \
  --lowmem \
  --lowmem-prefix "${TMPDIR_STEP1}/regenie_tmp_preds" \
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
echo "  Pred: ${PRED_FILE}"
echo "  Output: ${OUTDIR_STEP2}/regenie_step2_${PHENO_NAME}"
echo ""

regenie \
  --step 2 \
  --pgen "${GWAS_PATH}/UKBallchr" \
  --phenoFile "${PHENO_FILE}" \
  --covarFile "${COVAR_FILE}" \
  --pred "${PRED_FILE}" \
  --chr ${CHR_TEST} \
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
