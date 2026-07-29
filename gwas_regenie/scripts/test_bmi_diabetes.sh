#!/bin/bash
# Test REGENIE pipeline with BMI for diabetes stratum
# Quick test before running all hallmarks across all strata

set -e

echo "=========================================="
echo "REGENIE Pipeline Test: BMI (Diabetes Stratum)"
echo "=========================================="

SCRIPT_DIR="/n/groups/patel/sivateja/regenie_pipeline/scripts"
INPUT_BASE="/n/groups/patel/sivateja/regenie_pipeline/inputs"
TEST_OUTPUT_BASE="/n/scratch/users/s/st320/regenie_test"

# Test parameters
STRATUM="diabetes"
SAMPLE_TYPE="hallmarks_heldout"
PHENO_NAME="BMI"

# Create test output directory
mkdir -p "$TEST_OUTPUT_BASE"

# Check input files exist
PHENO_INPUT="${INPUT_BASE}/${STRATUM}/${SAMPLE_TYPE}/${PHENO_NAME}"
PHENO_FILE="${PHENO_INPUT}/pheno.txt"
COVAR_FILE="${PHENO_INPUT}/covar.txt"

if [ ! -f "$PHENO_FILE" ]; then
    echo "ERROR: Phenotype file not found: $PHENO_FILE"
    exit 1
fi

if [ ! -f "$COVAR_FILE" ]; then
    echo "ERROR: Covariate file not found: $COVAR_FILE"
    exit 1
fi

echo "Stratum: $STRATUM"
echo "Sample type: $SAMPLE_TYPE"
echo "Phenotype: $PHENO_NAME"
echo "Input directory: $PHENO_INPUT"
echo ""

# Count individuals
N_INDIV=$(tail -n +2 "$PHENO_FILE" | wc -l)
echo "Number of individuals: $N_INDIV"
echo ""

# Output directories
OUTDIR_STEP1="${TEST_OUTPUT_BASE}/${PHENO_NAME}/${STRATUM}/step1"
OUTDIR_STEP2="${TEST_OUTPUT_BASE}/${PHENO_NAME}/${STRATUM}/step2"

echo "Step 1 output: $OUTDIR_STEP1"
echo "Step 2 output: $OUTDIR_STEP2"
echo ""

# Run REGENIE test
echo "=== Running REGENIE test ==="
bash "${SCRIPT_DIR}/run_regenie_test.sh" \
    "$STRATUM" \
    "$SAMPLE_TYPE" \
    "$PHENO_NAME" \
    "$PHENO_INPUT" \
    "$OUTDIR_STEP1" \
    "$OUTDIR_STEP2"

echo ""
echo "=========================================="
echo "Test complete! Check outputs in:"
echo "  Step 1: $OUTDIR_STEP1"
echo "  Step 2: $OUTDIR_STEP2"
echo "=========================================="
