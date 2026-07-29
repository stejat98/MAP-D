#!/bin/bash
# Test REGENIE pipeline on a small subset
# Tests one protein and one hallmark with subset of individuals

set -e

echo "=========================================="
echo "REGENIE Pipeline Test"
echo "=========================================="

SCRIPT_DIR="/n/groups/patel/sivateja/regenie_pipeline/scripts"
INPUT_BASE="/n/groups/patel/sivateja/regenie_pipeline/inputs"
TEST_OUTPUT_BASE="/n/scratch/users/s/st320/regenie_test"

# Create test output directory
mkdir -p "$TEST_OUTPUT_BASE"

# Find one protein and one hallmark to test
STRATUM="full"
PROTEINS_DIR="${INPUT_BASE}/${STRATUM}/proteins_only"
HALLMARKS_DIR="${INPUT_BASE}/${STRATUM}/hallmarks_heldout"

# Get first available protein
PROTEIN_NAME=$(ls "$PROTEINS_DIR" 2>/dev/null | head -1)
if [ -z "$PROTEIN_NAME" ]; then
    echo "ERROR: No proteins found. Run generate_inputs.R first."
    exit 1
fi

# Get first available hallmark
HALLMARK_NAME=$(ls "$HALLMARKS_DIR" 2>/dev/null | head -1)
if [ -z "$HALLMARK_NAME" ]; then
    echo "ERROR: No hallmarks found. Run generate_inputs.R first."
    exit 1
fi

echo "Test phenotype 1 (protein): $PROTEIN_NAME"
echo "Test phenotype 2 (hallmark): $HALLMARK_NAME"
echo ""

# Create subset of individuals for testing (first 1000)
create_test_subset() {
    local INPUT_FILE="$1"
    local OUTPUT_FILE="$2"
    local SUBSET_SIZE="${3:-1000}"
    
    # Create subset (first N individuals)
    head -n $((SUBSET_SIZE + 1)) "$INPUT_FILE" > "$OUTPUT_FILE"
    echo "Created test subset: $OUTPUT_FILE ($SUBSET_SIZE individuals)"
}

# Test protein
echo "=== Testing Protein: $PROTEIN_NAME ==="
PROTEIN_INPUT="${PROTEINS_DIR}/${PROTEIN_NAME}"
PROTEIN_TEST_DIR="${TEST_OUTPUT_BASE}/${PROTEIN_NAME}/${STRATUM}/input"
mkdir -p "$PROTEIN_TEST_DIR"

# Create test subsets
create_test_subset "${PROTEIN_INPUT}/pheno.txt" "${PROTEIN_TEST_DIR}/pheno.txt" 1000
create_test_subset "${PROTEIN_INPUT}/covar.txt" "${PROTEIN_TEST_DIR}/covar.txt" 1000

# Run REGENIE test for protein (using scratch directory convention)
echo "Running REGENIE test for protein..."
bash "${SCRIPT_DIR}/run_regenie_test.sh" "$STRATUM" "proteins_only" "$PROTEIN_NAME" "$PROTEIN_TEST_DIR" "${TEST_OUTPUT_BASE}/${PROTEIN_NAME}/${STRATUM}/step1" "${TEST_OUTPUT_BASE}/${PROTEIN_NAME}/${STRATUM}/step2"

# Test hallmark
echo ""
echo "=== Testing Hallmark: $HALLMARK_NAME ==="
HALLMARK_INPUT="${HALLMARKS_DIR}/${HALLMARK_NAME}"
HALLMARK_TEST_DIR="${TEST_OUTPUT_BASE}/${HALLMARK_NAME}/${STRATUM}/input"
mkdir -p "$HALLMARK_TEST_DIR"

# Create test subsets
create_test_subset "${HALLMARK_INPUT}/pheno.txt" "${HALLMARK_TEST_DIR}/pheno.txt" 1000
create_test_subset "${HALLMARK_INPUT}/covar.txt" "${HALLMARK_TEST_DIR}/covar.txt" 1000

# Run REGENIE test for hallmark (using scratch directory convention)
echo "Running REGENIE test for hallmark..."
bash "${SCRIPT_DIR}/run_regenie_test.sh" "$STRATUM" "hallmarks_heldout" "$HALLMARK_NAME" "$HALLMARK_TEST_DIR" "${TEST_OUTPUT_BASE}/${HALLMARK_NAME}/${STRATUM}/step1" "${TEST_OUTPUT_BASE}/${HALLMARK_NAME}/${STRATUM}/step2"

echo ""
echo "=========================================="
echo "Test complete! Check outputs in:"
echo "  $TEST_OUTPUT_BASE"
echo "=========================================="
