#!/bin/bash
#SBATCH -c 8
#SBATCH -t 4:00:00
#SBATCH --mem=50G
#SBATCH -p medium
#SBATCH -o /path/to/project/regenie_pipeline/slurm/test_bmi_filter_%j.out
#SBATCH -e /path/to/project/regenie_pipeline/slurm/test_bmi_filter_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@example.com
#SBATCH --job-name=regenie_test_bmi

# Test REGENIE with BMI using the main pipeline script
# This ensures we're using the same script and filtered SNP list as the full pipeline

echo "=========================================="
echo "Testing REGENIE with BMI (using main pipeline script)"
echo "Verifying filtered SNP list configuration"
echo "=========================================="

SCRIPT_DIR="/path/to/project/regenie_pipeline/scripts"
REGENIE_SCRIPT="${SCRIPT_DIR}/run_regenie_array.sh"

# Test with full stratum first (largest sample, most likely to catch issues)
STRATUM="full"
SAMPLE_TYPE="hallmarks_heldout"
PHENO_NAME="BMI"

echo ""
echo "Testing: ${STRATUM} - ${SAMPLE_TYPE} - ${PHENO_NAME}"
echo ""

# Verify filtered SNP list exists
FILTERED_SNP_LIST="/path/to/project/regenie_pipeline/filtered_snps/chr22_maf0.001.snplist"
if [ ! -f "$FILTERED_SNP_LIST" ]; then
    echo "ERROR: Required filtered SNP list not found: $FILTERED_SNP_LIST"
    echo "  Run prefilter_snps.sh first"
    exit 1
fi

SNP_COUNT=$(wc -l < "$FILTERED_SNP_LIST")
echo "✓ Verified filtered SNP list exists"
echo "  File: $FILTERED_SNP_LIST"
echo "  SNP count: $SNP_COUNT"
echo ""

# Run the test using the main pipeline script
echo "Running REGENIE test..."
echo ""

"$REGENIE_SCRIPT" "$STRATUM" "$SAMPLE_TYPE" "$PHENO_NAME"

EXIT_CODE=$?

echo ""
echo "=========================================="
if [ $EXIT_CODE -eq 0 ]; then
    echo "✓✓✓ Test PASSED ✓✓✓"
    echo "  The pipeline is correctly configured with the filtered SNP list"
    echo "  Ready to submit full hallmarks pipeline"
else
    echo "✗✗✗ Test FAILED ✗✗✗"
    echo "  Exit code: $EXIT_CODE"
    echo "  Check logs for details"
fi
echo "=========================================="

exit $EXIT_CODE
