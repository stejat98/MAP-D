#!/bin/bash
#SBATCH -c 8
#SBATCH -t 12:00:00
#SBATCH --mem=50G
#SBATCH -p short
#SBATCH -o /n/groups/patel/sivateja/regenie_pipeline/slurm/test_all_strata_%j.out
#SBATCH -e /n/groups/patel/sivateja/regenie_pipeline/slurm/test_all_strata_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sivateja_tangirala@hms.harvard.edu
#SBATCH --job-name=regenie_test_all

# Test REGENIE jobs for all three strata to verify fixes work
# Using protein_x.1004 for consistency

SCRIPT_DIR="/n/groups/patel/sivateja/regenie_pipeline/scripts"
PROTEIN="protein_x.1004"

echo "=========================================="
echo "Testing REGENIE for all strata"
echo "Protein: $PROTEIN"
echo "=========================================="

# Test full stratum
echo ""
echo "--- Testing FULL stratum ---"
/n/groups/patel/sivateja/regenie_pipeline/scripts/run_regenie_array.sh full proteins_only $PROTEIN
FULL_EXIT=$?

if [ $FULL_EXIT -eq 0 ]; then
    echo "✓ FULL stratum test PASSED"
else
    echo "✗ FULL stratum test FAILED (exit code: $FULL_EXIT)"
fi

# Test prediabetes stratum  
echo ""
echo "--- Testing PREDIABETES stratum ---"
/n/groups/patel/sivateja/regenie_pipeline/scripts/run_regenie_array.sh prediabetes proteins_only $PROTEIN
PREDIAB_EXIT=$?

if [ $PREDIAB_EXIT -eq 0 ]; then
    echo "✓ PREDIABETES stratum test PASSED"
else
    echo "✗ PREDIABETES stratum test FAILED (exit code: $PREDIAB_EXIT)"
fi

# Test diabetes stratum
echo ""
echo "--- Testing DIABETES stratum ---"
/n/groups/patel/sivateja/regenie_pipeline/scripts/run_regenie_array.sh diabetes proteins_only $PROTEIN
DIAB_EXIT=$?

if [ $DIAB_EXIT -eq 0 ]; then
    echo "✓ DIABETES stratum test PASSED"
else
    echo "✗ DIABETES stratum test FAILED (exit code: $DIAB_EXIT)"
fi

# Summary
echo ""
echo "=========================================="
echo "Test Summary:"
echo "=========================================="
echo "Full stratum:        $([ $FULL_EXIT -eq 0 ] && echo 'PASSED' || echo 'FAILED')"
echo "Prediabetes stratum:  $([ $PREDIAB_EXIT -eq 0 ] && echo 'PASSED' || echo 'FAILED')"
echo "Diabetes stratum:     $([ $DIAB_EXIT -eq 0 ] && echo 'PASSED' || echo 'FAILED')"

if [ $FULL_EXIT -eq 0 ] && [ $PREDIAB_EXIT -eq 0 ] && [ $DIAB_EXIT -eq 0 ]; then
    echo ""
    echo "✓ All tests PASSED - ready to submit full batch!"
    exit 0
else
    echo ""
    echo "✗ Some tests FAILED - check logs before submitting full batch"
    exit 1
fi
