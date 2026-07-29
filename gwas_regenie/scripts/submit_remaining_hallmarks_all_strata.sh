#!/bin/bash
# Submit REGENIE jobs for remaining hallmarks (HDL, LDL, TRIG_HDL_RATIO, systolic_BP, diastolic_BP)
# across all strata (diabetes, prediabetes, full) using chromosome-specific configuration
# Follows the same pattern that worked successfully for BMI and HbA1c

SCRIPT_DIR="/n/groups/patel/sivateja/regenie_pipeline/scripts"

echo "=========================================="
echo "Submitting REGENIE jobs for remaining hallmarks"
echo "=========================================="
echo "Hallmarks: HDL, LDL, TRIG_HDL_RATIO, systolic_BP, diastolic_BP"
echo "Strata: diabetes, prediabetes, full"
echo "Configuration: Chromosome-specific (same as BMI/HbA1c)"
echo ""

# Submit diabetes first (shortest, highest priority partition)
echo "=========================================="
echo "1. Submitting diabetes stratum..."
echo "=========================================="
bash "${SCRIPT_DIR}/submit_diabetes_per_chr_remaining.sh"
DIAB_EXIT=$?

if [ $DIAB_EXIT -ne 0 ]; then
    echo "WARNING: Diabetes submission failed with exit code $DIAB_EXIT"
fi

echo ""
echo "=========================================="
echo "2. Submitting prediabetes stratum..."
echo "=========================================="
bash "${SCRIPT_DIR}/submit_prediabetes_per_chr_remaining.sh"
PREDIAB_EXIT=$?

if [ $PREDIAB_EXIT -ne 0 ]; then
    echo "WARNING: Prediabetes submission failed with exit code $PREDIAB_EXIT"
fi

echo ""
echo "=========================================="
echo "3. Submitting full stratum..."
echo "=========================================="
bash "${SCRIPT_DIR}/submit_full_per_chr_remaining.sh"
FULL_EXIT=$?

if [ $FULL_EXIT -ne 0 ]; then
    echo "WARNING: Full stratum submission failed with exit code $FULL_EXIT"
fi

echo ""
echo "=========================================="
echo "Job submission summary"
echo "=========================================="
if [ $DIAB_EXIT -eq 0 ] && [ $PREDIAB_EXIT -eq 0 ] && [ $FULL_EXIT -eq 0 ]; then
    echo "  ✓ All strata submitted successfully"
else
    echo "  ⚠ Some submissions had issues (check above for details)"
fi
echo ""
echo "Total jobs submitted:"
echo "  - Diabetes: 110 jobs (5 phenotypes × 22 chromosomes)"
echo "  - Prediabetes: 110 jobs (5 phenotypes × 22 chromosomes)"
echo "  - Full: 110 jobs (5 phenotypes × 22 chromosomes)"
echo "  - Total: 330 jobs"
echo ""
echo "Monitor jobs:"
echo "  squeue -u st320"
echo "  squeue -u st320 -p short | grep regenie.*perchr_remaining"
echo "  squeue -u st320 -p medium | grep regenie.*perchr_remaining"
echo ""
echo "Logs: /n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_*_perchr_remaining_*.out"
