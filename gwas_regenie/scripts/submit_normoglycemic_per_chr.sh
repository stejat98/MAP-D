#!/bin/bash
# Submit REGENIE jobs for normoglycemic BMI GWAS (held-out sample) split by chromosome

SCRIPT_DIR="/n/groups/patel/sivateja/regenie_pipeline/scripts"
REGENIE_SCRIPT="${SCRIPT_DIR}/run_regenie_array_per_chr.sh"
PHENO_LIST="${SCRIPT_DIR}/pheno_list_normoglycemic_per_chr.txt"

echo "=========================================="
echo "Submitting REGENIE jobs for Normoglycemic BMI GWAS"
echo "=========================================="
echo "Phenotype: BMI"
echo "Stratum: Normoglycemic (held-out, no proteomics)"
echo "Chromosomes: 1-22 (22 jobs total)"
echo ""

# Check pheno list exists
if [ ! -f "$PHENO_LIST" ]; then
    echo "ERROR: Pheno list file not found: $PHENO_LIST"
    echo "Run generate_inputs_normoglycemic.R first"
    exit 1
fi

# Check input files exist
INPUT_DIR="/n/groups/patel/sivateja/regenie_pipeline/inputs/normoglycemic/hallmarks_heldout/BMI"
if [ ! -f "${INPUT_DIR}/pheno.txt" ] || [ ! -f "${INPUT_DIR}/covar.txt" ]; then
    echo "ERROR: Input files not found in ${INPUT_DIR}"
    echo "Run generate_inputs_normoglycemic.R first"
    exit 1
fi

TOTAL_JOBS=$(wc -l < "$PHENO_LIST")
echo "Total jobs to submit: $TOTAL_JOBS"
echo "Input files: ${INPUT_DIR}/"
echo ""

# Verify sample size
SAMPLE_SIZE=$(wc -l < "${INPUT_DIR}/pheno.txt")
echo "Sample size: $((SAMPLE_SIZE - 1)) individuals"
echo ""

echo "=========================================="
echo "Submitting jobs..."
echo "=========================================="

sbatch \
  --array=1-${TOTAL_JOBS} \
  --time=2-00:00:00 \
  --cpus-per-task=8 \
  --mem=50G \
  --partition=medium \
  --job-name=regenie_normo_BMI \
  --output=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_normo_BMI_%A_%a.out \
  --error=/n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_normo_BMI_%A_%a.err \
  --export=PHENO_LIST_FILE="$PHENO_LIST" \
  "$REGENIE_SCRIPT"

if [ $? -eq 0 ]; then
    echo ""
    echo "✓ Normoglycemic BMI GWAS jobs submitted!"
    echo ""
    echo "Summary:"
    echo "  - Stratum: Normoglycemic (HbA1c < 39 mmol/mol)"
    echo "  - Sample: Held-out (no proteomics)"
    echo "  - Jobs: 22 (one per chromosome)"
    echo "  - Time limit: 2 days per job"
    echo ""
    echo "Monitor: squeue -u \$USER | grep regenie_normo"
    echo "Logs: slurm/regenie_normo_BMI_*.out"
fi
