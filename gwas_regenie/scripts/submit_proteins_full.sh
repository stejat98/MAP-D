#!/bin/bash
# Submit REGENIE jobs for the protein GWAS in the Olink (proteomics) sample.
# Full cohort only. Uses the same covariate approach as the hallmark GWAS.

PIPE="${REGENIE_PIPELINE_ROOT:-/path/to/project/regenie_pipeline}"
SCRIPT_DIR="${PIPE}/scripts"
REGENIE_SCRIPT="${SCRIPT_DIR}/run_regenie_array.sh"

PROTEINS_DIR="${PIPE}/inputs/full/proteins_only"

echo "=========================================="
echo "Submitting REGENIE jobs for proteins"
echo "=========================================="

PROTEINS=($(ls "$PROTEINS_DIR" | sort))
echo "Found ${#PROTEINS[@]} proteins"
echo ""

# Build the pheno list consumed by the array job
PHENO_LIST="${SCRIPT_DIR}/pheno_list_proteins.txt"
> "$PHENO_LIST"

JOB_COUNT=0
for PROTEIN in "${PROTEINS[@]}"; do
    INPUT_PATH="${PROTEINS_DIR}/${PROTEIN}"
    if [ -f "${INPUT_PATH}/pheno.txt" ] && [ -f "${INPUT_PATH}/covar.txt" ]; then
        echo "full|proteins_only|${PROTEIN}" >> "$PHENO_LIST"
        ((JOB_COUNT++))
    fi
done

echo "Total jobs to submit: $JOB_COUNT"
echo "Pheno list file: $PHENO_LIST"
echo ""

if [ "$JOB_COUNT" -eq 0 ]; then
    echo "Nothing to submit. Run generate_inputs.R first."
    exit 1
fi

# Step 1 ~8-10 h, Step 2 ~22-44 h across all chromosomes.
sbatch \
  --array=1-${JOB_COUNT} \
  --time=3-00:00:00 \
  --cpus-per-task=8 \
  --mem=50G \
  --partition=medium \
  --job-name=regenie_proteins_full \
  --output=${PIPE}/slurm/regenie_proteins_full_%A_%a.out \
  --error=${PIPE}/slurm/regenie_proteins_full_%A_%a.err \
  --export=PHENO_LIST_FILE="$PHENO_LIST" \
  "$REGENIE_SCRIPT"

if [ $? -eq 0 ]; then
    echo "  ✓ Protein jobs submitted"
else
    echo "  ✗ Failed to submit protein jobs"
    exit 1
fi

echo ""
echo "Monitor jobs:  squeue -u $USER -p medium"
echo "Logs:          ${PIPE}/slurm/"
echo ""
