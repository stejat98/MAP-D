#!/bin/bash
# Submit REGENIE for the 7 pre-spec Olink proteins (full | proteins_only), BMI in covar.
# Inputs: scripts/pheno_list_prespec7_protein_bmi_full.txt (from generate_regenie_inputs_prespec7_protein_bmi.R)
# Outputs: REGENIE_OUTPUT_ROOT / ${PHENO_NAME} / full / {step1,step2} (default scratch below).
#
# Usage: bash scripts/submit_regenie_prespec7_proteins_full.sh
#     or:  sbatch scripts/submit_regenie_prespec7_proteins_full.sh   # if the shebang job block is used

set -euo pipefail
SCRIPT_DIR="/path/to/project/regenie_pipeline/scripts"
REGENIE_SCRIPT="${SCRIPT_DIR}/run_regenie_array.sh"
PHENO_LIST_FILE="${PHENO_LIST_FILE:-${SCRIPT_DIR}/pheno_list_prespec7_protein_bmi_full.txt}"
PIPE="/path/to/project/regenie_pipeline"

if [ ! -f "$PHENO_LIST_FILE" ]; then
  echo "ERROR: Pheno list not found: $PHENO_LIST_FILE"
  echo "       Run: module load gcc/14.2.0 R/4.4.2 && Rscript ${SCRIPT_DIR}/generate_regenie_inputs_prespec7_protein_bmi.R"
  exit 1
fi

VERIFY_ROOT="${REGENIE_INPUT_ROOT:-$PIPE/inputs}"
while IFS= read -r line; do
  stratum="${line%%|*}"
  rest="${line#*|}"
  sample="${rest%%|*}"
  phen="${rest#*|}"
  pdir="${VERIFY_ROOT}/${stratum}/${sample}/${phen}"
  if [ ! -f "${pdir}/pheno.txt" ] || [ ! -f "${pdir}/covar.txt" ]; then
    echo "ERROR: Missing inputs under ${pdir}"
    exit 1
  fi
done < "$PHENO_LIST_FILE"

TOTAL_JOBS=$(wc -l < "$PHENO_LIST_FILE")
echo "Submitting ${TOTAL_JOBS} jobs (full | proteins_only | prespec7 proteins)"
echo "  PHENO_LIST_FILE=$PHENO_LIST_FILE"
echo "  First line: $(head -1 "$PHENO_LIST_FILE")"

# Match submit_proteins_full.sh: medium partition, 3d walltime
sbatch \
  --array=1-"${TOTAL_JOBS}" \
  --time=3-00:00:00 \
  --cpus-per-task=8 \
  --mem=50G \
  --partition=medium \
  --job-name=regenie_prespec7_pQTL \
  --output="${PIPE}/slurm/regenie_prespec7_protein_BMI_full_%A_%a.out" \
  --error="${PIPE}/slurm/regenie_prespec7_protein_BMI_full_%A_%a.err" \
  --export=PHENO_LIST_FILE="${PHENO_LIST_FILE}" \
  "$REGENIE_SCRIPT"

echo "Done. Check: squeue -u \$USER; logs under ${PIPE}/slurm/regenie_prespec7_protein_BMI_full_*"
