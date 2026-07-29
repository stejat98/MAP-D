#!/bin/bash
# TINY REGENIE test - uses bed format with only first 1000 variants
# Should run in 2-3 minutes

set -e

echo "=========================================="
echo "REGENIE Tiny Test - BMI (1000 variants)"
echo "=========================================="

# Configuration
STRATUM="full"
SAMPLE_TYPE="hallmarks_heldout"
PHENO_NAME="BMI"
CHR_TEST="22"

# Paths
INPUT_PATH="/n/groups/patel/sivateja/regenie_pipeline/inputs/${STRATUM}/${SAMPLE_TYPE}/${PHENO_NAME}"
OUTPUT_ROOT="/n/scratch/users/s/st320/regenie_test"
CONDA_ENV="/n/groups/patel/sivateja/.conda/envs/regenie_env"
GENO_PATH="/n/no_backup2/patel/uk_biobank/ukb_genetics/22881/bgen_converted"

# Output directories
OUTDIR_STEP1="${OUTPUT_ROOT}/${PHENO_NAME}/${STRATUM}/step1_tiny"
OUTDIR_STEP2="${OUTPUT_ROOT}/${PHENO_NAME}/${STRATUM}/step2_tiny"
TMPDIR_STEP1="${OUTDIR_STEP1}/tmpdir_${SLURM_JOB_ID:-$$}"

mkdir -p "${OUTDIR_STEP1}" "${OUTDIR_STEP2}" "${TMPDIR_STEP1}"

# Check input files
PHENO_FILE="${INPUT_PATH}/pheno.txt"
COVAR_FILE="${INPUT_PATH}/covar.txt"
BED_FILE="${GENO_PATH}/ukb${CHR_TEST}.bed"
BIM_FILE="${GENO_PATH}/ukb${CHR_TEST}.bim"
FAM_FILE="${GENO_PATH}/ukb${CHR_TEST}.fam"

if [ ! -f "$PHENO_FILE" ] || [ ! -f "$COVAR_FILE" ] || [ ! -f "$BED_FILE" ] || [ ! -f "$BIM_FILE" ] || [ ! -f "$FAM_FILE" ]; then
    echo "ERROR: Required input files not found"
    exit 1
fi

# Create minimal bed/bim/fam with only first 1000 variants
echo "Creating minimal genotype files (first 1000 variants)..."
TINY_BED="${TMPDIR_STEP1}/ukb${CHR_TEST}_tiny.bed"
TINY_BIM="${TMPDIR_STEP1}/ukb${CHR_TEST}_tiny.bim"
TINY_FAM="${TMPDIR_STEP1}/ukb${CHR_TEST}_tiny.fam"

# Copy first 1000 lines of bim and fam
head -1000 "$BIM_FILE" > "$TINY_BIM"
cp "$FAM_FILE" "$TINY_FAM"

# Extract first 1000 variants from bed file (bed is binary, need plink or similar)
# For now, use plink to subset
module load plink/1.90b6.26 2>/dev/null || echo "Plink not available, will use full bed"

# Create minimal phenotype/covariate (10000 individuals for better SNP coverage)
echo "Creating minimal phenotype files (10000 individuals)..."
MIN_PHENO="${TMPDIR_STEP1}/pheno_tiny.txt"
MIN_COVAR="${TMPDIR_STEP1}/covar_tiny.txt"

FIRST_10000_IIDS="${TMPDIR_STEP1}/first_10000_iids.txt"
head -10001 "$FAM_FILE" | tail -n +2 | awk '{print $2}' | sort -n > "$FIRST_10000_IIDS"

head -1 "$PHENO_FILE" > "$MIN_PHENO"
awk -v iids_file="$FIRST_10000_IIDS" '
    BEGIN {
        while ((getline line < iids_file) > 0) valid[line] = 1
        close(iids_file)
    }
    NR > 1 && $2 in valid { print; if (++count >= 10000) exit }
' "$PHENO_FILE" >> "$MIN_PHENO"

head -1 "$COVAR_FILE" > "$MIN_COVAR"
awk -v iids_file="$FIRST_10000_IIDS" '
    BEGIN {
        while ((getline line < iids_file) > 0) valid[line] = 1
        close(iids_file)
    }
    NR > 1 && $2 in valid { print; if (++count >= 10000) exit }
' "$COVAR_FILE" >> "$MIN_COVAR"

N_FILTERED=$(tail -n +2 "$MIN_PHENO" | wc -l)
echo "  Created files with $N_FILTERED individuals"

# Use plink to extract first 1000 variants
if command -v plink &> /dev/null; then
    echo "Extracting first 1000 variants using plink..."
    plink --bed "$BED_FILE" --bim "$BIM_FILE" --fam "$FAM_FILE" \
          --extract <(head -1000 "$BIM_FILE" | awk '{print $2}') \
          --make-bed --out "${TMPDIR_STEP1}/ukb${CHR_TEST}_tiny" \
          --silent --threads 2 2>/dev/null || {
        echo "Plink extraction failed, using full bed file (will be slower)"
        cp "$BED_FILE" "$TINY_BED"
        cp "$BIM_FILE" "$TINY_BIM"
        cp "$FAM_FILE" "$TINY_FAM"
    }
    TINY_BED="${TMPDIR_STEP1}/ukb${CHR_TEST}_tiny.bed"
    TINY_BIM="${TMPDIR_STEP1}/ukb${CHR_TEST}_tiny.bim"
    TINY_FAM="${TMPDIR_STEP1}/ukb${CHR_TEST}_tiny.fam"
else
    echo "Plink not available, using full bed file"
    TINY_BED="$BED_FILE"
    TINY_BIM="$BIM_FILE"
    TINY_FAM="$FAM_FILE"
fi

# Conda setup
export HOME=/n/groups/patel/sivateja
export CONDA_NO_PLUGINS=true
export CONDA_PKGS_DIRS=/n/groups/patel/sivateja/.conda/pkgs
export CONDA_ENVS_PATH=/n/groups/patel/sivateja/.conda/envs

module load conda/miniforge3/24.11.3-0
eval "$(conda shell.bash hook)"
conda activate "${CONDA_ENV}"

if ! command -v regenie &> /dev/null; then
    echo "ERROR: regenie command not found"
    exit 1
fi

echo ""
echo "REGENIE version: $(regenie --version)"
N_VARIANTS=$(wc -l < "$TINY_BIM")
echo "  Variants: $N_VARIANTS"
echo "  Individuals: $N_FILTERED"
echo ""

# Step 1
echo "=========================================="
echo "Step 1: Running null model"
echo "=========================================="

# Exclude known problematic SNPs
EXCLUDE_FILE="${TMPDIR_STEP1}/exclude_snps.txt"
cat > "$EXCLUDE_FILE" << EOF
rs587697622
rs587627608
rs587634452
EOF
echo "  Excluding $(wc -l < "$EXCLUDE_FILE") known problematic SNPs"

# Use position range to limit to first 1MB (fewer variants, common SNPs)
BP_START=16050000
BP_END=17050000

regenie \
  --step 1 \
  --bed "${TINY_BED%.bed}" \
  --phenoFile "${MIN_PHENO}" \
  --covarFile "${MIN_COVAR}" \
  --exclude "${EXCLUDE_FILE}" \
  --bp ${BP_START} ${BP_END} \
  --bsize 200 \
  --threads 2 \
  --lowmem \
  --lowmem-prefix "${TMPDIR_STEP1}/regenie_tmp_preds" \
  --force-step1 \
  --out "${OUTDIR_STEP1}/regenie_step1_${PHENO_NAME}" 2>&1 | tee "${OUTDIR_STEP1}/regenie_step1_${PHENO_NAME}.log"

STEP1_EXIT=${PIPESTATUS[0]}
if [ $STEP1_EXIT -ne 0 ]; then
    echo "ERROR: Step 1 failed"
    tail -20 "${OUTDIR_STEP1}/regenie_step1_${PHENO_NAME}.log"
    exit $STEP1_EXIT
fi

PRED_FILE="${OUTDIR_STEP1}/regenie_step1_${PHENO_NAME}_pred.list"
if [ ! -f "$PRED_FILE" ]; then
    echo "ERROR: Prediction file not created"
    exit 1
fi

echo "✓ Step 1 complete!"

# Step 2
echo ""
echo "=========================================="
echo "Step 2: Running association testing"
echo "=========================================="

regenie \
  --step 2 \
  --bed "${TINY_BED%.bed}" \
  --phenoFile "${MIN_PHENO}" \
  --covarFile "${MIN_COVAR}" \
  --pred "${PRED_FILE}" \
  --bp ${BP_START} ${BP_END} \
  --bsize 50 \
  --threads 2 \
  --out "${OUTDIR_STEP2}/regenie_step2_${PHENO_NAME}" 2>&1 | tee "${OUTDIR_STEP2}/regenie_step2_${PHENO_NAME}.log"

STEP2_EXIT=${PIPESTATUS[0]}
if [ $STEP2_EXIT -ne 0 ]; then
    echo "ERROR: Step 2 failed"
    tail -20 "${OUTDIR_STEP2}/regenie_step2_${PHENO_NAME}.log"
    exit $STEP2_EXIT
fi

OUTPUT_FILE="${OUTDIR_STEP2}/regenie_step2_${PHENO_NAME}.regenie.gz"
if [ -f "$OUTPUT_FILE" ]; then
    echo ""
    echo "✓✓✓ Test completed successfully! ✓✓✓"
    echo ""
    echo "Output: $OUTPUT_FILE"
    N_VARIANTS=$(zcat "$OUTPUT_FILE" | tail -n +2 | wc -l)
    echo "Variants tested: $N_VARIANTS"
    echo ""
    echo "First 5 results:"
    zcat "$OUTPUT_FILE" | head -6
else
    echo "WARNING: Output file not found"
    exit 1
fi
