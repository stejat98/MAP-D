#!/bin/bash
#SBATCH -c 8
#SBATCH -t 2-00:00:00
#SBATCH --mem=50G
#SBATCH -p medium
# Modified version for full strata: processes ONE chromosome per job (to avoid timeouts)
# Each chromosome takes ~13 hours, so 2-day limit is safe
#SBATCH -o /n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_full_perchr_%A_%a.out
#SBATCH -e /n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_full_perchr_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sivateja_tangirala@hms.harvard.edu

# REGENIE GWAS array job - PER CHROMOSOME VERSION (for full strata only)
# Runs Step 1 (null model) once per phenotype, then Step 2 for ONE chromosome only
# Format: stratum|sample_type|pheno|chromosome

# Parse command line arguments or use SLURM_ARRAY_TASK_ID
if [ -z "$1" ]; then
    # Get phenotype info from array task ID
    PHENO_LIST_FILE="${PHENO_LIST_FILE:-/n/groups/patel/sivateja/regenie_pipeline/scripts/pheno_list_full_per_chr.txt}"
    if [ ! -f "$PHENO_LIST_FILE" ]; then
        echo "ERROR: Phenotype list file not found: $PHENO_LIST_FILE"
        exit 1
    fi
    PHENO_LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$PHENO_LIST_FILE")
    STRATUM=$(echo "$PHENO_LINE" | cut -d'|' -f1)
    SAMPLE_TYPE=$(echo "$PHENO_LINE" | cut -d'|' -f2)
    PHENO_NAME=$(echo "$PHENO_LINE" | cut -d'|' -f3)
    TARGET_CHR=$(echo "$PHENO_LINE" | cut -d'|' -f4)
else
    # Direct specification (for testing)
    STRATUM="$1"
    SAMPLE_TYPE="$2"
    PHENO_NAME="$3"
    TARGET_CHR="$4"
fi

echo "=========================================="
echo "REGENIE GWAS - PER CHROMOSOME (Full Strata)"
echo "Stratum: $STRATUM"
echo "Sample type: $SAMPLE_TYPE"
echo "Phenotype: $PHENO_NAME"
echo "Chromosome: $TARGET_CHR"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID:-N/A}"
echo "=========================================="

# Paths
CONDA_ENV="${CONDA_ENV:-/n/groups/patel/sivateja/.conda/envs/regenie_env}"
GENO_PATH="/n/no_backup2/patel/uk_biobank/ukb_genetics/22881"
FAM_PATH="${GENO_PATH}/bgen_converted"
INPUT_PATH="/n/groups/patel/sivateja/regenie_pipeline/inputs/${STRATUM}/${SAMPLE_TYPE}/${PHENO_NAME}"
OUTPUT_ROOT="/n/scratch/users/s/st320/regenie"

# Per-phenotype output dirs & prefixes
OUTDIR_STEP1="${OUTPUT_ROOT}/${PHENO_NAME}/${STRATUM}/step1"
OUTDIR_STEP2="${OUTPUT_ROOT}/${PHENO_NAME}/${STRATUM}/step2"
TMPDIR_STEP1="${OUTDIR_STEP1}/tmpdir_${SLURM_JOB_ID:-$$}"

# Create output directories
mkdir -p "${OUTDIR_STEP1}" "${OUTDIR_STEP2}" "${TMPDIR_STEP1}"

# Check input files exist
PHENO_FILE="${INPUT_PATH}/pheno.txt"
COVAR_FILE="${INPUT_PATH}/covar.txt"

if [ ! -f "$PHENO_FILE" ]; then
    echo "ERROR: Phenotype file not found: $PHENO_FILE"
    exit 1
fi

if [ ! -f "$COVAR_FILE" ]; then
    echo "ERROR: Covariate file not found: $COVAR_FILE"
    exit 1
fi

# Initialize module system (required for SLURM jobs)
if [ -f /etc/profile.d/modules.sh ]; then
    source /etc/profile.d/modules.sh
fi

# Conda setup
export HOME=/n/groups/patel/sivateja
export CONDA_NO_PLUGINS=true
export CONDA_PKGS_DIRS=/n/groups/patel/sivateja/.conda/pkgs
export CONDA_ENVS_PATH=/n/groups/patel/sivateja/.conda/envs

module load conda/miniforge3/24.11.3-0
eval "$(conda shell.bash hook)"
conda activate "${CONDA_ENV}"

# Check if Step 1 is already completed
PRED_FILE="${OUTDIR_STEP1}/regenie_step1_${PHENO_NAME}_pred.list"
STEP1_COMPLETE=0

if [ -f "$PRED_FILE" ] && [ -s "$PRED_FILE" ]; then
    echo "[${PHENO_NAME}] Step 1 already completed, skipping..."
    echo "  Using existing prediction file: $PRED_FILE"
    STEP1_COMPLETE=1
fi

# Determine if binary trait
BINARY_FLAG=""
if [ -n "$BINARY_TRAIT" ] && [ "$BINARY_TRAIT" = "1" ]; then
    BINARY_FLAG="--bt"
fi

# --- Step 1: null model (only run if not completed) ---
# Use lock file to ensure Step 1 only runs once even if multiple chromosome jobs start simultaneously
STEP1_LOCK="${OUTDIR_STEP1}/step1_lock"
if [ $STEP1_COMPLETE -eq 0 ]; then
    # Try to acquire lock (only first job will succeed)
    if mkdir "$STEP1_LOCK" 2>/dev/null; then
        echo "[${PHENO_NAME}] Acquired Step 1 lock, running Step 1..."
        
        # Create sample file from FAM file
        SAMPLE_FILE="${TMPDIR_STEP1}/sample_file.sample"
        awk 'BEGIN {print "ID_1 ID_2 missing sex"; print "0 0 0 D"} {print $1, $2, "0", $5}' "${FAM_PATH}/ukb22.fam" > "$SAMPLE_FILE"
        
        # Use pre-filtered SNP list
        FILTERED_SNP_LIST="/n/groups/patel/sivateja/regenie_pipeline/filtered_snps/chr22_maf0.001.snplist"
        if [ ! -f "$FILTERED_SNP_LIST" ]; then
            FILTERED_SNP_LIST="/n/groups/patel/sivateja/regenie_pipeline/filtered_snps/chr22_mac50.snplist"
        fi
        if [ ! -f "$FILTERED_SNP_LIST" ]; then
            FILTERED_SNP_LIST="/n/groups/patel/sivateja/regenie_pipeline/filtered_snps/chr22_mac100.snplist"
        fi
        
        EXCLUDE_FILE="${TMPDIR_STEP1}/exclude_snps.txt"
        touch "$EXCLUDE_FILE"
        
        if [ ! -f "$FILTERED_SNP_LIST" ]; then
            echo "ERROR: Pre-filtered SNP list not found: $FILTERED_SNP_LIST"
            rmdir "$STEP1_LOCK" 2>/dev/null
            exit 1
        else
            echo "[${PHENO_NAME}] ✓ Using pre-filtered SNP list: $FILTERED_SNP_LIST"
            EXTRACT_FLAG="--extract ${FILTERED_SNP_LIST}"
        fi
        
        # Run Step 1 (simplified - full retry logic in original script)
        regenie \
          --step 1 \
          --bgen "${GENO_PATH}/ukb_imp_chr22_v3.bgen" \
          --sample "${SAMPLE_FILE}" \
          --phenoFile "${PHENO_FILE}" \
          --covarFile "${COVAR_FILE}" \
          ${BINARY_FLAG} \
          ${EXTRACT_FLAG} \
          $(if [ -s "$EXCLUDE_FILE" ]; then echo "--exclude ${EXCLUDE_FILE}"; fi) \
          --bsize 1000 \
          --threads 8 \
          --lowmem \
          --lowmem-prefix "${TMPDIR_STEP1}/regenie_tmp_preds" \
          --force-step1 \
          --out "${OUTDIR_STEP1}/regenie_step1_${PHENO_NAME}" 2>&1 | tee "${OUTDIR_STEP1}/regenie_step1_${PHENO_NAME}.log"
        
        STEP1_EXIT=${PIPESTATUS[0]}
        
        if [ $STEP1_EXIT -ne 0 ]; then
            echo "ERROR: REGENIE step 1 failed with exit code $STEP1_EXIT"
            rmdir "$STEP1_LOCK" 2>/dev/null
            exit $STEP1_EXIT
        fi
        
        # Release lock
        rmdir "$STEP1_LOCK" 2>/dev/null
        echo "[${PHENO_NAME}] Step 1 completed, lock released"
    else
        echo "[${PHENO_NAME}] Step 1 lock exists, waiting for another job to complete Step 1..."
        # Wait for Step 1 to complete (with timeout)
        MAX_WAIT=7200  # 2 hours max wait
        WAIT_COUNT=0
        while [ ! -f "$PRED_FILE" ] && [ $WAIT_COUNT -lt $MAX_WAIT ]; do
            sleep 30
            WAIT_COUNT=$((WAIT_COUNT + 30))
            if [ $((WAIT_COUNT % 300)) -eq 0 ]; then
                echo "[${PHENO_NAME}] Still waiting for Step 1... ($((WAIT_COUNT/60)) minutes)"
            fi
        done
        
        if [ ! -f "$PRED_FILE" ]; then
            echo "ERROR: Step 1 did not complete within timeout period (2 hours)"
            exit 1
        fi
        echo "[${PHENO_NAME}] Step 1 completed by another job, proceeding..."
    fi
fi

# Check that prediction file exists
if [ ! -f "$PRED_FILE" ]; then
    echo "ERROR: Prediction file not found: $PRED_FILE"
    exit 1
fi

# --- Step 2: association testing (ONE chromosome only) ---
echo "[${PHENO_NAME}] Running REGENIE step 2 for chromosome ${TARGET_CHR} only"
echo "  Input: ${PHENO_FILE}"
echo "  Pred: ${PRED_FILE}"
echo "  Output: ${OUTDIR_STEP2}/regenie_step2_${PHENO_NAME}_chr${TARGET_CHR}"

STEP2_FLAGS=""
if [ -n "$BINARY_FLAG" ]; then
    STEP2_FLAGS="--bt --firth --approx --pThresh 0.01"
fi

# Process only the target chromosome
CHR_BED="${FAM_PATH}/ukb${TARGET_CHR}.bed"
CHR_BIM="${FAM_PATH}/ukb${TARGET_CHR}.bim"
CHR_FAM="${FAM_PATH}/ukb${TARGET_CHR}.fam"

if [ ! -f "$CHR_BED" ] || [ ! -f "$CHR_BIM" ] || [ ! -f "$CHR_FAM" ]; then
    echo "ERROR: Chromosome ${TARGET_CHR} files not found"
    echo "  Expected: $CHR_BED, $CHR_BIM, $CHR_FAM"
    exit 1
fi

echo ""
echo "[${PHENO_NAME}] Processing chromosome ${TARGET_CHR}..."

# Use MAC threshold consistent with earlier analyses (default 5)
# All earlier successful jobs used default MAC 5
MIN_MAC_THRESHOLD=5  # Default REGENIE value, consistent with earlier analyses

regenie \
  --step 2 \
  --bed "${FAM_PATH}/ukb${TARGET_CHR}" \
  --phenoFile "${PHENO_FILE}" \
  --covarFile "${COVAR_FILE}" \
  --pred "${PRED_FILE}" \
  ${STEP2_FLAGS} \
  --minMAC ${MIN_MAC_THRESHOLD} \
  --bsize 400 \
  --threads 8 \
  --out "${OUTDIR_STEP2}/regenie_step2_${PHENO_NAME}_chr${TARGET_CHR}" 2>&1 | tee "${OUTDIR_STEP2}/regenie_step2_${PHENO_NAME}_chr${TARGET_CHR}.log"

CHR_EXIT=${PIPESTATUS[0]}
if [ $CHR_EXIT -ne 0 ]; then
    echo "ERROR: REGENIE step 2 failed for chromosome ${TARGET_CHR} with exit code $CHR_EXIT"
    exit $CHR_EXIT
else
    echo "[${PHENO_NAME}] Chromosome ${TARGET_CHR} complete"
fi

echo ""
echo "[${PHENO_NAME}] REGENIE GWAS complete for chromosome ${TARGET_CHR}"
echo "  Output file: ${OUTDIR_STEP2}/regenie_step2_${PHENO_NAME}_chr${TARGET_CHR}.regenie.gz"
