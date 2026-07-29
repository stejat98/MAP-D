#!/bin/bash
# Submit chunked MR jobs for cis and trans pQTLs

# Configuration
TOTAL_CHUNKS=8  # Split proteins into 8 chunks
SCRIPT_DIR="/n/groups/patel/sivateja/regenie_pipeline/scripts"
SLURM_DIR="/n/groups/patel/sivateja/regenie_pipeline/slurm"
OUTPUT_DIR="/n/groups/patel/sivateja/regenie_pipeline/results/twosampleMR/cis_trans_chunked"

# Create directories
mkdir -p "$SLURM_DIR"
mkdir -p "$OUTPUT_DIR"

echo "=== Submitting Chunked MR Jobs ==="
echo "Total chunks: $TOTAL_CHUNKS"
echo "QTL types: cis, trans"
echo "Total jobs: $((TOTAL_CHUNKS * 2))"
echo ""

# Submit jobs for each qtl_type and chunk
for qtl_type in cis trans; do
    echo "--- Submitting $qtl_type jobs ---"
    
    for chunk_id in $(seq 1 $TOTAL_CHUNKS); do
        job_name="MR_${qtl_type}_c${chunk_id}"
        
        sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=${job_name}
#SBATCH --partition=short
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH --output=${SLURM_DIR}/MR_${qtl_type}_chunk${chunk_id}_%j.out
#SBATCH --error=${SLURM_DIR}/MR_${qtl_type}_chunk${chunk_id}_%j.err

# Load modules
module load gcc/14.2.0
module load R/4.4.2

export R_LIBS_USER="/n/groups/patel/sivateja/.R/library"

cd /n/groups/patel/sivateja/regenie_pipeline

echo "Starting MR analysis: ${qtl_type} chunk ${chunk_id}/${TOTAL_CHUNKS}"
echo "Date: \$(date)"
echo ""

Rscript ${SCRIPT_DIR}/twosampleMR_cis_trans_chunked.R ${qtl_type} ${chunk_id} ${TOTAL_CHUNKS}

echo ""
echo "Completed: \$(date)"
EOF
        
        echo "  Submitted: $job_name (chunk $chunk_id/$TOTAL_CHUNKS)"
    done
    echo ""
done

echo "=== All jobs submitted ==="
echo ""
echo "Monitor with: squeue -u \$USER"
echo "Results will be saved to: $OUTPUT_DIR"
