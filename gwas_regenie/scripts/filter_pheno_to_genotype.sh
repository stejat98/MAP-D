#!/bin/bash
# Filter phenotype and covariate files to only include IIDs present in genotype file
# This ensures REGENIE can match individuals between files

set -e

if [ $# -lt 4 ]; then
    echo "Usage: $0 <pheno_file> <covar_file> <psam_file> <output_dir>"
    echo "  Creates filtered versions of pheno and covar files with only IIDs in psam"
    exit 1
fi

PHENO_FILE="$1"
COVAR_FILE="$2"
PSAM_FILE="$3"
OUTPUT_DIR="$4"

mkdir -p "$OUTPUT_DIR"

# Extract valid IIDs from psam file
VALID_IIDS="${OUTPUT_DIR}/valid_iids.txt"
awk -F'\t' 'NR>1 && $2>0 {print $2}' "$PSAM_FILE" | sort -n > "$VALID_IIDS"
N_VALID=$(wc -l < "$VALID_IIDS")
echo "Found $N_VALID valid IIDs in genotype file"

# Create filtered phenotype file
FILTERED_PHENO="${OUTPUT_DIR}/pheno_filtered.txt"
head -1 "$PHENO_FILE" > "$FILTERED_PHENO"
awk -v valid_iids="$VALID_IIDS" '
    BEGIN {
        while ((getline line < valid_iids) > 0) {
            valid[line] = 1
        }
        close(valid_iids)
    }
    NR > 1 && $2 in valid {
        print
    }
' "$PHENO_FILE" >> "$FILTERED_PHENO"

# Create filtered covariate file  
FILTERED_COVAR="${OUTPUT_DIR}/covar_filtered.txt"
head -1 "$COVAR_FILE" > "$FILTERED_COVAR"
awk -v valid_iids="$VALID_IIDS" '
    BEGIN {
        while ((getline line < valid_iids) > 0) {
            valid[line] = 1
        }
        close(valid_iids)
    }
    NR > 1 && $2 in valid {
        print
    }
' "$COVAR_FILE" >> "$FILTERED_COVAR"

N_PHENO=$(tail -n +2 "$FILTERED_PHENO" | wc -l)
N_COVAR=$(tail -n +2 "$FILTERED_COVAR" | wc -l)
echo "Filtered phenotype: $N_PHENO individuals"
echo "Filtered covariates: $N_COVAR individuals"

if [ $N_PHENO -eq 0 ]; then
    echo "ERROR: No individuals matched! Check IID format."
    exit 1
fi

echo "Output files:"
echo "  $FILTERED_PHENO"
echo "  $FILTERED_COVAR"
