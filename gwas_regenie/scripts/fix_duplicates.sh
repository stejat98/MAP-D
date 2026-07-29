#!/bin/bash
# Fix duplicate individuals in existing REGENIE phenotype and covariate files
# Keeps first occurrence of each duplicate individual

SCRIPT_DIR="/n/groups/patel/sivateja/regenie_pipeline/scripts"
INPUTS_DIR="/n/groups/patel/sivateja/regenie_pipeline/inputs"

echo "=========================================="
echo "Fixing duplicate individuals in REGENIE input files"
echo "=========================================="

# Process all phenotype files
find "$INPUTS_DIR" -name "pheno.txt" | while read pheno_file; do
    dir=$(dirname "$pheno_file")
    echo "Processing: $dir"
    
    # Check for duplicates
    n_duplicates=$(awk 'NR>1 {print $1, $2}' "$pheno_file" | sort | uniq -d | wc -l)
    
    if [ "$n_duplicates" -gt 0 ]; then
        echo "  Found $n_duplicates duplicate individuals, fixing..."
        
        # Create backup
        cp "$pheno_file" "${pheno_file}.backup"
        
        # Remove duplicates (keep first occurrence)
        # Use awk to track seen FID/IID pairs
        awk '
        BEGIN {OFS=" "}
        NR==1 {print; next}
        {
            key = $1 "_" $2
            if (!seen[key]) {
                print
                seen[key] = 1
            }
        }' "$pheno_file" > "${pheno_file}.tmp" && mv "${pheno_file}.tmp" "$pheno_file"
        
        # Fix corresponding covar.txt file
        covar_file="$dir/covar.txt"
        if [ -f "$covar_file" ]; then
            # Remove duplicates from covar file (keep first occurrence)
            cp "$covar_file" "${covar_file}.backup"
            awk '
            BEGIN {OFS=" "}
            NR==1 {print; next}
            {
                key = $1 "_" $2
                if (!seen[key]) {
                    print
                    seen[key] = 1
                }
            }' "$covar_file" > "${covar_file}.tmp" && mv "${covar_file}.tmp" "$covar_file"
            
            echo "  Fixed covar.txt"
        fi
        
        echo "  Fixed pheno.txt (removed duplicates)"
    else
        echo "  No duplicates found"
    fi
done

echo ""
echo "=========================================="
echo "Duplicate fixing complete!"
echo "=========================================="
echo "Backup files created with .backup extension"
echo "You can now resubmit your REGENIE jobs"
