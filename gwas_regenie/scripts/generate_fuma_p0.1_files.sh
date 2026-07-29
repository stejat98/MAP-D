#!/bin/bash
#
# Generate FUMA p-value < 0.1 filtered files
# Usage: ./generate_fuma_p0.1_files.sh [--force] [input_file.fuma.gz]
#
# If no input file is specified, processes all .fuma.gz files in results/
# that don't already have a p0.1 version.
#
# Options:
#   --force    Regenerate p0.1 files even if they already exist
#

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR="${SCRIPT_DIR}/../results"

FORCE=false
INPUT_FILE=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --force)
            FORCE=true
            shift
            ;;
        *)
            INPUT_FILE="$1"
            shift
            ;;
    esac
done

# Function to filter FUMA file by p-value < 0.1
filter_pvalue() {
    local input_file="$1"
    local output_file="${input_file%.fuma.gz}.p0.1.fuma.gz"
    
    # Check if output already exists (unless force flag is set)
    if [[ -f "$output_file" && "$FORCE" != "true" ]]; then
        echo "  Skipping: ${output_file} already exists (use --force to regenerate)"
        return 0
    fi
    
    echo "  Processing: $(basename "$input_file")"
    echo "    Output: $(basename "$output_file")"
    
    # Filter to P < 0.1 (column 6 is P-value in FUMA format)
    # NR==1 keeps the header, $6 < 0.1 filters by p-value
    zcat "$input_file" | awk 'NR==1 || $6 < 0.1' | gzip > "$output_file"
    
    # Report sizes
    local input_size=$(ls -lh "$input_file" | awk '{print $5}')
    local output_size=$(ls -lh "$output_file" | awk '{print $5}')
    local input_lines=$(zcat "$input_file" | wc -l)
    local output_lines=$(zcat "$output_file" | wc -l)
    
    echo "    Input:  $input_size ($((input_lines - 1)) variants)"
    echo "    Output: $output_size ($((output_lines - 1)) variants)"
    echo "    Reduction: $(echo "scale=1; (1 - ($output_lines / $input_lines)) * 100" | bc)%"
    echo ""
}

echo "=============================================="
echo "FUMA P-value Filter (P < 0.1)"
echo "=============================================="
echo ""

if [[ -n "$INPUT_FILE" ]]; then
    # Process single file
    if [[ ! -f "$INPUT_FILE" ]]; then
        echo "Error: File not found: $INPUT_FILE"
        exit 1
    fi
    filter_pvalue "$INPUT_FILE"
else
    # Process all .fuma.gz files in results directory
    echo "Scanning ${RESULTS_DIR} for FUMA files..."
    echo ""
    
    # Find all .fuma.gz files that are NOT p0.1 files
    files_found=0
    files_processed=0
    
    for fuma_file in "${RESULTS_DIR}"/*.fuma.gz; do
        # Skip if no files match
        [[ -e "$fuma_file" ]] || continue
        
        # Skip p0.1 files
        if [[ "$fuma_file" == *".p0.1.fuma.gz" ]]; then
            continue
        fi
        
        files_found=$((files_found + 1))
        filter_pvalue "$fuma_file"
        files_processed=$((files_processed + 1))
    done
    
    echo "=============================================="
    echo "Summary"
    echo "=============================================="
    echo "Files found: $files_found"
    echo "Files processed: $files_processed"
fi

echo ""
echo "Done!"
