#!/bin/bash

# =============================================================================
# MAP-D Analysis Runner Script
# Metabolic Atlas of Progression to Diabetes (MAP-D)
# =============================================================================

set -euo pipefail  # Exit on error, undefined variables, pipe failures

# Script configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_NAME="MAP-D"
LOG_DIR="${SCRIPT_DIR}/logs"
RESULTS_DIR="${SCRIPT_DIR}/results"

# Create necessary directories
mkdir -p "${LOG_DIR}"
mkdir -p "${RESULTS_DIR}"

# Logging function
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "${LOG_DIR}/analysis_$(date '+%Y%m%d').log"
}

# Error handling
error_exit() {
    log "ERROR: $1"
    exit 1
}

# Check if R is available
check_r_installation() {
    if ! command -v Rscript &> /dev/null; then
        error_exit "R is not installed or not in PATH"
    fi
    
    R_VERSION=$(Rscript --version 2>&1 | grep -oP 'R scripting front-end version \K[0-9]+\.[0-9]+\.[0-9]+')
    log "Found R version: ${R_VERSION}"
    
    # Check minimum R version (4.0.0)
    if ! Rscript -e "if (getRversion() < '4.0.0') quit(status=1)" 2>/dev/null; then
        error_exit "R version 4.0.0 or higher is required. Found: ${R_VERSION}"
    fi
}

# Install required packages
install_packages() {
    log "Checking and installing required R packages..."
    
    Rscript -e "
    source('requirements.R')
    install_mapd_packages()
    " || error_exit "Failed to install required packages"
    
    log "Package installation completed"
}

# Run tests
run_tests() {
    log "Running pipeline tests..."
    
    if [ -f "${SCRIPT_DIR}/tests/test_pipeline.R" ]; then
        Rscript "${SCRIPT_DIR}/tests/test_pipeline.R" || error_exit "Pipeline tests failed"
        log "All tests passed"
    else
        log "Warning: Test file not found, skipping tests"
    fi
}

# Main analysis function
run_analysis() {
    local full_pipeline=${1:-false}
    local skip_preprocessing=${2:-false}
    local output_dir=${3:-"${RESULTS_DIR}"}
    
    log "Starting ${PROJECT_NAME} analysis pipeline"
    log "Full pipeline: ${full_pipeline}"
    log "Skip preprocessing: ${skip_preprocessing}"
    log "Output directory: ${output_dir}"
    
    # Prepare R command arguments
    local r_args=""
    if [ "${full_pipeline}" = "true" ]; then
        r_args="${r_args} --full-pipeline"
    fi
    if [ "${skip_preprocessing}" = "true" ]; then
        r_args="${r_args} --skip-preprocessing"
    fi
    if [ "${output_dir}" != "${RESULTS_DIR}" ]; then
        r_args="${r_args} --output-dir=${output_dir}"
    fi
    
    # Run the main analysis
    log "Executing main analysis script..."
    cd "${SCRIPT_DIR}"
    
    if ! Rscript main_analysis.R ${r_args} 2>&1 | tee -a "${LOG_DIR}/analysis_$(date '+%Y%m%d').log"; then
        error_exit "Main analysis failed"
    fi
    
    log "Analysis completed successfully"
}

# Generate summary report
generate_report() {
    log "Generating analysis summary report..."
    
    # Check if results exist
    if [ ! -d "${RESULTS_DIR}" ] || [ -z "$(ls -A "${RESULTS_DIR}")" ]; then
        log "Warning: No results found to summarize"
        return
    fi
    
    # Create summary
    {
        echo "=================================="
        echo "${PROJECT_NAME} Analysis Summary"
        echo "Generated: $(date)"
        echo "=================================="
        echo ""
        
        echo "Analysis Configuration:"
        echo "- Script directory: ${SCRIPT_DIR}"
        echo "- Results directory: ${RESULTS_DIR}"
        echo "- Log directory: ${LOG_DIR}"
        echo ""
        
        echo "Output Files:"
        find "${RESULTS_DIR}" -type f -name "*.csv" -o -name "*.rds" -o -name "*.pdf" | \
            sort | while read -r file; do
            echo "- $(basename "${file}") ($(stat -f%z "${file}" 2>/dev/null || stat -c%s "${file}" 2>/dev/null) bytes)"
        done
        echo ""
        
        echo "Log Files:"
        find "${LOG_DIR}" -type f -name "*.log" | sort | while read -r file; do
            echo "- $(basename "${file}") ($(wc -l < "${file}") lines)"
        done
        
    } > "${RESULTS_DIR}/analysis_summary_$(date '+%Y%m%d_%H%M%S').txt"
    
    log "Summary report generated"
}

# Cleanup function
cleanup() {
    log "Cleaning up temporary files..."
    # Add any cleanup tasks here
    log "Cleanup completed"
}

# Print usage information
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Run the MAP-D (Metabolic Atlas of Progression to Diabetes) analysis pipeline.

OPTIONS:
    -f, --full-pipeline     Run complete pipeline from raw data
    -s, --skip-preprocessing Skip data preprocessing step
    -o, --output-dir DIR    Specify output directory (default: ./results)
    -t, --test-only         Run tests only, skip analysis
    -i, --install-only      Install packages only, skip analysis
    -h, --help              Show this help message

EXAMPLES:
    $0                      # Run standard analysis
    $0 -f                   # Run full pipeline from raw data
    $0 -s -o /path/to/out   # Skip preprocessing, custom output dir
    $0 -t                   # Run tests only
    $0 -i                   # Install packages only

REQUIREMENTS:
    - R (>= 4.0.0)
    - Required R packages (see requirements.R)
    - UK Biobank proteomics data
    - DrugBank database (optional)

For more information, see README.md and docs/METHODOLOGY.md
EOF
}

# Parse command line arguments
FULL_PIPELINE=false
SKIP_PREPROCESSING=false
OUTPUT_DIR="${RESULTS_DIR}"
TEST_ONLY=false
INSTALL_ONLY=false

while [[ $# -gt 0 ]]; do
    case $1 in
        -f|--full-pipeline)
            FULL_PIPELINE=true
            shift
            ;;
        -s|--skip-preprocessing)
            SKIP_PREPROCESSING=true
            shift
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -t|--test-only)
            TEST_ONLY=true
            shift
            ;;
        -i|--install-only)
            INSTALL_ONLY=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            log "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# Main execution
main() {
    log "Starting ${PROJECT_NAME} analysis runner"
    log "Working directory: ${SCRIPT_DIR}"
    
    # Set up trap for cleanup
    trap cleanup EXIT
    
    # Check R installation
    check_r_installation
    
    # Install packages if requested or if running analysis
    if [ "${INSTALL_ONLY}" = "true" ] || [ "${TEST_ONLY}" = "false" ]; then
        install_packages
    fi
    
    # Exit if install-only mode
    if [ "${INSTALL_ONLY}" = "true" ]; then
        log "Package installation completed. Exiting."
        exit 0
    fi
    
    # Run tests if requested or before analysis
    if [ "${TEST_ONLY}" = "true" ]; then
        run_tests
        log "Tests completed. Exiting."
        exit 0
    else
        run_tests
    fi
    
    # Run main analysis
    run_analysis "${FULL_PIPELINE}" "${SKIP_PREPROCESSING}" "${OUTPUT_DIR}"
    
    # Generate summary report
    generate_report
    
    log "${PROJECT_NAME} analysis pipeline completed successfully!"
    log "Results available in: ${OUTPUT_DIR}"
    log "Logs available in: ${LOG_DIR}"
}

# Execute main function
main "$@"
