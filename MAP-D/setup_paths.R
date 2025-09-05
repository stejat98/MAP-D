#!/usr/bin/env Rscript

# Quick setup script to configure your data paths
cat("=== MAP-D Path Configuration ===\n")

# Check what files exist in your system
cat("\nChecking for existing data files...\n")

# Common locations to check
possible_locations <- c(
  "/n/groups/patel/sivateja/UKB",
  "/n/groups/patel/uk_biobank", 
  "/Users/sivatejatang/HMS Dropbox",
  "~/Dropbox",
  getwd()
)

for (loc in possible_locations) {
  if (dir.exists(loc)) {
    cat("✓ Found directory:", loc, "\n")
    
    # Look for key files
    if (dir.exists(file.path(loc, "PEWAS_results"))) {
      cat("  - Has PEWAS_results subdirectory\n")
    }
    
    fst_files <- list.files(loc, pattern = "\\.fst$", recursive = TRUE)
    if (length(fst_files) > 0) {
      cat("  - Found", length(fst_files), ".fst files\n")
      cat("    Example:", head(fst_files, 3), "\n")
    }
    
    rds_files <- list.files(loc, pattern = "\\.rds$|.RDS$", recursive = TRUE)
    if (length(rds_files) > 0) {
      cat("  - Found", length(rds_files), ".rds files\n")
    }
  } else {
    cat("✗ Directory not found:", loc, "\n")
  }
}

# Check for the specific processed data file mentioned in your original scripts
key_files_to_check <- c(
  "/n/groups/patel/sivateja/UKB/all_cause_mort_data_exposures_full_plus_proteomics_plus_GLP_complications.fst",
  "/n/groups/patel/sivateja/UKB/PEWAS_results",
  "/n/groups/patel/uk_biobank/main_data_34521/ukb34521.fst",
  "/n/groups/patel/uk_biobank/olink_22881_52887/olink_data_22881.txt"
)

cat("\nChecking for key data files from your original scripts...\n")
for (file in key_files_to_check) {
  if (file.exists(file)) {
    cat("✓ Found:", file, "\n")
  } else {
    cat("✗ Not found:", file, "\n")
  }
}

# Suggest configuration
cat("\n=== Suggested Configuration ===\n")
cat("Based on your original scripts, update config.R with:\n\n")

cat("# If you have access to the original paths:\n")
cat("BASE_DATA_DIR <- \"/n/groups/patel/sivateja/UKB\"\n")
cat("RAW_DATA_DIR <- \"/n/groups/patel/uk_biobank\"\n")
cat("RESULTS_DIR <- file.path(BASE_DATA_DIR, \"PEWAS_results\")\n\n")

cat("# Or if you need to use local paths:\n")
cat("BASE_DATA_DIR <- \"", getwd(), "/data\"\n", sep = "")
cat("RAW_DATA_DIR <- \"", getwd(), "/data/raw\"\n", sep = "")
cat("RESULTS_DIR <- \"", getwd(), "/results\"\n\n", sep = "")

# Check current config
cat("=== Current Configuration ===\n")
if (file.exists("config.R")) {
  source("config.R")
  cat("BASE_DATA_DIR:", BASE_DATA_DIR, "\n")
  cat("RAW_DATA_DIR:", RAW_DATA_DIR, "\n") 
  cat("RESULTS_DIR:", RESULTS_DIR, "\n")
  cat("PROCESSED_DATA_FILE:", PROCESSED_DATA_FILE, "\n")
  
  # Check if these paths exist
  if (dir.exists(BASE_DATA_DIR)) {
    cat("✓ BASE_DATA_DIR exists\n")
  } else {
    cat("✗ BASE_DATA_DIR does not exist\n")
  }
  
  if (file.exists(PROCESSED_DATA_FILE)) {
    cat("✓ PROCESSED_DATA_FILE exists\n")
  } else {
    cat("✗ PROCESSED_DATA_FILE does not exist\n")
  }
} else {
  cat("config.R not found\n")
}

cat("\n=== Next Steps ===\n")
cat("1. Edit config.R to set the correct paths for your system\n")
cat("2. Run: Rscript test_quick.R (basic validation)\n")
cat("3. Run: Rscript test_with_synthetic_data.R (function testing)\n")
cat("4. If you have real data available, test with a subset\n")
cat("5. Run the full pipeline: source('main_analysis.R'); run_mapd_pipeline()\n")
