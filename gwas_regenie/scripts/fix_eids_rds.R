#!/usr/bin/env Rscript
# Helper script to fix/re-save the proteomics EIDs RDS file
# This may be needed if the file was created with a different R version

INPUT_FILE <- "/n/groups/patel/sivateja/olink_eids_for_proteins_gwas.RDS"
OUTPUT_FILE <- "/n/groups/patel/sivateja/olink_eids_for_proteins_gwas_fixed.RDS"

cat("Attempting to fix RDS file...\n")
cat("Input:", INPUT_FILE, "\n")

# Try different R versions or methods
cat("\nMethod 1: Try readRDS with current R version...\n")
tryCatch({
  eids <- readRDS(INPUT_FILE)
  cat(sprintf("  ✓ Success! Loaded %d EIDs\n", length(eids)))
  cat(sprintf("  Class: %s\n", paste(class(eids), collapse=", ")))
  
  # Re-save with version 2 (more compatible)
  saveRDS(eids, OUTPUT_FILE, version = 2)
  cat(sprintf("  ✓ Re-saved to: %s\n", OUTPUT_FILE))
  cat("\nYou can now use the fixed file or replace the original:\n")
  cat(sprintf("  mv %s %s\n", OUTPUT_FILE, INPUT_FILE))
  
}, error = function(e) {
  cat(sprintf("  ✗ Failed: %s\n", e$message))
  cat("\nThe file may need to be loaded with the R version that created it.\n")
  cat("Try:\n")
  cat("  1. Load the R version that created the file\n")
  cat("  2. Run: eids <- readRDS('", INPUT_FILE, "')\n", sep="")
  cat("  3. Run: saveRDS(eids, '", INPUT_FILE, "', version=2)\n", sep="")
})
