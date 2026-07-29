#!/usr/bin/env Rscript
# Preflight checks and sanity validations for REGENIE pipeline

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# Check if R.utils is needed for gzipped CSV
if (!requireNamespace("R.utils", quietly = TRUE)) {
  cat("WARNING: R.utils package not found. Attempting to install or use alternative method.\n")
}

# ============================================================================
# Configuration
# ============================================================================

MAIN_RDS <- "/path/to/project/UKB/PEWAS_results/data_plus_GLP_complications_glycemic_status_HbA1c_adjusted.RDS"
PROTEOMICS_EIDS_RDS <- "/path/to/project/olink_eids_for_proteins_gwas.RDS"
VALIDATED_PROTEINS_CSV <- "/path/to/project/UKB/merged_validated_proteins_w_EntrezGeneSymbol_w_protein_var_code_UKB.csv"
INPUT_BASE <- "/path/to/project/regenie_pipeline/inputs"

HALLMARK_TRAITS <- c("BMI", "HbA1c", "TRIG_HDL_RATIO")

# ============================================================================
# Helper function
# ============================================================================

`%notin%` <- Negate(`%in%`)

# ============================================================================
# Load and validate data
# ============================================================================

cat("=== PREFLIGHT CHECKS ===\n\n")

# 1. Check input files exist
cat("1. Checking input files...\n")
files_to_check <- list(
  "Main RDS" = MAIN_RDS,
  "Proteomics EIDs RDS" = PROTEOMICS_EIDS_RDS,
  "Validated proteins CSV" = VALIDATED_PROTEINS_CSV
)

all_files_exist <- TRUE
for (name in names(files_to_check)) {
  if (file.exists(files_to_check[[name]])) {
    cat(sprintf("   ✓ %s exists\n", name))
  } else {
    cat(sprintf("   ✗ %s NOT FOUND: %s\n", name, files_to_check[[name]]))
    all_files_exist <- FALSE
  }
}

if (!all_files_exist) {
  stop("ERROR: Some input files are missing")
}

# 2. Load and check data
cat("\n2. Loading and validating data...\n")

main_df <- readRDS(MAIN_RDS)
cat(sprintf("   Main dataframe: %d rows, %d columns\n", nrow(main_df), ncol(main_df)))

# Try to load proteomics EIDs (saved with save(), not saveRDS())
cat("   Loading proteomics EIDs (using load())...\n")
tryCatch({
  # Load the .RDS file (actually .RData format saved with save())
  load(PROTEOMICS_EIDS_RDS)
  
  # Check if the expected variable name exists
  if (exists("olink_eids_for_proteins_gwas")) {
    olink_eids <- olink_eids_for_proteins_gwas
    cat(sprintf("   ✓ Proteomics EIDs loaded: %d EIDs (class: %s)\n", length(olink_eids), class(olink_eids)[1]))
  } else {
    # Check what variables were loaded
    loaded_vars <- ls()
    cat("   ⚠ WARNING: Variable 'olink_eids_for_proteins_gwas' not found.\n")
    cat("   Variables loaded from file:", paste(loaded_vars, collapse=", "), "\n")
    cat("   Attempting to find EID vector...\n")
    
    # Try to find any vector variable
    for (var_name in loaded_vars) {
      var_obj <- get(var_name)
      if (is.vector(var_obj) || is.numeric(var_obj) || is.character(var_obj)) {
        olink_eids <<- var_obj
        cat(sprintf("   ✓ Using variable '%s': %d EIDs\n", var_name, length(olink_eids)))
        break
      }
    }
    
    if (!exists("olink_eids")) {
      olink_eids <<- NULL
    }
  }
}, error = function(e) {
  cat(sprintf("   ✗ ERROR loading proteomics EIDs: %s\n", e$message))
  olink_eids <<- NULL
})

# Try to read validated proteins CSV (handle gzipped files and encoding issues)
cat("   Reading validated proteins CSV...\n")
tryCatch({
  # Suppress warnings about column mismatches (common with messy CSVs)
  merged_validated_proteins <- suppressWarnings(fread(VALIDATED_PROTEINS_CSV, fill = TRUE))
  
  # Use Exposure column which contains the actual dataframe column names (e.g., "protein_x.12")
  # This matches the column names in the main dataframe, unlike Protein_UKB which has protein names
  if ("Exposure" %in% colnames(merged_validated_proteins)) {
    validated_proteins <- unique(merged_validated_proteins$Exposure)
    validated_proteins <- unique(validated_proteins[!is.na(validated_proteins) & validated_proteins != ""])
    cat(sprintf("   ✓ Validated proteins (from Exposure): %d\n", length(validated_proteins)))
  } else if ("Exposure_UKB" %in% colnames(merged_validated_proteins)) {
    validated_proteins <- unique(merged_validated_proteins$Exposure_UKB)
    validated_proteins <- unique(validated_proteins[!is.na(validated_proteins) & validated_proteins != ""])
    cat(sprintf("   ✓ Validated proteins (from Exposure_UKB): %d\n", length(validated_proteins)))
  } else if ("Protein_UKB" %in% colnames(merged_validated_proteins)) {
    validated_proteins <- unique(merged_validated_proteins$Protein_UKB)
    validated_proteins <- unique(validated_proteins[!is.na(validated_proteins) & validated_proteins != ""])
    cat(sprintf("   ✓ Validated proteins (from Protein_UKB): %d\n", length(validated_proteins)))
  } else if ("X" %in% colnames(merged_validated_proteins)) {
    # CSV has single column "X" with protein names
    validated_proteins <- unique(merged_validated_proteins$X)
    validated_proteins <- unique(validated_proteins[!is.na(validated_proteins) & validated_proteins != ""])
    cat(sprintf("   ✓ Validated proteins (from column 'X'): %d\n", length(validated_proteins)))
  } else {
    cat("   ⚠ WARNING: 'Exposure', 'Exposure_UKB', 'Protein_UKB', or 'X' column not found. Available columns:\n")
    cat(sprintf("     %s\n", paste(head(colnames(merged_validated_proteins), 10), collapse=", ")))
    validated_proteins <<- character(0)
  }
}, error = function(e) {
  cat(sprintf("   ✗ ERROR reading validated proteins CSV: %s\n", e$message))
  if (grepl("R.utils", e$message)) {
    cat("   Suggestion: Install R.utils package:\n")
    cat("     install.packages('R.utils')\n")
  }
  validated_proteins <<- character(0)
})

# 3. Check required columns
cat("\n3. Checking required columns...\n")
required_cols <- c("eid", "GlycemicStatus")
for (col in required_cols) {
  if (col %in% colnames(main_df)) {
    cat(sprintf("   ✓ Column '%s' exists\n", col))
  } else {
    stop(sprintf("ERROR: Required column '%s' not found", col))
  }
}

# 4. Check EID overlap
cat("\n4. Checking EID overlap...\n")
if (exists("olink_eids") && !is.null(olink_eids)) {
  overlap <- sum(main_df$eid %in% olink_eids)
  heldout <- sum(main_df$eid %notin% olink_eids)
  cat(sprintf("   Proteomics sample: %d individuals\n", overlap))
  cat(sprintf("   Held-out sample: %d individuals\n", heldout))
  cat(sprintf("   Total in main_df: %d individuals\n", nrow(main_df)))
  
  if (overlap < 10000) {
    cat("   ⚠ WARNING: Proteomics sample seems small (<10k)\n")
  }
  if (heldout < 100000) {
    cat("   ⚠ WARNING: Held-out sample seems small (<100k)\n")
  }
} else {
  cat("   ⚠ SKIPPED: Cannot check EID overlap (proteomics EIDs file not loaded)\n")
  cat("   Please check the file format (should be saved with save(), loaded with load())\n")
}
if (heldout < 100000) {
  cat("   WARNING: Held-out sample seems small (<100k)\n")
}

# 5. Check phenotype availability
cat("\n5. Checking phenotype availability...\n")

if (exists("validated_proteins") && length(validated_proteins) > 0) {
  available_proteins <- intersect(validated_proteins, colnames(main_df))
  missing_proteins <- setdiff(validated_proteins, colnames(main_df))
  cat(sprintf("   Validated proteins: %d available, %d missing\n", 
              length(available_proteins), length(missing_proteins)))
  if (length(missing_proteins) > 0 && length(missing_proteins) <= 20) {
    cat("   Missing proteins:", paste(head(missing_proteins, 10), collapse=", "), "\n")
  }
} else {
  cat("   ⚠ SKIPPED: Cannot check validated proteins (file not loaded)\n")
  available_proteins <- character(0)
}

available_hallmarks <- intersect(HALLMARK_TRAITS, colnames(main_df))
missing_hallmarks <- setdiff(HALLMARK_TRAITS, colnames(main_df))
cat(sprintf("   Hallmark traits: %d available, %d missing\n", 
            length(available_hallmarks), length(missing_hallmarks)))

if (length(missing_hallmarks) > 0) {
  cat("   Missing hallmarks:", paste(missing_hallmarks, collapse = ", "), "\n")
}

# 6. Check GlycemicStatus values
cat("\n6. Checking GlycemicStatus distribution...\n")
if ("GlycemicStatus" %in% colnames(main_df)) {
  status_table <- table(main_df$GlycemicStatus, useNA = "ifany")
  for (status in names(status_table)) {
    cat(sprintf("   %s: %d individuals\n", status, status_table[status]))
  }
}

# 7. Check generated input files
cat("\n7. Checking generated input files...\n")
STRATA <- c("full")
SAMPLE_TYPES <- c("proteins_only", "hallmarks_heldout")

total_phenos <- 0
for (stratum in STRATA) {
  for (sample_type in SAMPLE_TYPES) {
    input_dir <- file.path(INPUT_BASE, stratum, sample_type)
    if (dir.exists(input_dir)) {
      pheno_dirs <- list.dirs(input_dir, recursive = FALSE)
      pheno_count <- 0
      for (pheno_dir in pheno_dirs) {
        pheno_name <- basename(pheno_dir)
        pheno_file <- file.path(pheno_dir, "pheno.txt")
        covar_file <- file.path(pheno_dir, "covar.txt")
        if (file.exists(pheno_file) && file.exists(covar_file)) {
          pheno_count <- pheno_count + 1
          total_phenos <- total_phenos + 1
        }
      }
      cat(sprintf("   %s/%s: %d phenotypes\n", stratum, sample_type, pheno_count))
    } else {
      cat(sprintf("   %s/%s: directory not found (run generate_inputs.R first)\n", 
                  stratum, sample_type))
    }
  }
}

cat(sprintf("\n   Total phenotypes with input files: %d\n", total_phenos))

# 8. Check genotype files
cat("\n8. Checking genotype files...\n")
GWAS_PATH <- "/path/to/shared_data/UKB/gwas"
STEP1_PGEN <- file.path(GWAS_PATH, "ukb_nonimputed_snps.pgen")
STEP2_PGEN <- file.path(GWAS_PATH, "UKBallchr.pgen")

if (file.exists(STEP1_PGEN) || file.exists(paste0(STEP1_PGEN, ".zst"))) {
  cat("   ✓ Step 1 pgen file exists\n")
} else {
  cat("   ✗ Step 1 pgen file not found (check path)\n")
}

if (file.exists(STEP2_PGEN) || file.exists(paste0(STEP2_PGEN, ".zst"))) {
  cat("   ✓ Step 2 pgen file exists\n")
} else {
  cat("   ✗ Step 2 pgen file not found (check path)\n")
}

cat("\n=== PREFLIGHT CHECKS COMPLETE ===\n")
cat("If all checks passed, you can proceed with:\n")
cat("  1. Run generate_inputs.R to create input files\n")
cat("  2. Run submit_regenie.sh to submit GWAS jobs\n")
