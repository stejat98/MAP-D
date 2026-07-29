#!/usr/bin/env Rscript
# Generate REGENIE input files (pheno.txt, covar.txt) for each phenotype
# Enforces EID split: proteins in proteomics-only sample, hallmarks in held-out sample

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# ============================================================================
# Configuration
# ============================================================================

# Input files
MAIN_RDS <- "/n/groups/patel/sivateja/UKB/PEWAS_results/data_plus_GLP_complications_glycemic_status_HbA1c_adjusted.RDS"
PROTEOMICS_EIDS_RDS <- "/n/groups/patel/sivateja/olink_eids_for_proteins_gwas.RDS"
VALIDATED_PROTEINS_CSV <- "/n/groups/patel/sivateja/UKB/merged_validated_proteins_w_EntrezGeneSymbol_w_protein_var_code_UKB.csv"

# Output base directory (overridable for sensitivity analyses with different covariates)
OUTPUT_BASE <- Sys.getenv("REGENIE_INPUT_OUTPUT_BASE",
                          "/n/groups/patel/sivateja/regenie_pipeline/inputs")

# Hallmark traits (non-proteomics / held-out REGENIE phenotypes)
HALLMARK_TRAITS <- c("BMI", "HbA1c", "TRIG_HDL_RATIO")

# Full cohort only
STRATA <- list(
  full = list(filter = NULL, name = "full")
)

# Optional fast path (e.g. single hallmark): set REGENIE_SKIP_PROTEINS=1,
# REGENIE_ONLY_HALLMARKS=BMI
ONLY_STRATUM <- Sys.getenv("REGENIE_ONLY_STRATUM", "")
SKIP_PROTEINS <- Sys.getenv("REGENIE_SKIP_PROTEINS", "") == "1"
ONLY_HALLMARKS <- Sys.getenv("REGENIE_ONLY_HALLMARKS", "")

# ============================================================================
# Helper function: %notin%
# ============================================================================

`%notin%` <- Negate(`%in%`)

# ============================================================================
# Load data
# ============================================================================

cat("Loading main phenotype/covariate dataframe...\n")
main_df <- readRDS(MAIN_RDS)
cat(sprintf("  Loaded %d rows, %d columns\n", nrow(main_df), ncol(main_df)))

cat("Loading proteomics EIDs...\n")
# File was saved with save(), not saveRDS(), so use load()
load(PROTEOMICS_EIDS_RDS)
if (exists("olink_eids_for_proteins_gwas")) {
  olink_eids <- olink_eids_for_proteins_gwas
} else {
  # Try to find the variable
  loaded_vars <- ls()
  eid_var <- grep("eid|EID|olink", loaded_vars, ignore.case = TRUE, value = TRUE)
  if (length(eid_var) > 0) {
    olink_eids <- get(eid_var[1])
  } else {
    stop("ERROR: Variable 'olink_eids_for_proteins_gwas' not found after loading file")
  }
}
cat(sprintf("  Loaded %d proteomics EIDs\n", length(olink_eids)))

cat("Loading validated proteins list...\n")
merged_validated_proteins <- suppressWarnings(fread(VALIDATED_PROTEINS_CSV, fill = TRUE))
# Use Exposure column which contains the actual dataframe column names (e.g., "protein_x.12")
# This matches the column names in the main dataframe, unlike Protein_UKB which has protein names
if ("Exposure" %in% colnames(merged_validated_proteins)) {
  validated_proteins_unique_all_step1_step_2 <- unique(merged_validated_proteins$Exposure)
} else if ("Exposure_UKB" %in% colnames(merged_validated_proteins)) {
  validated_proteins_unique_all_step1_step_2 <- unique(merged_validated_proteins$Exposure_UKB)
} else if ("Protein_UKB" %in% colnames(merged_validated_proteins)) {
  validated_proteins_unique_all_step1_step_2 <- unique(merged_validated_proteins$Protein_UKB)
} else if ("X" %in% colnames(merged_validated_proteins)) {
  validated_proteins_unique_all_step1_step_2 <- unique(merged_validated_proteins$X)
} else {
  stop("ERROR: Neither 'Exposure', 'Exposure_UKB', 'Protein_UKB', nor 'X' column found in validated proteins CSV")
}
# Remove NA and empty values, and ensure unique (already done above, but double-check)
validated_proteins_unique_all_step1_step_2 <- unique(validated_proteins_unique_all_step1_step_2[!is.na(validated_proteins_unique_all_step1_step_2) & validated_proteins_unique_all_step1_step_2 != ""])
cat(sprintf("  Found %d unique validated proteins\n", length(validated_proteins_unique_all_step1_step_2)))

# Convert to data.table for efficiency
if (!is.data.table(main_df)) {
  main_df <- as.data.table(main_df)
}

# Ensure eid column exists
if (!"eid" %in% colnames(main_df)) {
  stop("ERROR: 'eid' column not found in main dataframe")
}

# Ensure GlycemicStatus column exists
if (!"GlycemicStatus" %in% colnames(main_df)) {
  stop("ERROR: 'GlycemicStatus' column not found in main dataframe")
}

# ============================================================================
# Load adjustments (covariates to include)
# ============================================================================

cat("\nLoading adjustments from adjustments_survival_analysis.Rdata...\n")
adjustments_file <- "/home/st320/UKB_PEWAS/adjustments_survival_analysis.Rdata"
if (!file.exists(adjustments_file)) {
  stop(sprintf("ERROR: Adjustments file not found: %s", adjustments_file))
}

# Load the adjustments variable (file was saved with save(), not saveRDS())
load(adjustments_file)
if (!exists("adjustments")) {
  # Try to find the variable
  loaded_vars <- ls()
  adj_var <- grep("adjustment|adj", loaded_vars, ignore.case = TRUE, value = TRUE)
  if (length(adj_var) > 0) {
    adjustments <- get(adj_var[1])
  } else {
    stop("ERROR: Variable 'adjustments' not found after loading file")
  }
}

# Add x.738 to adjustments
adjustments <- c(adjustments, "x.738")
adjustments <- unique(adjustments)  # Remove duplicates if x.738 was already present

cat(sprintf("  Loaded %d covariates from adjustments file\n", length(adjustments) - 1))  # -1 because we added x.738
cat(sprintf("  Added x.738, total: %d covariates\n", length(adjustments)))
cat("  First 10 covariates:", paste(head(adjustments, 10), collapse = ", "), "\n")

# ============================================================================
# Identify covariates (use only specified adjustments)
# ============================================================================

# Filter to only covariates that exist in main_df
covariate_cols <- intersect(adjustments, colnames(main_df))
missing_covars <- setdiff(adjustments, colnames(main_df))
if (length(missing_covars) > 0) {
  cat(sprintf("\nWARNING: %d covariates from adjustments file not found in main_df\n", length(missing_covars)))
  cat("  Missing covariates (first 10):", paste(head(missing_covars, 10), collapse = ", "), "\n")
}

cat(sprintf("\nUsing %d covariates (from adjustments file + x.738)\n", length(covariate_cols)))
cat("  First 10 covariates:", paste(head(covariate_cols, 10), collapse = ", "), "\n")

# ============================================================================
# Preflight checks
# ============================================================================

cat("\n=== PREFLIGHT CHECKS ===\n")

# Check overlap between main_df eids and proteomics eids
overlap_count <- sum(main_df$eid %in% olink_eids)
cat(sprintf("EID overlap: %d individuals in main_df are in proteomics sample\n", overlap_count))
cat(sprintf("Proteomics sample size: %d\n", length(olink_eids)))
cat(sprintf("Held-out sample size (approximate): %d\n", sum(main_df$eid %notin% olink_eids)))

# Check which validated proteins exist in main_df
available_proteins <- intersect(validated_proteins_unique_all_step1_step_2, colnames(main_df))
missing_proteins <- setdiff(validated_proteins_unique_all_step1_step_2, colnames(main_df))
cat(sprintf("\nValidated proteins: %d available, %d missing\n", length(available_proteins), length(missing_proteins)))
if (length(missing_proteins) > 0) {
  cat("  Missing proteins (first 10):", paste(head(missing_proteins, 10), collapse = ", "), "\n")
}

# Check which hallmarks exist
available_hallmarks <- intersect(HALLMARK_TRAITS, colnames(main_df))
missing_hallmarks <- setdiff(HALLMARK_TRAITS, colnames(main_df))
cat(sprintf("\nHallmark traits: %d available, %d missing\n", length(available_hallmarks), length(missing_hallmarks)))
if (length(missing_hallmarks) > 0) {
  cat("  Missing hallmarks:", paste(missing_hallmarks, collapse = ", "), "\n")
}

# ============================================================================
# Function to create REGENIE input files
# ============================================================================

create_regenie_inputs <- function(df_subset, phenotype_name, output_dir, stratum_name, sample_type) {
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Check if phenotype exists
  if (!phenotype_name %in% colnames(df_subset)) {
    return(FALSE)  # Don't print warning for every missing protein
  }

  covar_cols_this <- covariate_cols
  
  # Extract phenotype and eid only first (memory efficient)
  pheno_data <- df_subset[, .(eid, phenotype = get(phenotype_name))]
  
  # Remove rows with missing phenotype
  pheno_data <- pheno_data[!is.na(phenotype)]
  
  if (nrow(pheno_data) == 0) {
    return(FALSE)
  }
  
  # Get eids with non-missing phenotype
  valid_eids <- pheno_data$eid
  
  # Get covariates for individuals with non-missing phenotype (subset first)
  covar_data <- df_subset[eid %in% valid_eids, c("eid", covar_cols_this), with = FALSE]
  
  # Remove covariates with all missing values or constant values
  valid_covar_cols <- c()
  for (col in covar_cols_this) {
    if (col %in% colnames(covar_data)) {
      col_data <- covar_data[[col]]
      if (!all(is.na(col_data)) && length(unique(col_data[!is.na(col_data)])) > 1) {
        valid_covar_cols <- c(valid_covar_cols, col)
      }
    }
  }
  
  # Select valid covariates
  if (length(valid_covar_cols) > 0) {
    covar_data <- covar_data[, c("eid", valid_covar_cols), with = FALSE]
  } else {
    # If no valid covariates, create minimal covariate file with just eid
    covar_data <- covar_data[, .(eid)]
  }
  
  # Merge to ensure same individuals and order
  merged_data <- merge(pheno_data, covar_data, by = "eid", all.x = TRUE)
  
  # Remove duplicate individuals (keep first occurrence)
  # This is critical for REGENIE which requires unique FID/IID pairs
  if (any(duplicated(merged_data$eid))) {
    n_duplicates <- sum(duplicated(merged_data$eid))
    cat(sprintf("    WARNING: Found %d duplicate individuals, keeping first occurrence\n", n_duplicates))
    merged_data <- merged_data[!duplicated(merged_data$eid)]
  }
  
  # Write pheno.txt (REGENIE format: FID IID phenotype)
  pheno_file <- file.path(output_dir, "pheno.txt")
  pheno_output <- merged_data[, .(FID = eid, IID = eid, get("phenotype"))]
  setnames(pheno_output, "V3", phenotype_name)
  fwrite(pheno_output, file = pheno_file, sep = " ", quote = FALSE, na = "NA")
  cat(sprintf("    Wrote %s (%d individuals)\n", pheno_file, nrow(pheno_output)))
  
  # Write covar.txt (REGENIE format: FID IID cov1 cov2 ...)
  covar_file <- file.path(output_dir, "covar.txt")
  covar_cols <- setdiff(colnames(merged_data), c("eid", "phenotype"))
  if (length(covar_cols) > 0) {
    # Create FID and IID directly, then add covariates
    covar_output <- merged_data[, c("eid", covar_cols), with = FALSE]
    covar_output[, FID := eid]
    covar_output[, IID := eid]
    # Reorder columns: FID, IID, then covariates
    setcolorder(covar_output, c("FID", "IID", covar_cols))
    covar_output[, eid := NULL]  # Remove original eid column
  } else {
    # If no covariates, create minimal file with FID and IID only
    covar_output <- merged_data[, .(FID = eid, IID = eid)]
  }
  fwrite(covar_output, file = covar_file, sep = " ", quote = FALSE, na = "NA")
  cat(sprintf("    Wrote %s (%d covariates)\n", covar_file, max(0, ncol(covar_output) - 2)))
  
  return(TRUE)
}

# ============================================================================
# Generate inputs for each stratum and sample type
# ============================================================================

cat("\n=== GENERATING INPUT FILES ===\n")

for (stratum_name in names(STRATA)) {
  if (nzchar(ONLY_STRATUM) && stratum_name != ONLY_STRATUM) {
    next
  }
  stratum_info <- STRATA[[stratum_name]]
  cat(sprintf("\n--- Processing stratum: %s ---\n", stratum_name))
  
  # Filter by stratum if needed (use data.table syntax for memory efficiency)
  if (!is.null(stratum_info$filter)) {
    df_stratum <- main_df[GlycemicStatus == stratum_info$filter]
  } else {
    df_stratum <- main_df  # Don't copy, just reference
  }
  cat(sprintf("  Stratum size: %d individuals\n", nrow(df_stratum)))
  
  # Split by EID (use data.table syntax, avoid copying when possible)
  df_proteins_only <- df_stratum[eid %in% olink_eids]
  df_hallmarks_heldout <- df_stratum[eid %notin% olink_eids]
  
  # Clean up memory
  if (!is.null(stratum_info$filter)) {
    rm(df_stratum)
    gc()
  }
  
  cat(sprintf("  Proteins-only sample: %d individuals\n", nrow(df_proteins_only)))
  cat(sprintf("  Hallmarks-heldout sample: %d individuals\n", nrow(df_hallmarks_heldout)))
  
  # Process proteins in proteins-only sample
  proteins_output_dir <- file.path(OUTPUT_BASE, stratum_name, "proteins_only")
  if (isTRUE(SKIP_PROTEINS)) {
    cat("\n  Skipping proteins (REGENIE_SKIP_PROTEINS=1)\n")
    proteins_processed <- 0L
  } else {
    cat("\n  Processing proteins in proteins-only sample...\n")
    proteins_processed <- 0
    for (prot in available_proteins) {
      prot_dir <- file.path(proteins_output_dir, prot)
      if (create_regenie_inputs(df_proteins_only, prot, prot_dir, stratum_name, "proteins_only")) {
        proteins_processed <- proteins_processed + 1
      }
      if (proteins_processed %% 10 == 0) {
        gc(verbose = FALSE)
      }
    }
    cat(sprintf("  Processed %d/%d proteins\n", proteins_processed, length(available_proteins)))
  }
  
  # Process hallmarks in held-out sample
  cat("\n  Processing hallmarks in held-out sample...\n")
  hallmarks_output_dir <- file.path(OUTPUT_BASE, stratum_name, "hallmarks_heldout")
  
  hallmarks_run <- if (nzchar(ONLY_HALLMARKS)) {
    intersect(available_hallmarks, strsplit(ONLY_HALLMARKS, ",", fixed = TRUE)[[1]])
  } else {
    available_hallmarks
  }
  if (length(hallmarks_run) == 0L) {
    cat("\n  No hallmarks to run (check REGENIE_ONLY_HALLMARKS vs available names)\n")
  }

  hallmarks_processed <- 0
  for (hallmark in hallmarks_run) {
    hallmark_dir <- file.path(hallmarks_output_dir, hallmark)
    if (create_regenie_inputs(df_hallmarks_heldout, hallmark, hallmark_dir, stratum_name, "hallmarks_heldout")) {
      hallmarks_processed <- hallmarks_processed + 1
    }
    # Force garbage collection after each hallmark to free memory
    gc(verbose = FALSE)
  }
  cat(sprintf("  Processed %d/%d hallmarks\n", hallmarks_processed, length(hallmarks_run)))
  
  # Clean up large dataframes before next stratum
  rm(df_proteins_only, df_hallmarks_heldout)
  gc(verbose = FALSE)
}

cat("\n=== COMPLETE ===\n")
cat("All input files generated successfully!\n")
