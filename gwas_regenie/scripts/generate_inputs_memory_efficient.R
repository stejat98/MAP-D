#!/usr/bin/env Rscript
# Memory-efficient version: Process one stratum at a time, clear memory between

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# Configuration
MAIN_RDS <- "/path/to/project/UKB/PEWAS_results/data_plus_GLP_complications_glycemic_status_HbA1c_adjusted.RDS"
PROTEOMICS_EIDS_RDS <- "/path/to/project/olink_eids_for_proteins_gwas.RDS"
VALIDATED_PROTEINS_CSV <- "/path/to/project/UKB/merged_validated_proteins_w_EntrezGeneSymbol_w_protein_var_code_UKB.csv"
OUTPUT_BASE <- "/path/to/project/regenie_pipeline/inputs"

HALLMARK_TRAITS <- c("BMI", "HbA1c", "TRIG_HDL_RATIO")

STRATA <- list(
  full = list(filter = NULL, name = "full")
)

`%notin%` <- Negate(`%in%`)

# Load only what we need
cat("Loading data...\n")
main_df <- readRDS(MAIN_RDS)
if (!is.data.table(main_df)) {
  main_df <- as.data.table(main_df)
}

load(PROTEOMICS_EIDS_RDS)
if (exists("olink_eids_for_proteins_gwas")) {
  olink_eids <- olink_eids_for_proteins_gwas
} else {
  stop("ERROR: Variable 'olink_eids_for_proteins_gwas' not found")
}

merged_validated_proteins <- suppressWarnings(fread(VALIDATED_PROTEINS_CSV, fill = TRUE))
# Use Exposure column which contains the actual dataframe column names (e.g., "protein_x.12")
# This matches the column names in the main dataframe, unlike Protein_UKB which has protein names
if ("Exposure" %in% colnames(merged_validated_proteins)) {
  validated_proteins <- unique(merged_validated_proteins$Exposure)
} else if ("Exposure_UKB" %in% colnames(merged_validated_proteins)) {
  validated_proteins <- unique(merged_validated_proteins$Exposure_UKB)
} else if ("Protein_UKB" %in% colnames(merged_validated_proteins)) {
  validated_proteins <- unique(merged_validated_proteins$Protein_UKB)
} else {
  stop("ERROR: No protein column found (Exposure, Exposure_UKB, or Protein_UKB)")
}
# Remove NA and empty values, and ensure unique (already done above, but double-check)
validated_proteins <- unique(validated_proteins[!is.na(validated_proteins) & validated_proteins != ""])

# Identify available proteins and hallmarks
available_proteins <- intersect(validated_proteins, colnames(main_df))
available_hallmarks <- intersect(HALLMARK_TRAITS, colnames(main_df))

cat(sprintf("Available: %d proteins, %d hallmarks\n", length(available_proteins), length(available_hallmarks)))

# Identify covariates (do this once)
exclude_cols <- c("eid", "GlycemicStatus", HALLMARK_TRAITS, validated_proteins)
covariate_cols <- setdiff(colnames(main_df), exclude_cols)

# Function to create files
create_regenie_inputs <- function(df_subset, phenotype_name, output_dir) {
  if (!phenotype_name %in% colnames(df_subset)) {
    return(FALSE)
  }
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Get phenotype and eid only
  pheno_data <- df_subset[, .(eid, phenotype = get(phenotype_name))]
  pheno_data <- pheno_data[!is.na(phenotype)]
  
  if (nrow(pheno_data) == 0) {
    return(FALSE)
  }
  
  valid_eids <- pheno_data$eid

  covar_cols_this <- covariate_cols
  
  # Get covariates for valid individuals only
  covar_data <- df_subset[eid %in% valid_eids, c("eid", covar_cols_this), with = FALSE]
  
  # Remove constant/missing covariates
  valid_covar_cols <- c()
  for (col in covar_cols_this) {
    if (col %in% colnames(covar_data)) {
      col_data <- covar_data[[col]]
      if (!all(is.na(col_data)) && length(unique(col_data[!is.na(col_data)])) > 1) {
        valid_covar_cols <- c(valid_covar_cols, col)
      }
    }
  }
  
  if (length(valid_covar_cols) > 0) {
    covar_data <- covar_data[, c("eid", valid_covar_cols), with = FALSE]
  } else {
    covar_data <- covar_data[, .(eid)]
  }
  
  # Merge
  merged_data <- merge(pheno_data, covar_data, by = "eid", all.x = TRUE)
  
  # Write pheno.txt
  pheno_file <- file.path(output_dir, "pheno.txt")
  pheno_output <- merged_data[, .(FID = eid, IID = eid, get("phenotype"))]
  setnames(pheno_output, "V3", phenotype_name)
  fwrite(pheno_output, file = pheno_file, sep = " ", quote = FALSE, na = "NA")
  
  # Write covar.txt
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
    covar_output <- merged_data[, .(FID = eid, IID = eid)]
  }
  fwrite(covar_output, file = covar_file, sep = " ", quote = FALSE, na = "NA")
  
  # Clean up
  rm(pheno_data, covar_data, merged_data, pheno_output, covar_output)
  gc(verbose = FALSE)
  
  return(TRUE)
}

# Process one stratum at a time
cat("\n=== GENERATING INPUT FILES ===\n")

for (stratum_name in names(STRATA)) {
  stratum_info <- STRATA[[stratum_name]]
  cat(sprintf("\n--- Processing stratum: %s ---\n", stratum_name))
  
  # Filter stratum
  if (!is.null(stratum_info$filter)) {
    df_stratum <- main_df[GlycemicStatus == stratum_info$filter]
  } else {
    df_stratum <- main_df
  }
  
  # Split by EID
  df_proteins_only <- df_stratum[eid %in% olink_eids]
  df_hallmarks_heldout <- df_stratum[eid %notin% olink_eids]
  
  cat(sprintf("  Proteins-only: %d, Held-out: %d\n", nrow(df_proteins_only), nrow(df_hallmarks_heldout)))
  
  # Process proteins
  cat("  Processing proteins...\n")
  proteins_output_dir <- file.path(OUTPUT_BASE, stratum_name, "proteins_only")
  proteins_processed <- 0
  
  for (i in seq_along(available_proteins)) {
    prot <- available_proteins[i]
    prot_dir <- file.path(proteins_output_dir, prot)
    if (create_regenie_inputs(df_proteins_only, prot, prot_dir)) {
      proteins_processed <- proteins_processed + 1
    }
    
    # Aggressive GC every 5 proteins
    if (i %% 5 == 0) {
      gc(verbose = FALSE)
    }
  }
  cat(sprintf("  Processed %d/%d proteins\n", proteins_processed, length(available_proteins)))
  
  # Clean up proteins dataframe
  rm(df_proteins_only)
  gc(verbose = FALSE)
  
  # Process hallmarks
  cat("  Processing hallmarks...\n")
  hallmarks_output_dir <- file.path(OUTPUT_BASE, stratum_name, "hallmarks_heldout")
  hallmarks_processed <- 0
  
  for (hallmark in available_hallmarks) {
    hallmark_dir <- file.path(hallmarks_output_dir, hallmark)
    if (create_regenie_inputs(df_hallmarks_heldout, hallmark, hallmark_dir)) {
      hallmarks_processed <- hallmarks_processed + 1
    }
    gc(verbose = FALSE)
  }
  cat(sprintf("  Processed %d/%d hallmarks\n", hallmarks_processed, length(available_hallmarks)))
  
  # Clean up before next stratum
  rm(df_hallmarks_heldout)
  if (!is.null(stratum_info$filter)) {
    rm(df_stratum)
  }
  gc(verbose = FALSE)
  
  cat(sprintf("  Completed stratum: %s\n", stratum_name))
}

cat("\n=== COMPLETE ===\n")
