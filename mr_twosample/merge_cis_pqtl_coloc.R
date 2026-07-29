#!/usr/bin/env Rscript
#
# Merge Cis-pQTL MR results with colocalization data
#
# Merges:
# 1. Table_MR_Cis_pQTL_HbA1c_Full.csv with merged_step1_2_ukb_T2D_coloc.csv
# 2. Table_MR_Cis_pQTL_BMI_Full.csv with merged_step1_2_ukb_obesity_coloc.csv
#
# Removes colocalization-specific columns before merging and relabels appropriately.
#
# Usage: module load gcc/14.2.0 R/4.4.2 && Rscript merge_cis_pqtl_coloc.R

# Define base paths
base_path <- "/n/groups/patel/sivateja/regenie_pipeline"
mr_path <- file.path(base_path, "results/twosampleMR/supplemental_tables")
output_path <- file.path(base_path, "results/twosampleMR/supplemental_tables")

# Columns to remove from coloc files
coloc_cols_to_remove <- c(
  "nsnps",
  "PP.H0.abf",
  "PP.H1.abf",
  "PP.H2.abf",
  "PP.H3.abf",
  "PP.H4.abf"
)

#' Load and clean colocalization data
#'
#' @param filepath Path to the colocalization CSV file
#' @param coloc_sig_col Name of the coloc significance column
#' @return Cleaned dataframe with coloc-specific columns removed
load_and_clean_coloc <- function(filepath, coloc_sig_col) {
  message("Loading: ", basename(filepath))
  
  df <- read.csv(filepath, stringsAsFactors = FALSE, check.names = FALSE)
  
  # Columns to remove (including file-specific sig column)
  cols_to_remove <- c(coloc_cols_to_remove, coloc_sig_col)
  
  # Also remove unnamed index columns (X, empty string, etc.)
  index_cols <- names(df)[grepl("^X$|^$|^Unnamed", names(df))]
  cols_to_remove <- c(cols_to_remove, index_cols)
  
  # Get columns that actually exist
  existing_cols_to_remove <- intersect(cols_to_remove, names(df))
  
  # Remove columns
  df_cleaned <- df[, !(names(df) %in% existing_cols_to_remove), drop = FALSE]
  
  message("  Original columns: ", ncol(df))
  message("  Removed columns: ", paste(existing_cols_to_remove, collapse = ", "))
  message("  Remaining columns: ", ncol(df_cleaned))
  
  return(df_cleaned)
}

#' Load MR Cis-pQTL data
#'
#' @param filepath Path to the MR CSV file
#' @return MR dataframe
load_mr_data <- function(filepath) {
  message("Loading: ", basename(filepath))
  
  df <- read.csv(filepath, stringsAsFactors = FALSE, check.names = FALSE)
  
  message("  Columns: ", paste(names(df), collapse = ", "))
  message("  Rows: ", nrow(df))
  
  return(df)
}

#' Rename MR columns with Cis_pQTL prefix
#'
#' @param df MR dataframe
#' @param merge_key Column name to keep unchanged
#' @return Dataframe with renamed columns
rename_mr_columns <- function(df, merge_key = "Protein") {
  new_names <- names(df)
  
  for (i in seq_along(new_names)) {
    col <- new_names[i]
    if (col != merge_key) {
      # Clean column name and add prefix
      clean_name <- gsub(" ", "_", col)
      clean_name <- gsub("<", "lt", clean_name)
      clean_name <- gsub("-", "_", clean_name)
      new_names[i] <- paste0("Cis_pQTL_", clean_name)
    }
  }
  
  names(df) <- new_names
  return(df)
}

#' Merge MR and colocalization datasets
#'
#' @param mr_df MR Cis-pQTL data
#' @param coloc_df Cleaned colocalization data
#' @param mr_key Column name in MR data to merge on
#' @param coloc_key Column name in coloc data to merge on
#' @param description Description for logging
#' @return Merged dataframe
merge_datasets <- function(mr_df, coloc_df, mr_key = "Protein", 
                           coloc_key = "code_UKB", description = "") {
  
  # Rename MR columns with Cis_pQTL prefix
  mr_renamed <- rename_mr_columns(mr_df, merge_key = mr_key)
  
  # Perform left join using base R merge (keep all MR results)
  merged <- merge(mr_renamed, coloc_df, 
                  by.x = mr_key, by.y = coloc_key, 
                  all.x = TRUE, sort = FALSE)
  
  # Count matches (non-NA in first coloc column)
  coloc_cols <- setdiff(names(coloc_df), coloc_key)
  if (length(coloc_cols) > 0) {
    matched <- sum(!is.na(merged[[coloc_cols[1]]]))
  } else {
    matched <- 0
  }
  
  message("\nMerge results (", description, "):")
  message("  MR records: ", nrow(mr_df))
  message("  Coloc records: ", nrow(coloc_df))
  message("  Merged records: ", nrow(merged))
  message("  Matched records: ", matched)
  
  return(merged)
}

# ==============================================================================
# Main execution
# ==============================================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("Processing HbA1c/T2D merge")
message(paste(rep("=", 60), collapse = ""))

# Load HbA1c MR data
mr_hba1c <- load_mr_data(
  file.path(mr_path, "Table_MR_Cis_pQTL_HbA1c_Full.csv")
)

# Load and clean T2D coloc data
coloc_t2d <- load_and_clean_coloc(
  file.path(base_path, "merged_step1_2_ukb_T2D_coloc.csv"),
  coloc_sig_col = "T2D_coloc_sig"
)

# Merge HbA1c with T2D coloc
merged_hba1c <- merge_datasets(
  mr_hba1c, coloc_t2d,
  mr_key = "Protein",
  coloc_key = "code_UKB",
  description = "HbA1c/T2D"
)

# Save HbA1c merged data
output_hba1c <- file.path(output_path, "Table_MR_Cis_pQTL_HbA1c_Full_merged_coloc.csv")
write.csv(merged_hba1c, output_hba1c, row.names = FALSE)
message("\nSaved to: ", output_hba1c)
message("Final columns: ", paste(names(merged_hba1c), collapse = ", "))

message("\n", paste(rep("=", 60), collapse = ""))
message("Processing BMI/Obesity merge")
message(paste(rep("=", 60), collapse = ""))

# Load BMI MR data
mr_bmi <- load_mr_data(
  file.path(mr_path, "Table_MR_Cis_pQTL_BMI_Full.csv")
)

# Load and clean Obesity coloc data
coloc_obesity <- load_and_clean_coloc(
  file.path(base_path, "merged_step1_2_ukb_obesity_coloc.csv"),
  coloc_sig_col = "coloc_sig"
)

# Merge BMI with Obesity coloc
merged_bmi <- merge_datasets(
  mr_bmi, coloc_obesity,
  mr_key = "Protein",
  coloc_key = "code_UKB",
  description = "BMI/Obesity"
)

# Save BMI merged data
output_bmi <- file.path(output_path, "Table_MR_Cis_pQTL_BMI_Full_merged_coloc.csv")
write.csv(merged_bmi, output_bmi, row.names = FALSE)
message("\nSaved to: ", output_bmi)
message("Final columns: ", paste(names(merged_bmi), collapse = ", "))

message("\n", paste(rep("=", 60), collapse = ""))
message("Merge operations completed successfully!")
message(paste(rep("=", 60), collapse = ""))
