#!/usr/bin/env Rscript
# MR Summary Statistics Table for BMI and HbA1c
# Creates table with Overall, Cis-only, and Trans-only QTL results
# Using base R only

# Define paths
results_dir <- "/path/to/project/regenie_pipeline/results/twosampleMR"
cis_trans_dir <- file.path(results_dir, "cis_trans_chunked")

cat("Reading cis/trans MR results...\n")

# Read and combine cis chunks
cis_files <- list.files(cis_trans_dir, pattern = "MR_cis_chunk.*\\.csv", full.names = TRUE)
cis_list <- lapply(cis_files, read.csv)
cis_df <- do.call(rbind, cis_list)

# Read and combine trans chunks
trans_files <- list.files(cis_trans_dir, pattern = "MR_trans_chunk.*\\.csv", full.names = TRUE)
trans_list <- lapply(trans_files, read.csv)
trans_df <- do.call(rbind, trans_list)

# Filter for BMI and HbA1c only
cis_df <- cis_df[cis_df$phenotype %in% c("BMI", "HbA1c"), ]
trans_df <- trans_df[trans_df$phenotype %in% c("BMI", "HbA1c"), ]

cat("=== Data shapes ===\n")
cat(sprintf("Cis MR (BMI/HbA1c): %d rows\n", nrow(cis_df)))
cat(sprintf("Trans MR (BMI/HbA1c): %d rows\n", nrow(trans_df)))
cat("\n")

# Overall MR files
overall_files <- list(
  BMI_full = "MR_proteins_BMI_full.csv",
  HbA1c_full = "MR_proteins_HbA1c_full.csv",
  TRIG_HDL_RATIO_full = "MR_proteins_TRIG_HDL_RATIO_full.csv"
)

# Function to get summary stats for cis/trans
# For cis (single SNP): use "Wald ratio"
# For trans (multiple SNPs): use "Inverse variance weighted"
get_summary <- function(df, pheno, strat, qtl_type = NULL) {
  subset <- df[df$phenotype == pheno & df$stratum == strat, ]
  
  # Determine method based on qtl_type
  if (!is.null(qtl_type) && "qtl_type" %in% names(subset)) {
    subset <- subset[subset$qtl_type == qtl_type, ]
  }
  
  # For cis, we typically have Wald ratio (single SNP)
  # For trans, we use IVW when nsnp > 1
  if ("method" %in% names(subset)) {
    # Keep IVW for multi-SNP or Wald ratio for single-SNP
    subset <- subset[subset$method %in% c("Inverse variance weighted", "Wald ratio"), ]
    # Prefer IVW if both exist for same protein, else keep Wald ratio
    # For simplicity, just take first method per protein
    if (nrow(subset) > 0) {
      subset <- do.call(rbind, by(subset, subset$protein, function(x) {
        if (any(x$method == "Inverse variance weighted")) {
          x[x$method == "Inverse variance weighted", ][1, ]
        } else {
          x[x$method == "Wald ratio", ][1, ]
        }
      }))
    }
  }
  
  n_proteins <- length(unique(subset$protein))
  n_sig <- sum(subset$pval < 0.05, na.rm = TRUE)
  n_bonf <- ifelse(n_proteins > 0, sum(subset$pval < 0.05/n_proteins, na.rm = TRUE), 0)
  
  return(list(n_proteins = n_proteins, n_sig_p05 = n_sig, n_bonf = n_bonf))
}

# Create summary table
summary_data <- data.frame()

for (phenotype in c("BMI", "HbA1c")) {
  for (stratum in c("full")) {
    # Overall (from original MR files)
    overall_key <- paste0(phenotype, "_", stratum)
    overall_file <- file.path(results_dir, overall_files[[overall_key]])
    
    if (file.exists(overall_file)) {
      overall_df <- read.csv(overall_file)
      ivw_df <- overall_df[overall_df$method == "Inverse variance weighted", ]
      n_proteins_overall <- length(unique(ivw_df$exposure))
      n_sig_overall <- sum(ivw_df$pval < 0.05, na.rm = TRUE)
      n_bonf_overall <- ifelse(n_proteins_overall > 0, sum(ivw_df$pval < 0.05/n_proteins_overall, na.rm = TRUE), 0)
    } else {
      n_proteins_overall <- n_sig_overall <- n_bonf_overall <- 0
    }
    
    # Cis-only
    cis_stats <- get_summary(cis_df, phenotype, stratum, qtl_type = "cis")
    
    # Trans-only  
    trans_stats <- get_summary(trans_df, phenotype, stratum, qtl_type = "trans")
    
    row <- data.frame(
      Phenotype = phenotype,
      Stratum = stratum,
      Overall_N_Proteins = n_proteins_overall,
      Overall_N_Sig_p05 = n_sig_overall,
      Overall_N_Bonf = n_bonf_overall,
      Cis_N_Proteins = cis_stats$n_proteins,
      Cis_N_Sig_p05 = cis_stats$n_sig_p05,
      Cis_N_Bonf = cis_stats$n_bonf,
      Trans_N_Proteins = trans_stats$n_proteins,
      Trans_N_Sig_p05 = trans_stats$n_sig_p05,
      Trans_N_Bonf = trans_stats$n_bonf
    )
    
    summary_data <- rbind(summary_data, row)
  }
}

cat("=== MR Summary Statistics: BMI and HbA1c ===\n\n")
print(summary_data)
cat("\n")

# Save to CSV
output_path <- file.path(results_dir, "MR_summary_BMI_HbA1c_cis_trans.csv")
write.csv(summary_data, output_path, row.names = FALSE)
cat(sprintf("Summary saved to: %s\n", output_path))
