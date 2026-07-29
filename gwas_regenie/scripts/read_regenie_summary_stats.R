#!/usr/bin/env Rscript
# Script to read REGENIE summary statistics files
# Compatible with UK Biobank GWAS format standards

library(data.table)
library(dplyr)

# Function to read REGENIE summary statistics
read_regenie_summary_stats <- function(file_path, 
                                       add_pvalue = TRUE,
                                       convert_to_ukbb_format = FALSE) {
  # Read the gzipped file
  cat("Reading file:", file_path, "\n")
  df <- fread(file_path, sep = " ", header = TRUE, 
              stringsAsFactors = FALSE, data.table = FALSE)
  
  cat("Loaded", nrow(df), "variants\n")
  
  # Convert LOG10P to P-value if requested
  if (add_pvalue && "LOG10P" %in% colnames(df)) {
    df$P <- 10^(-df$LOG10P)
    cat("Added P-value column (from LOG10P)\n")
  }
  
  # Convert to UK Biobank format if requested
  if (convert_to_ukbb_format) {
    df <- convert_to_ukbb_format(df)
  }
  
  return(df)
}

# Function to convert REGENIE format to UK Biobank format
convert_to_ukbb_format <- function(df) {
  # UK Biobank standard column names:
  # CHR, BP, SNP, A1, A2, FRQ, N, BETA, SE, P
  
  df_ukbb <- df %>%
    rename(
      CHR = CHROM,
      BP = GENPOS,
      SNP = ID,
      A1 = ALLELE1,  # Effect allele
      A2 = ALLELE0,  # Reference allele
      FRQ = A1FREQ,
      P = P  # Should already exist if add_pvalue=TRUE
    ) %>%
    select(CHR, BP, SNP, A1, A2, FRQ, N, BETA, SE, P, CHISQ, LOG10P)
  
  cat("Converted to UK Biobank format\n")
  return(df_ukbb)
}

# Function to filter variants (common QC steps)
filter_summary_stats <- function(df, 
                                 min_maf = 0.01,
                                 max_maf = 0.99,
                                 min_n = 100,
                                 remove_missing = TRUE) {
  n_before <- nrow(df)
  
  if (remove_missing) {
    df <- df %>%
      filter(!is.na(BETA), !is.na(SE), !is.na(P),
             BETA != "NA", SE != "NA", P != "NA")
  }
  
  if ("A1FREQ" %in% colnames(df)) {
    df <- df %>%
      filter(A1FREQ >= min_maf, A1FREQ <= max_maf)
  }
  
  if ("N" %in% colnames(df)) {
    df <- df %>%
      filter(N >= min_n)
  }
  
  n_after <- nrow(df)
  cat("Filtered:", n_before, "->", n_after, "variants (removed", 
      n_before - n_after, ")\n")
  
  return(df)
}

# Example usage
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 1) {
    cat("Usage: Rscript read_regenie_summary_stats.R <file_path> [--ukbb] [--filter]\n")
    cat("\nOptions:\n")
    cat("  --ukbb    : Convert to UK Biobank format\n")
    cat("  --filter  : Apply standard QC filters\n")
    cat("\nExample:\n")
    cat("  Rscript read_regenie_summary_stats.R BMI_diabetes_all_chr.regenie.gz --ukbb --filter\n")
    quit(status = 1)
  }
  
  file_path <- args[1]
  convert_ukbb <- "--ukbb" %in% args
  apply_filter <- "--filter" %in% args
  
  # Read the file
  df <- read_regenie_summary_stats(file_path, 
                                   add_pvalue = TRUE,
                                   convert_to_ukbb_format = convert_ukbb)
  
  # Apply filters if requested
  if (apply_filter) {
    df <- filter_summary_stats(df)
  }
  
  # Print summary
  cat("\n=== Summary Statistics ===\n")
  cat("Total variants:", nrow(df), "\n")
  cat("Chromosomes:", paste(sort(unique(df$CHROM)), collapse = ", "), "\n")
  cat("Sample size range:", min(df$N, na.rm = TRUE), "-", 
      max(df$N, na.rm = TRUE), "\n")
  cat("P-value range:", min(df$P, na.rm = TRUE), "-", 
      max(df$P, na.rm = TRUE), "\n")
  
  # Show first few rows
  cat("\n=== First 5 variants ===\n")
  print(head(df, 5))
}
