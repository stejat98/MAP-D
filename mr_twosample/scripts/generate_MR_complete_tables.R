#!/usr/bin/env Rscript
# Generate Complete MR Tables - ALL protein associations
# BMI and HbA1c across strata: full, prediabetes, diabetes
# Separate tables for: All pQTLs, Cis-only, Trans-only

# ==============================================================================
# SETUP
# ==============================================================================

results_dir <- "/n/groups/patel/sivateja/regenie_pipeline/results/twosampleMR"
output_dir <- file.path(results_dir, "supplemental_tables")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("=============================================================\n")
cat("Generating Complete MR Tables (All Protein Associations)\n")
cat("=============================================================\n\n")

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

format_pval <- function(p) {
  sapply(p, function(x) {
    if (is.na(x)) return("NA")
    if (x < 0.001) return(sprintf("%.2e", x))
    if (x < 0.01) return(sprintf("%.3f", x))
    return(sprintf("%.3f", x))
  })
}

# ==============================================================================
# READ ALL MR RESULTS
# ==============================================================================

cat("Reading MR results files...\n")

# Overall MR files
mr_files <- list(
  BMI_full = "MR_proteins_BMI_full.csv",
  BMI_prediabetes = "MR_proteins_BMI_prediabetes.csv",
  BMI_diabetes = "MR_proteins_BMI_diabetes.csv",
  HbA1c_full = "MR_proteins_HbA1c_full.csv",
  HbA1c_prediabetes = "MR_proteins_HbA1c_prediabetes.csv",
  HbA1c_diabetes = "MR_proteins_HbA1c_diabetes.csv"
)

mr_data <- list()
for (key in names(mr_files)) {
  filepath <- file.path(results_dir, mr_files[[key]])
  if (file.exists(filepath)) {
    mr_data[[key]] <- read.csv(filepath)
    cat(sprintf("  Read %s: %d rows\n", key, nrow(mr_data[[key]])))
  }
}

# Read cis/trans MR results
cis_trans_dir <- file.path(results_dir, "cis_trans_chunked")
cis_files <- list.files(cis_trans_dir, pattern = "MR_cis_chunk.*\\.csv", full.names = TRUE)
trans_files <- list.files(cis_trans_dir, pattern = "MR_trans_chunk.*\\.csv", full.names = TRUE)

cis_df <- do.call(rbind, lapply(cis_files, read.csv))
trans_df <- do.call(rbind, lapply(trans_files, read.csv))

# Filter for BMI and HbA1c
cis_df <- cis_df[cis_df$phenotype %in% c("BMI", "HbA1c"), ]
trans_df <- trans_df[trans_df$phenotype %in% c("BMI", "HbA1c"), ]

cat(sprintf("\nCis-only results (BMI/HbA1c): %d rows\n", nrow(cis_df)))
cat(sprintf("Trans-only results (BMI/HbA1c): %d rows\n\n", nrow(trans_df)))

# ==============================================================================
# FUNCTION: FORMAT AND SAVE COMPLETE TABLE
# ==============================================================================

create_complete_table <- function(df, protein_col, phenotype, stratum, qtl_type, n_total_proteins) {
  
  # Calculate Bonferroni threshold
  bonf_threshold <- 0.05 / n_total_proteins
  
  # Format stratum display
  stratum_display <- switch(stratum,
    "full" = "Full cohort",
    "prediabetes" = "Prediabetes",
    "diabetes" = "T2D",
    stratum
  )
  
  # Create formatted table
  result <- data.frame(
    Protein = df[[protein_col]],
    N_SNPs = df$nsnp,
    Beta = sprintf("%.4f", df$b),
    SE = sprintf("%.4f", df$se),
    P_value = format_pval(df$pval),
    P_value_numeric = df$pval,
    OR = sprintf("%.3f", df$OR),
    OR_95CI_Lower = sprintf("%.3f", df$OR_lci95),
    OR_95CI_Upper = sprintf("%.3f", df$OR_uci95),
    Nominal_Sig = ifelse(df$pval < 0.05, "Yes", "No"),
    Bonferroni_Sig = ifelse(df$pval < bonf_threshold, "Yes", "No"),
    stringsAsFactors = FALSE
  )
  
  # Sort by p-value
  result <- result[order(result$P_value_numeric), ]
  
  # Remove numeric p-value column (used only for sorting)
  result$P_value_numeric <- NULL
  
  # Rename columns for publication
  colnames(result) <- c(
    "Protein", "N SNPs", "Beta", "SE", "P-value",
    "OR", "OR 95% CI Lower", "OR 95% CI Upper",
    "P<0.05", "Bonferroni Significant"
  )
  
  return(result)
}

# ==============================================================================
# GENERATE TABLES: ALL pQTLs (OVERALL)
# ==============================================================================

cat("Generating tables for ALL pQTLs (Overall)...\n")

for (phenotype in c("BMI", "HbA1c")) {
  for (stratum in c("full", "prediabetes", "diabetes")) {
    key <- paste0(phenotype, "_", stratum)
    
    if (!(key %in% names(mr_data))) next
    
    df <- mr_data[[key]]
    
    # Use IVW method only (primary analysis)
    ivw_df <- df[df$method == "Inverse variance weighted", ]
    n_proteins <- length(unique(ivw_df$exposure))
    
    if (nrow(ivw_df) == 0) next
    
    # Create complete table
    complete_table <- create_complete_table(
      ivw_df, "exposure", phenotype, stratum, "all", n_proteins
    )
    
    # Create filename
    stratum_label <- switch(stratum,
      "full" = "Full",
      "prediabetes" = "Prediabetes", 
      "diabetes" = "T2D"
    )
    
    filename <- sprintf("Table_MR_All_pQTL_%s_%s.csv", phenotype, stratum_label)
    filepath <- file.path(output_dir, filename)
    
    write.csv(complete_table, filepath, row.names = FALSE)
    cat(sprintf("  Saved: %s (%d proteins)\n", filename, nrow(complete_table)))
  }
}

# ==============================================================================
# GENERATE TABLES: CIS-ONLY pQTLs
# ==============================================================================

cat("\nGenerating tables for CIS-ONLY pQTLs...\n")

for (phenotype in c("BMI", "HbA1c")) {
  for (stratum in c("full", "prediabetes", "diabetes")) {
    
    sub_df <- cis_df[cis_df$phenotype == phenotype & cis_df$stratum == stratum, ]
    
    if (nrow(sub_df) == 0) next
    
    n_proteins <- length(unique(sub_df$protein))
    
    # Create complete table
    complete_table <- create_complete_table(
      sub_df, "protein", phenotype, stratum, "cis", n_proteins
    )
    
    # Create filename
    stratum_label <- switch(stratum,
      "full" = "Full",
      "prediabetes" = "Prediabetes",
      "diabetes" = "T2D"
    )
    
    filename <- sprintf("Table_MR_Cis_pQTL_%s_%s.csv", phenotype, stratum_label)
    filepath <- file.path(output_dir, filename)
    
    write.csv(complete_table, filepath, row.names = FALSE)
    cat(sprintf("  Saved: %s (%d proteins)\n", filename, nrow(complete_table)))
  }
}

# ==============================================================================
# GENERATE TABLES: TRANS-ONLY pQTLs
# ==============================================================================

cat("\nGenerating tables for TRANS-ONLY pQTLs...\n")

for (phenotype in c("BMI", "HbA1c")) {
  for (stratum in c("full", "prediabetes", "diabetes")) {
    
    sub_df <- trans_df[trans_df$phenotype == phenotype & trans_df$stratum == stratum, ]
    
    if (nrow(sub_df) == 0) next
    
    # Deduplicate - keep one result per protein (prefer IVW over Wald ratio)
    if ("method" %in% names(sub_df)) {
      sub_dedup <- do.call(rbind, by(sub_df, sub_df$protein, function(x) {
        if (any(x$method == "Inverse variance weighted")) {
          x[x$method == "Inverse variance weighted", ][1, ]
        } else {
          x[1, ]
        }
      }))
    } else {
      sub_dedup <- sub_df
    }
    
    n_proteins <- length(unique(sub_dedup$protein))
    
    # Create complete table
    complete_table <- create_complete_table(
      sub_dedup, "protein", phenotype, stratum, "trans", n_proteins
    )
    
    # Create filename
    stratum_label <- switch(stratum,
      "full" = "Full",
      "prediabetes" = "Prediabetes",
      "diabetes" = "T2D"
    )
    
    filename <- sprintf("Table_MR_Trans_pQTL_%s_%s.csv", phenotype, stratum_label)
    filepath <- file.path(output_dir, filename)
    
    write.csv(complete_table, filepath, row.names = FALSE)
    cat(sprintf("  Saved: %s (%d proteins)\n", filename, nrow(complete_table)))
  }
}

# ==============================================================================
# GENERATE COMBINED MASTER TABLES (All strata in one file per phenotype/QTL type)
# ==============================================================================

cat("\nGenerating combined master tables...\n")

# Function to create combined table across strata
create_master_table <- function(phenotype, qtl_type) {
  
  combined_rows <- list()
  
  for (stratum in c("full", "prediabetes", "diabetes")) {
    
    if (qtl_type == "all") {
      key <- paste0(phenotype, "_", stratum)
      if (!(key %in% names(mr_data))) next
      df <- mr_data[[key]]
      df <- df[df$method == "Inverse variance weighted", ]
      protein_col <- "exposure"
    } else if (qtl_type == "cis") {
      df <- cis_df[cis_df$phenotype == phenotype & cis_df$stratum == stratum, ]
      protein_col <- "protein"
    } else if (qtl_type == "trans") {
      df <- trans_df[trans_df$phenotype == phenotype & trans_df$stratum == stratum, ]
      # Deduplicate
      if ("method" %in% names(df) && nrow(df) > 0) {
        df <- do.call(rbind, by(df, df$protein, function(x) {
          if (any(x$method == "Inverse variance weighted")) {
            x[x$method == "Inverse variance weighted", ][1, ]
          } else {
            x[1, ]
          }
        }))
      }
      protein_col <- "protein"
    }
    
    if (nrow(df) == 0) next
    
    n_proteins <- length(unique(df[[protein_col]]))
    bonf_threshold <- 0.05 / n_proteins
    
    stratum_display <- switch(stratum,
      "full" = "Full cohort",
      "prediabetes" = "Prediabetes",
      "diabetes" = "T2D"
    )
    
    row_data <- data.frame(
      Stratum = stratum_display,
      Protein = df[[protein_col]],
      N_SNPs = df$nsnp,
      Beta = sprintf("%.4f", df$b),
      SE = sprintf("%.4f", df$se),
      P_value = format_pval(df$pval),
      P_value_numeric = df$pval,
      OR = sprintf("%.3f", df$OR),
      OR_95CI_Lower = sprintf("%.3f", df$OR_lci95),
      OR_95CI_Upper = sprintf("%.3f", df$OR_uci95),
      Nominal_Sig = ifelse(df$pval < 0.05, "Yes", "No"),
      Bonferroni_Sig = ifelse(df$pval < bonf_threshold, "Yes", "No"),
      Bonf_Threshold = sprintf("%.2e", bonf_threshold),
      stringsAsFactors = FALSE
    )
    
    combined_rows[[stratum]] <- row_data
  }
  
  if (length(combined_rows) == 0) return(NULL)
  
  combined <- do.call(rbind, combined_rows)
  combined <- combined[order(combined$Stratum, combined$P_value_numeric), ]
  combined$P_value_numeric <- NULL
  rownames(combined) <- NULL
  
  colnames(combined) <- c(
    "Stratum", "Protein", "N SNPs", "Beta", "SE", "P-value",
    "OR", "OR 95% CI Lower", "OR 95% CI Upper",
    "P<0.05", "Bonferroni Sig", "Bonf Threshold"
  )
  
  return(combined)
}

# Generate master tables
for (phenotype in c("BMI", "HbA1c")) {
  for (qtl_type in c("all", "cis", "trans")) {
    
    master_table <- create_master_table(phenotype, qtl_type)
    
    if (is.null(master_table)) next
    
    qtl_label <- switch(qtl_type,
      "all" = "All_pQTL",
      "cis" = "Cis_pQTL",
      "trans" = "Trans_pQTL"
    )
    
    filename <- sprintf("Table_MR_%s_%s_All_Strata.csv", qtl_label, phenotype)
    filepath <- file.path(output_dir, filename)
    
    write.csv(master_table, filepath, row.names = FALSE)
    cat(sprintf("  Saved: %s (%d total rows)\n", filename, nrow(master_table)))
  }
}

# ==============================================================================
# GENERATE XLSX-STYLE SUMMARY INDEX
# ==============================================================================

cat("\nGenerating table index...\n")

index_data <- data.frame(
  Table_Name = c(
    # Individual tables - All pQTLs
    "Table_MR_All_pQTL_BMI_Full.csv",
    "Table_MR_All_pQTL_BMI_Prediabetes.csv",
    "Table_MR_All_pQTL_BMI_T2D.csv",
    "Table_MR_All_pQTL_HbA1c_Full.csv",
    "Table_MR_All_pQTL_HbA1c_Prediabetes.csv",
    "Table_MR_All_pQTL_HbA1c_T2D.csv",
    # Individual tables - Cis pQTLs
    "Table_MR_Cis_pQTL_BMI_Full.csv",
    "Table_MR_Cis_pQTL_BMI_Prediabetes.csv",
    "Table_MR_Cis_pQTL_BMI_T2D.csv",
    "Table_MR_Cis_pQTL_HbA1c_Full.csv",
    "Table_MR_Cis_pQTL_HbA1c_Prediabetes.csv",
    "Table_MR_Cis_pQTL_HbA1c_T2D.csv",
    # Individual tables - Trans pQTLs
    "Table_MR_Trans_pQTL_BMI_Full.csv",
    "Table_MR_Trans_pQTL_BMI_Prediabetes.csv",
    "Table_MR_Trans_pQTL_BMI_T2D.csv",
    "Table_MR_Trans_pQTL_HbA1c_Full.csv",
    "Table_MR_Trans_pQTL_HbA1c_Prediabetes.csv",
    "Table_MR_Trans_pQTL_HbA1c_T2D.csv",
    # Master tables
    "Table_MR_All_pQTL_BMI_All_Strata.csv",
    "Table_MR_All_pQTL_HbA1c_All_Strata.csv",
    "Table_MR_Cis_pQTL_BMI_All_Strata.csv",
    "Table_MR_Cis_pQTL_HbA1c_All_Strata.csv",
    "Table_MR_Trans_pQTL_BMI_All_Strata.csv",
    "Table_MR_Trans_pQTL_HbA1c_All_Strata.csv"
  ),
  Description = c(
    # All pQTLs
    "All protein-BMI MR results (Full cohort) - All pQTLs",
    "All protein-BMI MR results (Prediabetes) - All pQTLs",
    "All protein-BMI MR results (T2D) - All pQTLs",
    "All protein-HbA1c MR results (Full cohort) - All pQTLs",
    "All protein-HbA1c MR results (Prediabetes) - All pQTLs",
    "All protein-HbA1c MR results (T2D) - All pQTLs",
    # Cis pQTLs
    "All protein-BMI MR results (Full cohort) - Cis-pQTLs only",
    "All protein-BMI MR results (Prediabetes) - Cis-pQTLs only",
    "All protein-BMI MR results (T2D) - Cis-pQTLs only",
    "All protein-HbA1c MR results (Full cohort) - Cis-pQTLs only",
    "All protein-HbA1c MR results (Prediabetes) - Cis-pQTLs only",
    "All protein-HbA1c MR results (T2D) - Cis-pQTLs only",
    # Trans pQTLs
    "All protein-BMI MR results (Full cohort) - Trans-pQTLs only",
    "All protein-BMI MR results (Prediabetes) - Trans-pQTLs only",
    "All protein-BMI MR results (T2D) - Trans-pQTLs only",
    "All protein-HbA1c MR results (Full cohort) - Trans-pQTLs only",
    "All protein-HbA1c MR results (Prediabetes) - Trans-pQTLs only",
    "All protein-HbA1c MR results (T2D) - Trans-pQTLs only",
    # Master tables
    "Combined BMI MR results across all strata - All pQTLs",
    "Combined HbA1c MR results across all strata - All pQTLs",
    "Combined BMI MR results across all strata - Cis-pQTLs only",
    "Combined HbA1c MR results across all strata - Cis-pQTLs only",
    "Combined BMI MR results across all strata - Trans-pQTLs only",
    "Combined HbA1c MR results across all strata - Trans-pQTLs only"
  ),
  stringsAsFactors = FALSE
)

write.csv(index_data, file.path(output_dir, "TABLE_INDEX.csv"), row.names = FALSE)
cat("  Saved: TABLE_INDEX.csv\n")

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n=============================================================\n")
cat("COMPLETE MR TABLES GENERATED SUCCESSFULLY\n")
cat("=============================================================\n")
cat(sprintf("Output directory: %s\n\n", output_dir))

cat("Individual tables (per phenotype/stratum/QTL type):\n")
cat("  - 6 tables for All pQTLs (BMI/HbA1c x 3 strata)\n")
cat("  - 6 tables for Cis-only pQTLs (BMI/HbA1c x 3 strata)\n")
cat("  - 6 tables for Trans-only pQTLs (BMI/HbA1c x 3 strata)\n\n")

cat("Master tables (combined across strata):\n")
cat("  - 2 tables for All pQTLs (BMI, HbA1c)\n")
cat("  - 2 tables for Cis-only pQTLs (BMI, HbA1c)\n")
cat("  - 2 tables for Trans-only pQTLs (BMI, HbA1c)\n\n")

cat("Table columns:\n")
cat("  Protein, N SNPs, Beta, SE, P-value, OR, 95% CI, P<0.05, Bonferroni Sig\n")
cat("=============================================================\n")
