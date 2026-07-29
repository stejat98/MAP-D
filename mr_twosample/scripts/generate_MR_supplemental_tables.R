#!/usr/bin/env Rscript
# Generate Nature Metabolism Style Supplemental Tables for MR Analysis
# BMI and HbA1c across strata: full, prediabetes, diabetes
# Includes: proportions, top hits, formatted tables

# ==============================================================================
# SETUP
# ==============================================================================

results_dir <- "/n/groups/patel/sivateja/regenie_pipeline/results/twosampleMR"
output_dir <- file.path(results_dir, "supplemental_tables")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("=============================================================\n")
cat("Generating Nature Metabolism Style Supplemental Tables for MR\n")
cat("=============================================================\n\n")

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

# Format p-values for publication
format_pval <- function(p) {
  sapply(p, function(x) {
    if (is.na(x)) return("NA")
    if (x < 0.001) return(sprintf("%.2e", x))
    if (x < 0.01) return(sprintf("%.3f", x))
    return(sprintf("%.2f", x))
  })
}

# Format effect estimates (beta or OR) with CI
format_estimate <- function(est, lo, hi, digits = 2) {
  sprintf("%.*f (%.*f, %.*f)", digits, est, digits, lo, digits, hi)
}

# Calculate proportion with percentage
calc_proportion <- function(n, total) {
  if (total == 0) return("0/0 (0.0%)")
  pct <- (n / total) * 100
  sprintf("%d/%d (%.1f%%)", n, total, pct)
}

# ==============================================================================
# READ AND COMBINE MR RESULTS
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

# Read all MR data
mr_data <- list()
for (key in names(mr_files)) {
  filepath <- file.path(results_dir, mr_files[[key]])
  if (file.exists(filepath)) {
    mr_data[[key]] <- read.csv(filepath)
    cat(sprintf("  Read %s: %d rows\n", key, nrow(mr_data[[key]])))
  } else {
    cat(sprintf("  WARNING: File not found: %s\n", filepath))
  }
}

# Read cis/trans MR results
cis_trans_dir <- file.path(results_dir, "cis_trans_chunked")
cis_files <- list.files(cis_trans_dir, pattern = "MR_cis_chunk.*\\.csv", full.names = TRUE)
trans_files <- list.files(cis_trans_dir, pattern = "MR_trans_chunk.*\\.csv", full.names = TRUE)

cis_df <- do.call(rbind, lapply(cis_files, read.csv))
trans_df <- do.call(rbind, lapply(trans_files, read.csv))

# Filter for BMI and HbA1c only
cis_df <- cis_df[cis_df$phenotype %in% c("BMI", "HbA1c"), ]
trans_df <- trans_df[trans_df$phenotype %in% c("BMI", "HbA1c"), ]

cat(sprintf("\nCis-only results (BMI/HbA1c): %d rows\n", nrow(cis_df)))
cat(sprintf("Trans-only results (BMI/HbA1c): %d rows\n\n", nrow(trans_df)))

# ==============================================================================
# TABLE S1: SUMMARY TABLE WITH PROPORTIONS
# ==============================================================================

cat("Creating Table S1: Summary Statistics with Proportions...\n")

summary_rows <- list()

for (phenotype in c("BMI", "HbA1c")) {
  for (stratum in c("full", "prediabetes", "diabetes")) {
    key <- paste0(phenotype, "_", stratum)
    
    # Overall MR (IVW method)
    if (key %in% names(mr_data)) {
      df <- mr_data[[key]]
      ivw_df <- df[df$method == "Inverse variance weighted", ]
      n_total <- length(unique(ivw_df$exposure))
      n_sig_p05 <- sum(ivw_df$pval < 0.05, na.rm = TRUE)
      n_bonf <- sum(ivw_df$pval < 0.05/n_total, na.rm = TRUE)
    } else {
      n_total <- n_sig_p05 <- n_bonf <- 0
    }
    
    # Cis MR (Wald ratio for single SNP)
    cis_sub <- cis_df[cis_df$phenotype == phenotype & cis_df$stratum == stratum, ]
    n_cis_total <- length(unique(cis_sub$protein))
    n_cis_sig <- sum(cis_sub$pval < 0.05, na.rm = TRUE)
    n_cis_bonf <- ifelse(n_cis_total > 0, sum(cis_sub$pval < 0.05/n_cis_total, na.rm = TRUE), 0)
    
    # Trans MR (IVW or Wald ratio)
    trans_sub <- trans_df[trans_df$phenotype == phenotype & trans_df$stratum == stratum, ]
    # Get one result per protein (prefer IVW over Wald ratio)
    if (nrow(trans_sub) > 0 && "method" %in% names(trans_sub)) {
      trans_dedup <- do.call(rbind, by(trans_sub, trans_sub$protein, function(x) {
        if (any(x$method == "Inverse variance weighted")) {
          x[x$method == "Inverse variance weighted", ][1, ]
        } else {
          x[1, ]
        }
      }))
    } else {
      trans_dedup <- trans_sub
    }
    n_trans_total <- length(unique(trans_dedup$protein))
    n_trans_sig <- sum(trans_dedup$pval < 0.05, na.rm = TRUE)
    n_trans_bonf <- ifelse(n_trans_total > 0, sum(trans_dedup$pval < 0.05/n_trans_total, na.rm = TRUE), 0)
    
    # Stratum display name
    stratum_display <- switch(stratum,
      "full" = "Full cohort",
      "prediabetes" = "Prediabetes",
      "diabetes" = "T2D"
    )
    
    row <- data.frame(
      Phenotype = phenotype,
      Stratum = stratum_display,
      # Overall
      Overall_N = n_total,
      Overall_Sig_p05 = calc_proportion(n_sig_p05, n_total),
      Overall_Bonf = calc_proportion(n_bonf, n_total),
      # Cis
      Cis_N = n_cis_total,
      Cis_Sig_p05 = calc_proportion(n_cis_sig, n_cis_total),
      Cis_Bonf = calc_proportion(n_cis_bonf, n_cis_total),
      # Trans
      Trans_N = n_trans_total,
      Trans_Sig_p05 = calc_proportion(n_trans_sig, n_trans_total),
      Trans_Bonf = calc_proportion(n_trans_bonf, n_trans_total),
      stringsAsFactors = FALSE
    )
    
    summary_rows[[length(summary_rows) + 1]] <- row
  }
}

summary_table <- do.call(rbind, summary_rows)

# Rename columns for publication
colnames(summary_table) <- c(
  "Phenotype", "Stratum",
  "Overall N", "Overall P<0.05", "Overall Bonferroni",
  "Cis-pQTL N", "Cis P<0.05", "Cis Bonferroni",
  "Trans-pQTL N", "Trans P<0.05", "Trans Bonferroni"
)

# Save Table S1
write.csv(summary_table, file.path(output_dir, "Table_S1_MR_Summary_Statistics.csv"), row.names = FALSE)
cat("  Saved: Table_S1_MR_Summary_Statistics.csv\n")

# ==============================================================================
# TABLE S2: TOP HITS - SIGNIFICANT MR RESULTS
# ==============================================================================

cat("\nCreating Table S2: Top MR Hits...\n")

get_top_hits <- function(df, key, n_top = 20) {
  ivw_df <- df[df$method == "Inverse variance weighted", ]
  ivw_df <- ivw_df[order(ivw_df$pval), ]
  
  # Get total N for Bonferroni
  n_total <- length(unique(ivw_df$exposure))
  bonf_threshold <- 0.05 / n_total
  
  # Select top N or all significant
  sig_df <- ivw_df[ivw_df$pval < 0.05, ]
  top_df <- head(ivw_df, n_top)
  
  if (nrow(top_df) == 0) return(NULL)
  
  # Parse phenotype and stratum from key
  parts <- strsplit(key, "_")[[1]]
  phenotype <- parts[1]
  stratum <- parts[2]
  
  stratum_display <- switch(stratum,
    "full" = "Full cohort",
    "prediabetes" = "Prediabetes",
    "diabetes" = "T2D"
  )
  
  result <- data.frame(
    Phenotype = phenotype,
    Stratum = stratum_display,
    Protein = top_df$exposure,
    N_SNPs = top_df$nsnp,
    Beta = sprintf("%.3f", top_df$b),
    SE = sprintf("%.3f", top_df$se),
    OR = sprintf("%.2f", top_df$OR),
    OR_95CI = sprintf("(%.2f, %.2f)", top_df$OR_lci95, top_df$OR_uci95),
    P_value = format_pval(top_df$pval),
    Bonferroni_Sig = ifelse(top_df$pval < bonf_threshold, "Yes", "No"),
    stringsAsFactors = FALSE
  )
  
  return(result)
}

# Collect top hits for all phenotypes/strata
top_hits_list <- list()
for (key in names(mr_data)) {
  hits <- get_top_hits(mr_data[[key]], key)
  if (!is.null(hits)) {
    top_hits_list[[key]] <- hits
  }
}

top_hits_table <- do.call(rbind, top_hits_list)
rownames(top_hits_table) <- NULL

# Rename columns for publication
colnames(top_hits_table) <- c(
  "Phenotype", "Stratum", "Protein", "N SNPs",
  "Beta", "SE", "OR", "95% CI", "P-value", "Bonferroni Significant"
)

write.csv(top_hits_table, file.path(output_dir, "Table_S2_MR_Top_Hits.csv"), row.names = FALSE)
cat("  Saved: Table_S2_MR_Top_Hits.csv\n")

# ==============================================================================
# TABLE S3: CIS-PQTL TOP HITS
# ==============================================================================

cat("\nCreating Table S3: Cis-pQTL Top Hits...\n")

cis_top_list <- list()

for (phenotype in c("BMI", "HbA1c")) {
  for (stratum in c("full", "prediabetes", "diabetes")) {
    sub_df <- cis_df[cis_df$phenotype == phenotype & cis_df$stratum == stratum, ]
    if (nrow(sub_df) == 0) next
    
    sub_df <- sub_df[order(sub_df$pval), ]
    n_total <- length(unique(sub_df$protein))
    bonf_threshold <- 0.05 / n_total
    
    top_df <- head(sub_df, 20)
    
    stratum_display <- switch(stratum,
      "full" = "Full cohort",
      "prediabetes" = "Prediabetes",
      "diabetes" = "T2D"
    )
    
    result <- data.frame(
      Phenotype = phenotype,
      Stratum = stratum_display,
      Protein = top_df$protein,
      N_SNPs = top_df$nsnp,
      Beta = sprintf("%.3f", top_df$b),
      SE = sprintf("%.3f", top_df$se),
      OR = sprintf("%.2f", top_df$OR),
      OR_95CI = sprintf("(%.2f, %.2f)", top_df$OR_lci95, top_df$OR_uci95),
      P_value = format_pval(top_df$pval),
      Bonferroni_Sig = ifelse(top_df$pval < bonf_threshold, "Yes", "No"),
      stringsAsFactors = FALSE
    )
    
    cis_top_list[[paste(phenotype, stratum, sep = "_")]] <- result
  }
}

cis_top_table <- do.call(rbind, cis_top_list)
rownames(cis_top_table) <- NULL

colnames(cis_top_table) <- c(
  "Phenotype", "Stratum", "Protein", "N SNPs",
  "Beta", "SE", "OR", "95% CI", "P-value", "Bonferroni Significant"
)

write.csv(cis_top_table, file.path(output_dir, "Table_S3_Cis_pQTL_Top_Hits.csv"), row.names = FALSE)
cat("  Saved: Table_S3_Cis_pQTL_Top_Hits.csv\n")

# ==============================================================================
# TABLE S4: TRANS-PQTL TOP HITS
# ==============================================================================

cat("\nCreating Table S4: Trans-pQTL Top Hits...\n")

trans_top_list <- list()

for (phenotype in c("BMI", "HbA1c")) {
  for (stratum in c("full", "prediabetes", "diabetes")) {
    sub_df <- trans_df[trans_df$phenotype == phenotype & trans_df$stratum == stratum, ]
    if (nrow(sub_df) == 0) next
    
    # Deduplicate - prefer IVW
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
    
    sub_dedup <- sub_dedup[order(sub_dedup$pval), ]
    n_total <- length(unique(sub_dedup$protein))
    bonf_threshold <- 0.05 / n_total
    
    top_df <- head(sub_dedup, 20)
    
    stratum_display <- switch(stratum,
      "full" = "Full cohort",
      "prediabetes" = "Prediabetes",
      "diabetes" = "T2D"
    )
    
    method_col <- if ("method" %in% names(top_df)) top_df$method else rep("IVW", nrow(top_df))
    
    result <- data.frame(
      Phenotype = phenotype,
      Stratum = stratum_display,
      Protein = top_df$protein,
      Method = method_col,
      N_SNPs = top_df$nsnp,
      Beta = sprintf("%.3f", top_df$b),
      SE = sprintf("%.3f", top_df$se),
      OR = sprintf("%.2f", top_df$OR),
      OR_95CI = sprintf("(%.2f, %.2f)", top_df$OR_lci95, top_df$OR_uci95),
      P_value = format_pval(top_df$pval),
      Bonferroni_Sig = ifelse(top_df$pval < bonf_threshold, "Yes", "No"),
      stringsAsFactors = FALSE
    )
    
    trans_top_list[[paste(phenotype, stratum, sep = "_")]] <- result
  }
}

trans_top_table <- do.call(rbind, trans_top_list)
rownames(trans_top_table) <- NULL

colnames(trans_top_table) <- c(
  "Phenotype", "Stratum", "Protein", "MR Method", "N SNPs",
  "Beta", "SE", "OR", "95% CI", "P-value", "Bonferroni Significant"
)

write.csv(trans_top_table, file.path(output_dir, "Table_S4_Trans_pQTL_Top_Hits.csv"), row.names = FALSE)
cat("  Saved: Table_S4_Trans_pQTL_Top_Hits.csv\n")

# ==============================================================================
# TABLE S5: DETAILED BONFERRONI-SIGNIFICANT RESULTS
# ==============================================================================

cat("\nCreating Table S5: Bonferroni-Significant Results (Detailed)...\n")

bonf_results_list <- list()

for (key in names(mr_data)) {
  df <- mr_data[[key]]
  ivw_df <- df[df$method == "Inverse variance weighted", ]
  n_total <- length(unique(ivw_df$exposure))
  bonf_threshold <- 0.05 / n_total
  
  bonf_df <- ivw_df[ivw_df$pval < bonf_threshold, ]
  if (nrow(bonf_df) == 0) next
  
  parts <- strsplit(key, "_")[[1]]
  phenotype <- parts[1]
  stratum <- parts[2]
  
  stratum_display <- switch(stratum,
    "full" = "Full cohort",
    "prediabetes" = "Prediabetes",
    "diabetes" = "T2D"
  )
  
  result <- data.frame(
    Phenotype = phenotype,
    Stratum = stratum_display,
    Protein = bonf_df$exposure,
    N_SNPs = bonf_df$nsnp,
    Beta = sprintf("%.4f", bonf_df$b),
    SE = sprintf("%.4f", bonf_df$se),
    OR = sprintf("%.3f", bonf_df$OR),
    OR_95CI_Lower = sprintf("%.3f", bonf_df$OR_lci95),
    OR_95CI_Upper = sprintf("%.3f", bonf_df$OR_uci95),
    P_value = format_pval(bonf_df$pval),
    Bonf_Threshold = sprintf("%.2e", bonf_threshold),
    stringsAsFactors = FALSE
  )
  
  bonf_results_list[[key]] <- result
}

if (length(bonf_results_list) > 0) {
  bonf_table <- do.call(rbind, bonf_results_list)
  rownames(bonf_table) <- NULL
  
  colnames(bonf_table) <- c(
    "Phenotype", "Stratum", "Protein", "N SNPs",
    "Beta", "SE", "OR", "OR 95% CI Lower", "OR 95% CI Upper",
    "P-value", "Bonferroni Threshold"
  )
  
  write.csv(bonf_table, file.path(output_dir, "Table_S5_Bonferroni_Significant_Results.csv"), row.names = FALSE)
  cat("  Saved: Table_S5_Bonferroni_Significant_Results.csv\n")
} else {
  cat("  No Bonferroni-significant results found.\n")
}

# ==============================================================================
# TABLE S6: CROSS-STRATUM COMPARISON
# ==============================================================================

cat("\nCreating Table S6: Cross-Stratum Comparison...\n")

# Get unique proteins across all analyses
all_proteins <- unique(unlist(lapply(mr_data, function(df) {
  ivw <- df[df$method == "Inverse variance weighted", ]
  unique(ivw$exposure)
})))

cross_stratum_list <- list()

for (phenotype in c("BMI", "HbA1c")) {
  for (protein in all_proteins) {
    row_data <- list(Phenotype = phenotype, Protein = protein)
    has_data <- FALSE
    
    for (stratum in c("full", "prediabetes", "diabetes")) {
      key <- paste0(phenotype, "_", stratum)
      if (key %in% names(mr_data)) {
        df <- mr_data[[key]]
        ivw_df <- df[df$method == "Inverse variance weighted" & df$exposure == protein, ]
        
        if (nrow(ivw_df) > 0) {
          has_data <- TRUE
          row_data[[paste0(stratum, "_beta")]] <- sprintf("%.3f", ivw_df$b[1])
          row_data[[paste0(stratum, "_pval")]] <- format_pval(ivw_df$pval[1])
          row_data[[paste0(stratum, "_sig")]] <- ifelse(ivw_df$pval[1] < 0.05, "*", "")
        } else {
          row_data[[paste0(stratum, "_beta")]] <- "NA"
          row_data[[paste0(stratum, "_pval")]] <- "NA"
          row_data[[paste0(stratum, "_sig")]] <- ""
        }
      }
    }
    
    # Only include if significant in at least one stratum
    if (has_data) {
      is_sig_any <- any(sapply(c("full_sig", "prediabetes_sig", "diabetes_sig"), function(x) {
        !is.null(row_data[[x]]) && row_data[[x]] == "*"
      }))
      
      if (is_sig_any) {
        cross_stratum_list[[length(cross_stratum_list) + 1]] <- as.data.frame(row_data, stringsAsFactors = FALSE)
      }
    }
  }
}

if (length(cross_stratum_list) > 0) {
  cross_stratum_table <- do.call(rbind, cross_stratum_list)
  
  colnames(cross_stratum_table) <- c(
    "Phenotype", "Protein",
    "Full Beta", "Full P", "Full Sig",
    "Prediabetes Beta", "Prediabetes P", "Prediabetes Sig",
    "T2D Beta", "T2D P", "T2D Sig"
  )
  
  write.csv(cross_stratum_table, file.path(output_dir, "Table_S6_Cross_Stratum_Comparison.csv"), row.names = FALSE)
  cat("  Saved: Table_S6_Cross_Stratum_Comparison.csv\n")
} else {
  cat("  No significant cross-stratum results found.\n")
}

# ==============================================================================
# FORMATTED TEXT TABLE (Reader-Friendly)
# ==============================================================================

cat("\nCreating formatted text summary...\n")

sink(file.path(output_dir, "MR_Supplemental_Tables_Summary.txt"))

cat("================================================================================\n")
cat("SUPPLEMENTARY TABLES: Two-Sample Mendelian Randomization Analysis\n")
cat("BMI and HbA1c Across Glycemic Strata\n")
cat("================================================================================\n\n")

cat("TABLE S1. Summary of MR Analyses\n")
cat("--------------------------------------------------------------------------------\n")
cat("Two-sample MR was performed to assess the causal effect of circulating proteins\n")
cat("on BMI and HbA1c across three glycemic strata: full cohort, prediabetes, and T2D.\n")
cat("Genetic instruments (pQTLs) were used to proxy protein levels.\n\n")

cat("                                    OVERALL                CIS-pQTL              TRANS-pQTL\n")
cat("Phenotype    Stratum         N    P<0.05    Bonf      N    P<0.05    Bonf      N    P<0.05    Bonf\n")
cat("--------------------------------------------------------------------------------\n")

for (i in 1:nrow(summary_table)) {
  row <- summary_table[i, ]
  cat(sprintf("%-10s   %-14s %3s   %-12s %-10s %3s   %-12s %-10s %3s   %-12s %s\n",
              row[[1]], row[[2]],
              row[[3]], row[[4]], row[[5]],
              row[[6]], row[[7]], row[[8]],
              row[[9]], row[[10]], row[[11]]))
}

cat("\n--------------------------------------------------------------------------------\n")
cat("Notes:\n")
cat("  - Overall: Analysis using all available pQTLs per protein\n")
cat("  - Cis-pQTL: Analysis restricted to cis-acting variants (within 1Mb of gene)\n")
cat("  - Trans-pQTL: Analysis restricted to trans-acting variants (>1Mb from gene)\n")
cat("  - P<0.05: Number (proportion) significant at nominal threshold\n")
cat("  - Bonf: Number (proportion) significant after Bonferroni correction\n")
cat("  - IVW method used for multi-SNP analyses; Wald ratio for single-SNP\n\n")

cat("================================================================================\n")
cat("TABLE S2-S5. Detailed results available in CSV format:\n")
cat("  - Table_S2: Top 20 MR hits per phenotype/stratum (overall pQTLs)\n")
cat("  - Table_S3: Top 20 cis-pQTL MR hits per phenotype/stratum\n")
cat("  - Table_S4: Top 20 trans-pQTL MR hits per phenotype/stratum\n")
cat("  - Table_S5: All Bonferroni-significant results with full statistics\n")
cat("  - Table_S6: Cross-stratum comparison of nominally significant proteins\n")
cat("================================================================================\n")

sink()

cat("  Saved: MR_Supplemental_Tables_Summary.txt\n")

# ==============================================================================
# SUMMARY OUTPUT
# ==============================================================================

cat("\n=============================================================\n")
cat("SUPPLEMENTAL TABLES GENERATED SUCCESSFULLY\n")
cat("=============================================================\n")
cat(sprintf("Output directory: %s\n\n", output_dir))
cat("Files created:\n")
cat("  1. Table_S1_MR_Summary_Statistics.csv       - Summary with proportions\n")
cat("  2. Table_S2_MR_Top_Hits.csv                 - Top overall MR hits\n")
cat("  3. Table_S3_Cis_pQTL_Top_Hits.csv           - Top cis-pQTL hits\n")
cat("  4. Table_S4_Trans_pQTL_Top_Hits.csv         - Top trans-pQTL hits\n")
cat("  5. Table_S5_Bonferroni_Significant_Results.csv - Bonferroni hits\n")
cat("  6. Table_S6_Cross_Stratum_Comparison.csv    - Cross-stratum comparison\n")
cat("  7. MR_Supplemental_Tables_Summary.txt       - Reader-friendly summary\n")
cat("=============================================================\n")
