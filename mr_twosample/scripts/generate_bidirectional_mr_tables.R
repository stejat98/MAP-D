#!/usr/bin/env Rscript
# =============================================================================
# Generate Nature Metabolism-Ready Bidirectional MR Tables
# For: BMI, HbA1c, TRIG_HDL_RATIO
#
# Output per phenotype:
#   1. Full bidirectional table (all proteins, unified columns)
#   2. Summary statistics table
#
# Output combined:
#   3. Cross-phenotype summary
# =============================================================================

setwd("/path/to/project/regenie_pipeline")
output_dir <- "results/twosampleMR/supplemental_tables"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("================================================================================\n")
cat("  Generating Nature Metabolism Bidirectional MR Tables\n")
cat("================================================================================\n\n")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================
format_pval <- function(p) {
  sapply(p, function(x) {
    if (is.na(x)) return(NA_character_)
    if (x == 0) return("< 1e-300")
    if (x < 0.001) return(sprintf("%.2e", x))
    return(sprintf("%.4f", x))
  })
}

format_beta <- function(b, digits = 4) {
  sapply(b, function(x) {
    if (is.na(x)) return(NA_character_)
    sprintf("%.*f", digits, x)
  })
}

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading bidirectional MR results...\n\n")

# --- BMI ---
bmi <- read.csv("results/twosampleMR/bidirectional_bmi_proteins/bidirectional_MR_BMI_all_proteins.csv")
cat(sprintf("  BMI: %d proteins\n", nrow(bmi)))
# Standardize column names
bmi_std <- data.frame(
  Protein            = bmi$Protein,
  Fwd_N_SNPs         = bmi$Fwd_N_SNPs,
  Fwd_Beta           = as.numeric(bmi$Fwd_Beta),
  Fwd_SE             = as.numeric(bmi$Fwd_SE),
  Fwd_P              = as.numeric(bmi$Fwd_P),
  Fwd_Bonf           = bmi$Fwd_Bonf,
  Rev_N_SNPs         = as.numeric(bmi$Rev_N_SNPs),
  Rev_Beta           = as.numeric(bmi$Rev_Beta),
  Rev_SE             = as.numeric(bmi$Rev_SE),
  Rev_P              = as.numeric(bmi$Rev_P),
  Rev_Q              = as.numeric(bmi$Q),
  Rev_Q_P            = as.numeric(bmi$Q_P),
  Rev_Egger_Beta     = as.numeric(bmi$Egger_Beta),
  Rev_Egger_Intercept   = as.numeric(bmi$Egger_Intercept),
  Rev_Egger_Intercept_P = as.numeric(bmi$Egger_Intercept_P),
  Rev_P_BH           = as.numeric(bmi$Rev_P_BH),
  Rev_Bonf           = bmi$Rev_Bonf,
  Fwd_Dir            = bmi$Fwd_Dir,
  Rev_Dir            = bmi$Rev_Dir,
  Concordant         = bmi$Concordant,
  stringsAsFactors   = FALSE
)

# --- HbA1c ---
hba1c <- read.csv("results/twosampleMR/bidirectional_hba1c_proteins/bidirectional_MR_HbA1c_all_proteins.csv")
cat(sprintf("  HbA1c: %d proteins\n", nrow(hba1c)))
hba1c_std <- data.frame(
  Protein            = hba1c$Protein,
  Fwd_N_SNPs         = hba1c$Fwd_N_SNPs,
  Fwd_Beta           = as.numeric(hba1c$Fwd_Beta),
  Fwd_SE             = as.numeric(hba1c$Fwd_SE),
  Fwd_P              = as.numeric(hba1c$Fwd_P),
  Fwd_Bonf           = hba1c$Fwd_Bonf,
  Rev_N_SNPs         = as.numeric(hba1c$N_SNPs),
  Rev_Beta           = as.numeric(hba1c$IVW_Beta),
  Rev_SE             = as.numeric(hba1c$IVW_SE),
  Rev_P              = as.numeric(hba1c$IVW_P),
  Rev_Q              = as.numeric(hba1c$Q),
  Rev_Q_P            = as.numeric(hba1c$Q_P),
  Rev_Egger_Beta     = as.numeric(hba1c$Egger_Beta),
  Rev_Egger_Intercept   = as.numeric(hba1c$Egger_Intercept),
  Rev_Egger_Intercept_P = as.numeric(hba1c$Egger_Intercept_P),
  Rev_P_BH           = as.numeric(hba1c$IVW_P_BH),
  Rev_Bonf           = hba1c$IVW_Bonf,
  Fwd_Dir            = hba1c$Fwd_Dir,
  Rev_Dir            = hba1c$Rev_Dir,
  Concordant         = hba1c$Concordant,
  stringsAsFactors   = FALSE
)

# --- TRIG_HDL_RATIO ---
trig <- read.csv("results/twosampleMR/bidirectional_trig_hdl_ratio_proteins/bidirectional_MR_TRIG_HDL_RATIO_all_proteins.csv")
cat(sprintf("  TRIG_HDL_RATIO: %d proteins\n", nrow(trig)))
trig_std <- data.frame(
  Protein            = trig$Protein,
  Fwd_N_SNPs         = trig$Fwd_N_SNPs,
  Fwd_Beta           = as.numeric(trig$Fwd_Beta),
  Fwd_SE             = as.numeric(trig$Fwd_SE),
  Fwd_P              = as.numeric(trig$Fwd_P),
  Fwd_Bonf           = trig$Fwd_Bonf,
  Rev_N_SNPs         = as.numeric(trig$Rev_N_SNPs),
  Rev_Beta           = as.numeric(trig$Rev_Beta),
  Rev_SE             = as.numeric(trig$Rev_SE),
  Rev_P              = as.numeric(trig$Rev_P),
  Rev_Q              = as.numeric(trig$Q),
  Rev_Q_P            = as.numeric(trig$Q_P),
  Rev_Egger_Beta     = as.numeric(trig$Egger_Beta),
  Rev_Egger_Intercept   = as.numeric(trig$Egger_Intercept),
  Rev_Egger_Intercept_P = as.numeric(trig$Egger_Intercept_P),
  Rev_P_BH           = as.numeric(trig$Rev_P_BH),
  Rev_Bonf           = trig$Rev_Bonf,
  Fwd_Dir            = trig$Fwd_Dir,
  Rev_Dir            = trig$Rev_Dir,
  Concordant         = trig$Concordant,
  stringsAsFactors   = FALSE
)

# =============================================================================
# GENERATE PUBLICATION TABLE FOR EACH PHENOTYPE
# =============================================================================

generate_pub_table <- function(df, phenotype_label, phenotype_short) {

  n_proteins <- nrow(df)
  n_fwd_tested <- sum(!is.na(df$Fwd_P))
  n_rev_tested <- sum(!is.na(df$Rev_P))
  bonf_fwd <- 0.05 / n_fwd_tested
  bonf_rev <- 0.05 / n_rev_tested

  cat(sprintf("\n--- %s ---\n", phenotype_label))
  cat(sprintf("  Proteins: %d (Fwd tested: %d, Rev tested: %d)\n",
              n_proteins, n_fwd_tested, n_rev_tested))
  cat(sprintf("  Bonf thresholds: Fwd=%.2e, Rev=%.2e\n", bonf_fwd, bonf_rev))

  # Publication-ready table
  pub <- data.frame(
    Protein = df$Protein,

    # Forward: Protein -> Phenotype (cis-pQTL MR)
    `Fwd: N Instruments`   = df$Fwd_N_SNPs,
    `Fwd: Beta`            = format_beta(df$Fwd_Beta),
    `Fwd: SE`              = format_beta(df$Fwd_SE),
    `Fwd: P-value`         = format_pval(df$Fwd_P),
    `Fwd: Direction`       = df$Fwd_Dir,
    `Fwd: Bonferroni Sig`  = df$Fwd_Bonf,

    # Reverse: Phenotype -> Protein (IVW MR)
    `Rev: N Instruments`   = df$Rev_N_SNPs,
    `Rev: IVW Beta`        = format_beta(df$Rev_Beta),
    `Rev: IVW SE`          = format_beta(df$Rev_SE),
    `Rev: IVW P-value`     = format_pval(df$Rev_P),
    `Rev: Direction`       = df$Rev_Dir,
    `Rev: Bonferroni Sig`  = df$Rev_Bonf,
    `Rev: FDR P-value`     = format_pval(df$Rev_P_BH),

    # Heterogeneity
    `Rev: Cochran Q`       = format_beta(df$Rev_Q, 2),
    `Rev: Cochran Q P`     = format_pval(df$Rev_Q_P),

    # MR-Egger sensitivity
    `Rev: Egger Beta`            = format_beta(df$Rev_Egger_Beta),
    `Rev: Egger Intercept`       = format_beta(df$Rev_Egger_Intercept),
    `Rev: Egger Intercept P`     = format_pval(df$Rev_Egger_Intercept_P),

    # Concordance
    `Directional Concordance` = df$Concordant,

    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  # Sort: both-sig first (by fwd P), then fwd-only sig, then rev-only sig, then rest
  fwd_sig <- !is.na(df$Fwd_Bonf) & df$Fwd_Bonf == "Yes"
  rev_sig <- !is.na(df$Rev_Bonf) & df$Rev_Bonf == "Yes"
  both_sig <- fwd_sig & rev_sig
  fwd_only <- fwd_sig & !rev_sig
  rev_only <- !fwd_sig & rev_sig
  neither <- !fwd_sig & !rev_sig

  pub$sort_group <- ifelse(both_sig, 1, ifelse(fwd_only, 2, ifelse(rev_only, 3, 4)))
  pub$sort_p <- ifelse(!is.na(df$Fwd_P), df$Fwd_P, 1)
  pub <- pub[order(pub$sort_group, pub$sort_p), ]
  pub$sort_group <- NULL
  pub$sort_p <- NULL

  # Summary counts
  n_fwd_bonf <- sum(fwd_sig, na.rm = TRUE)
  n_rev_bonf <- sum(rev_sig, na.rm = TRUE)
  n_both <- sum(both_sig, na.rm = TRUE)
  n_concordant_both <- sum(both_sig & !is.na(df$Concordant) & df$Concordant == TRUE, na.rm = TRUE)
  n_fwd_only <- sum(fwd_only, na.rm = TRUE)
  n_rev_only <- sum(rev_only, na.rm = TRUE)

  cat(sprintf("  Forward Bonf sig: %d\n", n_fwd_bonf))
  cat(sprintf("  Reverse Bonf sig: %d\n", n_rev_bonf))
  cat(sprintf("  Both sig: %d (concordant: %d, discordant: %d)\n",
              n_both, n_concordant_both, n_both - n_concordant_both))
  cat(sprintf("  Forward-only sig: %d\n", n_fwd_only))
  cat(sprintf("  Reverse-only sig: %d\n", n_rev_only))

  # Save full table
  outfile <- file.path(output_dir,
    sprintf("Table_Bidirectional_MR_%s_Full.csv", phenotype_short))
  write.csv(pub, outfile, row.names = FALSE)
  cat(sprintf("  Saved: %s\n", outfile))

  # ---- Significant-only table ----
  sig_rows <- fwd_sig | rev_sig
  if (sum(sig_rows) > 0) {
    pub_sig <- pub[1:sum(sig_rows), ]  # already sorted with sig first
    outfile_sig <- file.path(output_dir,
      sprintf("Table_Bidirectional_MR_%s_Significant.csv", phenotype_short))
    write.csv(pub_sig, outfile_sig, row.names = FALSE)
    cat(sprintf("  Saved: %s (%d rows)\n", outfile_sig, nrow(pub_sig)))
  }

  # Return summary for cross-phenotype table
  list(
    phenotype = phenotype_label,
    n_proteins = n_proteins,
    n_fwd_tested = n_fwd_tested,
    n_rev_tested = n_rev_tested,
    n_fwd_bonf = n_fwd_bonf,
    n_rev_bonf = n_rev_bonf,
    n_both = n_both,
    n_concordant_both = n_concordant_both,
    n_fwd_only = n_fwd_only,
    n_rev_only = n_rev_only,
    bonf_fwd = bonf_fwd,
    bonf_rev = bonf_rev,
    both_sig_df = df[both_sig, ],
    fwd_only_df = df[fwd_only, ]
  )
}

# Run for all three
bmi_res <- generate_pub_table(bmi_std, "BMI", "BMI")
hba1c_res <- generate_pub_table(hba1c_std, "HbA1c", "HbA1c")
trig_res <- generate_pub_table(trig_std, "TG/HDL-C Ratio", "TRIG_HDL_RATIO")

# =============================================================================
# CROSS-PHENOTYPE SUMMARY TABLE
# =============================================================================

cat("\n\n--- Cross-Phenotype Summary ---\n")

summary_df <- data.frame(
  Phenotype = c("BMI", "HbA1c", "TG/HDL-C Ratio"),
  `Forward: Instrument Source` = rep("UKB cis-pQTLs (split-sample)", 3),
  `Reverse: Instrument Source` = c(
    sprintf("UKB BMI GWAS (%d SNPs)", bmi_std$Rev_N_SNPs[!is.na(bmi_std$Rev_N_SNPs)][1]),
    sprintf("UKB HbA1c GWAS (%d SNPs)", hba1c_std$Rev_N_SNPs[!is.na(hba1c_std$Rev_N_SNPs)][1]),
    sprintf("UKB TG/HDL GWAS (%d SNPs)", trig_std$Rev_N_SNPs[!is.na(trig_std$Rev_N_SNPs)][1])
  ),
  `Proteins Tested` = c(bmi_res$n_proteins, hba1c_res$n_proteins, trig_res$n_proteins),
  `Fwd Bonferroni Threshold` = sprintf("%.2e", c(bmi_res$bonf_fwd, hba1c_res$bonf_fwd, trig_res$bonf_fwd)),
  `Fwd Bonferroni Significant` = c(bmi_res$n_fwd_bonf, hba1c_res$n_fwd_bonf, trig_res$n_fwd_bonf),
  `Rev Bonferroni Significant` = c(bmi_res$n_rev_bonf, hba1c_res$n_rev_bonf, trig_res$n_rev_bonf),
  `Both Directions Significant` = c(bmi_res$n_both, hba1c_res$n_both, trig_res$n_both),
  `Both Sig + Concordant` = c(bmi_res$n_concordant_both, hba1c_res$n_concordant_both, trig_res$n_concordant_both),
  `Forward-Only Significant` = c(bmi_res$n_fwd_only, hba1c_res$n_fwd_only, trig_res$n_fwd_only),
  `Reverse-Only Significant` = c(bmi_res$n_rev_only, hba1c_res$n_rev_only, trig_res$n_rev_only),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

outfile_summary <- file.path(output_dir, "Table_Bidirectional_MR_Summary.csv")
write.csv(summary_df, outfile_summary, row.names = FALSE)
cat(sprintf("  Saved: %s\n", outfile_summary))

# =============================================================================
# KEY FINDINGS TABLE: Both-directions significant proteins
# =============================================================================

cat("\n--- Key Findings: Both Directions Bonferroni Significant ---\n")

make_key_findings <- function(res, pheno_label) {
  d <- res$both_sig_df
  if (nrow(d) == 0) return(NULL)
  data.frame(
    Phenotype = pheno_label,
    Protein = d$Protein,
    `Fwd Beta (Protein -> Phenotype)` = format_beta(d$Fwd_Beta),
    `Fwd P` = format_pval(d$Fwd_P),
    `Rev Beta (Phenotype -> Protein)` = format_beta(d$Rev_Beta),
    `Rev P` = format_pval(d$Rev_P),
    `Egger Intercept` = format_beta(d$Rev_Egger_Intercept),
    `Egger Intercept P` = format_pval(d$Rev_Egger_Intercept_P),
    Concordant = d$Concordant,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

key_both <- rbind(
  make_key_findings(bmi_res, "BMI"),
  make_key_findings(hba1c_res, "HbA1c"),
  make_key_findings(trig_res, "TG/HDL-C Ratio")
)

if (!is.null(key_both) && nrow(key_both) > 0) {
  outfile_key_both <- file.path(output_dir, "Table_Bidirectional_MR_Both_Sig.csv")
  write.csv(key_both, outfile_key_both, row.names = FALSE)
  cat(sprintf("  Saved: %s (%d proteins)\n", outfile_key_both, nrow(key_both)))
}

# =============================================================================
# KEY FINDINGS TABLE: Forward-only significant (causal protein -> phenotype)
# =============================================================================

cat("\n--- Key Findings: Forward-Only Significant (Protein -> Phenotype) ---\n")

make_fwd_only <- function(res, pheno_label) {
  d <- res$fwd_only_df
  if (nrow(d) == 0) return(NULL)
  data.frame(
    Phenotype = pheno_label,
    Protein = d$Protein,
    `Fwd N SNPs` = d$Fwd_N_SNPs,
    `Fwd Beta` = format_beta(d$Fwd_Beta),
    `Fwd SE` = format_beta(d$Fwd_SE),
    `Fwd P` = format_pval(d$Fwd_P),
    `Rev Beta` = format_beta(d$Rev_Beta),
    `Rev P` = format_pval(d$Rev_P),
    `Interpretation` = "Protein likely causal for phenotype (no reverse effect)",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

key_fwd <- rbind(
  make_fwd_only(bmi_res, "BMI"),
  make_fwd_only(hba1c_res, "HbA1c"),
  make_fwd_only(trig_res, "TG/HDL-C Ratio")
)

if (!is.null(key_fwd) && nrow(key_fwd) > 0) {
  outfile_key_fwd <- file.path(output_dir, "Table_Bidirectional_MR_Fwd_Only_Sig.csv")
  write.csv(key_fwd, outfile_key_fwd, row.names = FALSE)
  cat(sprintf("  Saved: %s (%d proteins)\n", outfile_key_fwd, nrow(key_fwd)))
}

# =============================================================================
# PRINT FINAL SUMMARY
# =============================================================================

cat("\n================================================================================\n")
cat("  OUTPUT FILES\n")
cat("================================================================================\n\n")

cat("Full bidirectional tables (all proteins):\n")
cat(sprintf("  1. %s/Table_Bidirectional_MR_BMI_Full.csv\n", output_dir))
cat(sprintf("  2. %s/Table_Bidirectional_MR_HbA1c_Full.csv\n", output_dir))
cat(sprintf("  3. %s/Table_Bidirectional_MR_TRIG_HDL_RATIO_Full.csv\n", output_dir))

cat("\nSignificant-only tables:\n")
cat(sprintf("  4. %s/Table_Bidirectional_MR_BMI_Significant.csv\n", output_dir))
cat(sprintf("  5. %s/Table_Bidirectional_MR_HbA1c_Significant.csv\n", output_dir))
cat(sprintf("  6. %s/Table_Bidirectional_MR_TRIG_HDL_RATIO_Significant.csv\n", output_dir))

cat("\nSummary tables:\n")
cat(sprintf("  7. %s/Table_Bidirectional_MR_Summary.csv\n", output_dir))
cat(sprintf("  8. %s/Table_Bidirectional_MR_Both_Sig.csv\n", output_dir))
cat(sprintf("  9. %s/Table_Bidirectional_MR_Fwd_Only_Sig.csv\n", output_dir))

cat("\n================================================================================\n")
cat("                          DONE\n")
cat("================================================================================\n")
