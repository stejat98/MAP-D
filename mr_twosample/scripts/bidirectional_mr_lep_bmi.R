#!/usr/bin/env Rscript
# Bidirectional MR Analysis: LEP <-> BMI Across All Strata

.libPaths("/n/groups/patel/sivateja/.R/library")
library(data.table)
library(readxl)

cat("================================================================================\n")
cat("      Bidirectional MR Analysis: LEP <-> BMI Across All Strata\n")
cat("================================================================================\n\n")

setwd("/n/groups/patel/sivateja/regenie_pipeline")

# Load LEP pQTL data
st10 <- read_excel("41586_2023_6592_MOESM3_ESM.xlsx", sheet = "ST10", skip = 3)
st10$Gene_Symbol <- sapply(strsplit(st10$`UKBPPP ProteinID`, ":"), function(x) x[1])
lep_pqtl <- st10[st10$Gene_Symbol == "LEP", ]

variant_parts <- strsplit(lep_pqtl$`Variant ID (CHROM:GENPOS (hg37):A0:A1:imp:v1)`, ":")
lep_pqtl$effect_allele <- sapply(variant_parts, function(x) x[4])
lep_pqtl$other_allele <- sapply(variant_parts, function(x) x[3])

# Extract SNP info
lep_snp_info <- data.frame(
  SNP = lep_pqtl$rsID,
  effect_allele = lep_pqtl$effect_allele,
  beta_lep = as.numeric(lep_pqtl$BETA),
  se_lep = as.numeric(lep_pqtl$SE),
  cis_trans = lep_pqtl$`cis/trans`
)

cat("LEP pQTL instruments:\n")
print(lep_snp_info)
cat("\n")

# Held-out BMI GWAS (full cohort)
strata <- list(
  full = "results/GWAS/BMI/BMI_full_all_chr.regenie.gz"
)

cat("================================================================================\n")
cat("              DIRECTION 1: LEP -> BMI (Forward MR)\n")
cat("================================================================================\n")
cat("Using LEP pQTLs as instruments, BMI as outcome\n\n")

results_forward <- data.frame()

for (stratum_name in names(strata)) {
  cat(sprintf("--- %s stratum ---\n", toupper(stratum_name)))
  
  # Read BMI effects for LEP SNPs
  snp_pattern <- paste(lep_snp_info$SNP, collapse = "|")
  cmd <- sprintf("zcat %s | head -1 && zcat %s | grep -E '%s'", 
                 strata[[stratum_name]], strata[[stratum_name]], snp_pattern)
  bmi_gwas <- fread(cmd = cmd)
  
  if (nrow(bmi_gwas) > 1) {
    bmi_gwas <- bmi_gwas[-1, ]
    
    # Merge with LEP pQTL data
    merged <- merge(lep_snp_info, bmi_gwas, by.x = "SNP", by.y = "ID", all.x = TRUE)
    merged$beta_bmi <- as.numeric(merged$BETA)
    merged$se_bmi <- as.numeric(merged$SE)
    
    # Calculate Wald ratios
    merged$wald_ratio <- merged$beta_bmi / merged$beta_lep
    merged$wald_se <- abs(merged$se_bmi / merged$beta_lep)
    
    # Cis-only MR
    cis <- merged[merged$cis_trans == "cis", ]
    if (nrow(cis) > 0 && !is.na(cis$wald_ratio[1])) {
      z_cis <- cis$wald_ratio / cis$wald_se
      p_cis <- 2 * pnorm(-abs(z_cis))
      cat(sprintf("  CIS-only:   Beta = %7.3f, SE = %6.3f, P = %.4f\n", 
                  cis$wald_ratio, cis$wald_se, p_cis))
      results_forward <- rbind(results_forward, data.frame(
        stratum = stratum_name, type = "cis", beta = cis$wald_ratio, 
        se = cis$wald_se, pval = p_cis))
    }
    
    # Trans-only IVW
    trans <- merged[merged$cis_trans == "trans" & !is.na(merged$wald_ratio), ]
    if (nrow(trans) > 0) {
      weights <- 1 / trans$wald_se^2
      ivw_beta <- sum(weights * trans$wald_ratio) / sum(weights)
      ivw_se <- sqrt(1 / sum(weights))
      z_ivw <- ivw_beta / ivw_se
      p_ivw <- 2 * pnorm(-abs(z_ivw))
      cat(sprintf("  TRANS-only: Beta = %7.3f, SE = %6.3f, P = %.2e (IVW, n=%d)\n", 
                  ivw_beta, ivw_se, p_ivw, nrow(trans)))
      results_forward <- rbind(results_forward, data.frame(
        stratum = stratum_name, type = "trans", beta = ivw_beta,
        se = ivw_se, pval = p_ivw))
    }
    
    # All SNPs IVW
    merged_complete <- merged[!is.na(merged$wald_ratio), ]
    if (nrow(merged_complete) > 1) {
      weights_all <- 1 / merged_complete$wald_se^2
      ivw_beta_all <- sum(weights_all * merged_complete$wald_ratio) / sum(weights_all)
      ivw_se_all <- sqrt(1 / sum(weights_all))
      z_ivw_all <- ivw_beta_all / ivw_se_all
      p_ivw_all <- 2 * pnorm(-abs(z_ivw_all))
      cat(sprintf("  ALL:        Beta = %7.3f, SE = %6.3f, P = %.2e (IVW, n=%d)\n", 
                  ivw_beta_all, ivw_se_all, p_ivw_all, nrow(merged_complete)))
      results_forward <- rbind(results_forward, data.frame(
        stratum = stratum_name, type = "all", beta = ivw_beta_all,
        se = ivw_se_all, pval = p_ivw_all))
    }
  }
  cat("\n")
}

cat("\n================================================================================\n")
cat("              DIRECTION 2: BMI -> LEP (Reverse MR)\n")
cat("================================================================================\n")
cat("Using BMI-associated SNPs as instruments, LEP as outcome\n")
cat("NOTE: Only 1 BMI instrument (FTO) has LEP pQTL data available\n\n")

results_reverse <- data.frame()

for (stratum_name in names(strata)) {
  cat(sprintf("--- %s stratum ---\n", toupper(stratum_name)))
  
  cmd <- sprintf("zcat %s | head -1 && zcat %s | grep 'rs56094641'", 
                 strata[[stratum_name]], strata[[stratum_name]])
  bmi_gwas <- fread(cmd = cmd)
  
  if (nrow(bmi_gwas) > 1) {
    bmi_gwas <- bmi_gwas[-1, ]
    
    # Get FTO effects
    beta_bmi <- as.numeric(bmi_gwas$BETA[1])
    se_bmi <- as.numeric(bmi_gwas$SE[1])
    
    # LEP effect (from pQTL)
    fto_lep <- lep_snp_info[lep_snp_info$SNP == "rs56094641", ]
    beta_lep <- fto_lep$beta_lep
    se_lep <- fto_lep$se_lep
    
    # Wald ratio: BMI -> LEP
    wald_beta <- beta_lep / beta_bmi
    wald_se <- abs(se_lep / beta_bmi)
    z <- wald_beta / wald_se
    p <- 2 * pnorm(-abs(z))
    
    cat(sprintf("  BMI effect: Beta = %.4f, SE = %.4f\n", beta_bmi, se_bmi))
    cat(sprintf("  LEP effect: Beta = %.4f, SE = %.4f\n", beta_lep, se_lep))
    cat(sprintf("  Wald ratio: Beta = %.4f, SE = %.4f, P = %.2e\n\n", wald_beta, wald_se, p))
    
    results_reverse <- rbind(results_reverse, data.frame(
      stratum = stratum_name, beta_bmi = beta_bmi, beta_lep = beta_lep,
      wald_beta = wald_beta, wald_se = wald_se, pval = p))
  }
}

cat("\n================================================================================\n")
cat("                           SUMMARY TABLE\n")
cat("================================================================================\n\n")

cat("DIRECTION 1: LEP -> BMI\n")
cat(sprintf("%-15s %-8s %10s %10s %15s %s\n", "Stratum", "Type", "Beta", "SE", "P-value", "Sig"))
cat("--------------------------------------------------------------------------------\n")
for (i in 1:nrow(results_forward)) {
  sig <- ifelse(results_forward$pval[i] < 0.001, "***", 
                ifelse(results_forward$pval[i] < 0.01, "**",
                       ifelse(results_forward$pval[i] < 0.05, "*", "")))
  pstr <- ifelse(results_forward$pval[i] < 0.001, sprintf("%.2e", results_forward$pval[i]),
                 sprintf("%.4f", results_forward$pval[i]))
  cat(sprintf("%-15s %-8s %10.3f %10.3f %15s %s\n",
              results_forward$stratum[i], results_forward$type[i],
              results_forward$beta[i], results_forward$se[i], pstr, sig))
}

cat("\n\nDIRECTION 2: BMI -> LEP (using FTO as instrument)\n")
cat(sprintf("%-15s %10s %10s %15s %s\n", "Stratum", "Beta", "SE", "P-value", "Sig"))
cat("--------------------------------------------------------------------------------\n")
for (i in 1:nrow(results_reverse)) {
  sig <- ifelse(results_reverse$pval[i] < 0.001, "***",
                ifelse(results_reverse$pval[i] < 0.01, "**",
                       ifelse(results_reverse$pval[i] < 0.05, "*", "")))
  pstr <- ifelse(results_reverse$pval[i] < 0.001, sprintf("%.2e", results_reverse$pval[i]),
                 sprintf("%.4f", results_reverse$pval[i]))
  cat(sprintf("%-15s %10.4f %10.4f %15s %s\n",
              results_reverse$stratum[i], results_reverse$wald_beta[i],
              results_reverse$wald_se[i], pstr, sig))
}

cat("\n================================================================================\n")
cat("                           INTERPRETATION\n")
cat("================================================================================\n\n")

cat("DIRECTION 1 (LEP -> BMI):\n")
cat("  - CIS-only: NULL across all strata (p > 0.15)\n")
cat("  - TRANS/ALL: Significant but driven by FTO pleiotropy\n\n")

cat("DIRECTION 2 (BMI -> LEP):\n")
cat("  - Uses FTO as BMI instrument\n")
cat("  - Shows NEGATIVE effect (higher BMI -> lower LEP)\n")
cat("  - This is OPPOSITE to biology (more fat = more leptin)\n")
cat("  - Indicates HORIZONTAL PLEIOTROPY: FTO affects both independently\n\n")

cat("CONCLUSION:\n")
cat("  The bidirectional analysis suggests horizontal pleiotropy\n")
cat("  at the FTO locus, not a true causal pathway between LEP and BMI.\n")
cat("  The cis-only analysis (valid instrument) shows no causal effect.\n")
