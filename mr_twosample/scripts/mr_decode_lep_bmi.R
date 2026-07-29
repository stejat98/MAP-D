#!/usr/bin/env Rscript
# =============================================================================
# Mendelian Randomization: DECODE cis-pQTLs for LEP -> BMI (UKB full stratum)
#
# Exposure: Leptin (LEP) cis-pQTLs from DECODE Iceland (SomaScan)
#   File: Proteomics_SMP_PC0_8484_24_LEP_Leptin_10032022.txt.gz
#   LEP gene: chr7:128,241,278-128,257,629 (hg38)
#   Cis window: +/- 1Mb
#
# Outcome: BMI full stratum (UKB regenie pipeline)
#   File: BMI_full_all_chr.regenie.gz
# =============================================================================

.libPaths("/path/to/project/.R/library")

library(TwoSampleMR)
library(data.table)
library(dplyr)
library(ggplot2)

cat("================================================================================\n")
cat("  MR Analysis: DECODE cis-pQTLs for LEP -> BMI (UKB full)\n")
cat("================================================================================\n\n")

setwd("/path/to/project/regenie_pipeline")

output_dir <- "results/twosampleMR/decode_lep_bmi"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. Load DECODE LEP pQTL data - extract cis region
# =============================================================================
cat("Step 1: Loading DECODE LEP pQTL data (cis region)...\n")

decode_file <- "/path/to/shared_data/DECODE/pQTL/final_somascan_smp/Proteomics_SMP_PC0_8484_24_LEP_Leptin_10032022.txt.gz"

# LEP gene coordinates (from somascan_protein_map.tsv)
lep_chr <- "chr7"
lep_start <- 128241278
lep_end   <- 128257629
cis_window <- 1e6  # 1 Mb

cis_start <- lep_start - cis_window
cis_end   <- lep_end + cis_window

cat(sprintf("  LEP gene: %s:%d-%d\n", lep_chr, lep_start, lep_end))
cat(sprintf("  Cis window: %s:%d-%d (+/- %g Mb)\n", lep_chr, cis_start, cis_end, cis_window / 1e6))

# Read full file and filter to cis region
# Using awk for efficiency given large file size
cmd <- sprintf(
  "zcat %s | awk 'NR==1 || ($1==\"chr7\" && $2>=%d && $2<=%d)'",
  decode_file, cis_start, cis_end
)
cis_data <- fread(cmd = cmd)

cat(sprintf("  Total variants in cis region: %d\n", nrow(cis_data)))

# Filter to genome-wide significant
cis_sig <- cis_data[Pval < 5e-8, ]
cat(sprintf("  Genome-wide significant (P < 5e-8): %d\n", nrow(cis_sig)))

if (nrow(cis_sig) == 0) {
  stop("No genome-wide significant cis-pQTLs found for LEP in DECODE data!")
}

# Show top hits
cat("\n  Top 10 cis-pQTLs by p-value:\n")
top_hits <- cis_sig[order(Pval)][1:min(10, nrow(cis_sig)),
  .(rsids, Chrom, Pos, effectAllele, otherAllele, Beta, SE, Pval, ImpMAF)]
print(top_hits)
cat("\n")

# =============================================================================
# 2. Prepare exposure data for TwoSampleMR
# =============================================================================
cat("Step 2: Preparing exposure data...\n")

# Filter to variants with rsIDs (needed for clumping and matching)
cis_sig_rs <- cis_sig[grepl("^rs", rsids), ]
cat(sprintf("  Significant cis-pQTLs with rsIDs: %d\n", nrow(cis_sig_rs)))

# For variants with multiple rsIDs, take the first one
cis_sig_rs[, SNP := sapply(strsplit(rsids, ","), function(x) x[1])]

# Format as exposure data
exposure_raw <- data.frame(
  SNP = cis_sig_rs$SNP,
  beta = cis_sig_rs$Beta,
  se = cis_sig_rs$SE,
  effect_allele = cis_sig_rs$effectAllele,
  other_allele = cis_sig_rs$otherAllele,
  eaf = cis_sig_rs$ImpMAF,
  pval = cis_sig_rs$Pval,
  samplesize = cis_sig_rs$N,
  exposure = "LEP_DECODE_cis",
  id.exposure = "LEP_DECODE_cis",
  chr = gsub("chr", "", cis_sig_rs$Chrom),
  pos = cis_sig_rs$Pos,
  stringsAsFactors = FALSE
)

exposure_dat <- format_data(
  exposure_raw,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  samplesize_col = "samplesize",
  phenotype_col = "exposure",
  chr_col = "chr",
  pos_col = "pos"
)

cat(sprintf("  Formatted exposure SNPs: %d\n", nrow(exposure_dat)))

# =============================================================================
# 3. Clump instruments (LD-based pruning)
# =============================================================================
cat("\nStep 3: Clumping instruments (r2 < 0.001, window = 10Mb)...\n")

# Try clumping with TwoSampleMR (uses 1000 Genomes EUR LD reference)
tryCatch({
  exposure_clumped <- clump_data(
    exposure_dat,
    clump_r2 = 0.001,
    clump_kb = 10000,
    pop = "EUR"
  )
  cat(sprintf("  Independent instruments after clumping: %d\n", nrow(exposure_clumped)))
}, error = function(e) {
  cat(sprintf("  Clumping via API failed: %s\n", e$message))
  cat("  Falling back to local p-value-based pruning (keeping top hit per 500kb window)...\n")

  # Manual pruning: sort by p-value, keep non-overlapping 500kb windows
  exposure_dat_sorted <- exposure_dat[order(exposure_dat$pval.exposure), ]
  keep <- logical(nrow(exposure_dat_sorted))
  keep[1] <- TRUE
  kept_positions <- exposure_dat_sorted$pos.exposure[1]

  for (i in 2:nrow(exposure_dat_sorted)) {
    cur_pos <- exposure_dat_sorted$pos.exposure[i]
    if (all(abs(cur_pos - kept_positions) > 500000)) {
      keep[i] <- TRUE
      kept_positions <- c(kept_positions, cur_pos)
    }
  }

  exposure_clumped <<- exposure_dat_sorted[keep, ]
  cat(sprintf("  Independent instruments after distance pruning: %d\n", nrow(exposure_clumped)))
})

cat("\n  Instruments used:\n")
print(exposure_clumped[, c("SNP", "beta.exposure", "se.exposure", "pval.exposure",
                           "effect_allele.exposure", "other_allele.exposure", "eaf.exposure")])
cat("\n")

# =============================================================================
# 4. Load BMI full outcome data
# =============================================================================
cat("Step 4: Loading BMI full GWAS outcome data...\n")

bmi_file <- "results/GWAS/BMI/BMI_full_all_chr.regenie.gz"

# Only extract the SNPs we need
snp_list <- exposure_clumped$SNP
snp_pattern <- paste(snp_list, collapse = "|")

cmd_bmi <- sprintf(
  "zcat %s | head -1 && zcat %s | grep -wE '%s'",
  bmi_file, bmi_file, snp_pattern
)

bmi_gwas <- tryCatch({
  fread(cmd = cmd_bmi)
}, error = function(e) {
  cat("  Direct grep failed, trying full scan...\n")
  bmi_full <- fread(cmd = sprintf("zcat %s", bmi_file))
  bmi_full[ID %in% snp_list, ]
})

cat(sprintf("  BMI variants matching exposure instruments: %d\n", nrow(bmi_gwas)))

# Format as outcome
outcome_raw <- data.frame(
  SNP = bmi_gwas$ID,
  beta = as.numeric(bmi_gwas$BETA),
  se = as.numeric(bmi_gwas$SE),
  effect_allele = bmi_gwas$ALLELE1,
  other_allele = bmi_gwas$ALLELE0,
  eaf = as.numeric(bmi_gwas$A1FREQ),
  pval = 10^(-as.numeric(bmi_gwas$LOG10P)),
  samplesize = as.numeric(bmi_gwas$N),
  outcome = "BMI_full_UKB",
  id.outcome = "BMI_full_UKB",
  stringsAsFactors = FALSE
)

outcome_dat <- format_data(
  outcome_raw,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  samplesize_col = "samplesize",
  phenotype_col = "outcome"
)

cat(sprintf("  Formatted outcome SNPs: %d\n", nrow(outcome_dat)))

# =============================================================================
# 5. Harmonize exposure and outcome
# =============================================================================
cat("\nStep 5: Harmonizing exposure and outcome data...\n")

dat_harmonized <- harmonise_data(
  exposure_dat = exposure_clumped,
  outcome_dat = outcome_dat,
  action = 2
)

dat_harmonized <- dat_harmonized[dat_harmonized$mr_keep, ]
cat(sprintf("  Harmonized SNPs kept for MR: %d\n", nrow(dat_harmonized)))

if (nrow(dat_harmonized) == 0) {
  stop("No harmonized SNPs available for MR analysis!")
}

# Save harmonized data
write.csv(dat_harmonized,
  file.path(output_dir, "harmonized_decode_LEP_BMI_full.csv"),
  row.names = FALSE
)

cat("\n  Harmonized data:\n")
print(dat_harmonized[, c("SNP", "effect_allele.exposure", "other_allele.exposure",
                         "beta.exposure", "se.exposure", "pval.exposure",
                         "beta.outcome", "se.outcome", "pval.outcome")])

# =============================================================================
# 6. Run MR analysis
# =============================================================================
cat("\nStep 6: Running MR analysis...\n")

if (nrow(dat_harmonized) == 1) {
  cat("  Single instrument -> Wald ratio\n")
  mr_results <- mr(dat_harmonized, method_list = c("mr_wald_ratio"))
} else if (nrow(dat_harmonized) == 2) {
  cat("  Two instruments -> IVW and Wald ratio\n")
  mr_results <- mr(dat_harmonized, method_list = c("mr_ivw", "mr_wald_ratio"))
} else {
  cat("  Multiple instruments -> Full MR methods\n")
  mr_results <- mr(dat_harmonized, method_list = c(
    "mr_ivw",
    "mr_egger_regression",
    "mr_weighted_median",
    "mr_weighted_mode",
    "mr_wald_ratio"
  ))
}

cat("\n================================================================================\n")
cat("                          MR RESULTS\n")
cat("================================================================================\n\n")

cat(sprintf("  Exposure: LEP (Leptin) - DECODE Iceland cis-pQTLs\n"))
cat(sprintf("  Outcome:  BMI full stratum (UKB)\n"))
cat(sprintf("  N instruments: %d\n\n", nrow(dat_harmonized)))

print(mr_results[, c("method", "nsnp", "b", "se", "pval")])
cat("\n")

# Save results
write.csv(mr_results,
  file.path(output_dir, "MR_decode_LEP_BMI_full.csv"),
  row.names = FALSE
)

# =============================================================================
# 7. Sensitivity analyses (if >1 SNP)
# =============================================================================
if (nrow(dat_harmonized) > 1) {
  cat("Step 7: Running sensitivity analyses...\n\n")

  # Heterogeneity
  cat("--- Heterogeneity ---\n")
  het <- mr_heterogeneity(dat_harmonized)
  print(het)
  write.csv(het, file.path(output_dir, "heterogeneity_decode_LEP_BMI_full.csv"),
    row.names = FALSE)
  cat("\n")

  # Pleiotropy (MR-Egger intercept)
  if (nrow(dat_harmonized) >= 3) {
    cat("--- MR-Egger Pleiotropy Test ---\n")
    pleio <- mr_pleiotropy_test(dat_harmonized)
    print(pleio)
    write.csv(pleio, file.path(output_dir, "pleiotropy_decode_LEP_BMI_full.csv"),
      row.names = FALSE)
    cat("\n")
  }

  # Leave-one-out
  cat("--- Leave-one-out ---\n")
  loo <- mr_leaveoneout(dat_harmonized)
  print(loo)
  write.csv(loo, file.path(output_dir, "leaveoneout_decode_LEP_BMI_full.csv"),
    row.names = FALSE)
  cat("\n")

  # Single SNP analysis
  cat("--- Single SNP Analysis ---\n")
  single <- mr_singlesnp(dat_harmonized)
  print(single)
  write.csv(single, file.path(output_dir, "singlesnp_decode_LEP_BMI_full.csv"),
    row.names = FALSE)

  # =============================================================================
  # 8. Generate plots
  # =============================================================================
  cat("\nStep 8: Generating plots...\n")

  # Scatter plot
  p_scatter <- mr_scatter_plot(mr_results, dat_harmonized)
  ggsave(file.path(output_dir, "scatter_decode_LEP_BMI_full.png"),
    p_scatter[[1]], width = 8, height = 6, dpi = 150)
  ggsave(file.path(output_dir, "scatter_decode_LEP_BMI_full.pdf"),
    p_scatter[[1]], width = 8, height = 6)
  cat("  Saved scatter plot\n")

  # Forest plot
  p_forest <- mr_forest_plot(single)
  ggsave(file.path(output_dir, "forest_decode_LEP_BMI_full.png"),
    p_forest[[1]], width = 10, height = 6, dpi = 150)
  ggsave(file.path(output_dir, "forest_decode_LEP_BMI_full.pdf"),
    p_forest[[1]], width = 10, height = 6)
  cat("  Saved forest plot\n")

  # Leave-one-out plot
  p_loo <- mr_leaveoneout_plot(loo)
  ggsave(file.path(output_dir, "leaveoneout_decode_LEP_BMI_full.png"),
    p_loo[[1]], width = 10, height = 6, dpi = 150)
  ggsave(file.path(output_dir, "leaveoneout_decode_LEP_BMI_full.pdf"),
    p_loo[[1]], width = 10, height = 6)
  cat("  Saved leave-one-out plot\n")

  # Funnel plot
  p_funnel <- mr_funnel_plot(single)
  ggsave(file.path(output_dir, "funnel_decode_LEP_BMI_full.png"),
    p_funnel[[1]], width = 8, height = 6, dpi = 150)
  ggsave(file.path(output_dir, "funnel_decode_LEP_BMI_full.pdf"),
    p_funnel[[1]], width = 8, height = 6)
  cat("  Saved funnel plot\n")

} else {
  cat("\nStep 7: Single instrument - computing Wald ratio details...\n\n")

  snp <- dat_harmonized$SNP[1]
  beta_exp <- dat_harmonized$beta.exposure[1]
  se_exp   <- dat_harmonized$se.exposure[1]
  beta_out <- dat_harmonized$beta.outcome[1]
  se_out   <- dat_harmonized$se.outcome[1]

  wald_beta <- beta_out / beta_exp
  wald_se   <- abs(se_out / beta_exp)
  wald_z    <- wald_beta / wald_se
  wald_p    <- 2 * pnorm(-abs(wald_z))

  cat(sprintf("  Instrument: %s\n", snp))
  cat(sprintf("  Exposure (LEP):  beta = %.4f, SE = %.4f\n", beta_exp, se_exp))
  cat(sprintf("  Outcome  (BMI):  beta = %.4f, SE = %.4f\n", beta_out, se_out))
  cat(sprintf("  Wald ratio:      beta = %.4f, SE = %.4f, Z = %.3f, P = %.4e\n",
              wald_beta, wald_se, wald_z, wald_p))

  # F-statistic for instrument strength
  f_stat <- (beta_exp / se_exp)^2
  cat(sprintf("  F-statistic:     %.1f %s\n", f_stat,
              ifelse(f_stat >= 10, "(strong instrument)", "(WEAK instrument!)")))
}

# =============================================================================
# 9. Summary
# =============================================================================
cat("\n================================================================================\n")
cat("                          SUMMARY\n")
cat("================================================================================\n\n")

cat("Data sources:\n")
cat("  Exposure: DECODE Iceland SomaScan pQTLs (N = 35,692)\n")
cat("  Outcome:  UKB BMI full stratum (regenie)\n")
cat(sprintf("  Cis window: chr7:%d-%d (+/- 1Mb around LEP)\n", cis_start, cis_end))
cat(sprintf("  N instruments: %d\n\n", nrow(dat_harmonized)))

cat("Key MR result (IVW or Wald ratio):\n")
primary <- mr_results[mr_results$method %in% c("Inverse variance weighted", "Wald ratio"), ]
for (i in 1:nrow(primary)) {
  cat(sprintf("  %s: beta = %.4f, SE = %.4f, P = %.4e (nsnp = %d)\n",
              primary$method[i], primary$b[i], primary$se[i], primary$pval[i], primary$nsnp[i]))
}

cat(sprintf("\nResults saved to: %s/\n", output_dir))

cat("\n================================================================================\n")
cat("                          DONE\n")
cat("================================================================================\n")
