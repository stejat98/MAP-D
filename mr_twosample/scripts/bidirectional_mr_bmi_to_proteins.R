#!/usr/bin/env Rscript
# =============================================================================
# Bidirectional MR: BMI -> Protein (LEP, LEPR, SHBG)
#
# Design: UKB BMI GWAS instruments -> DECODE pQTL outcomes
#   - No sample overlap (UKB vs Iceland)
#   - Exposure: BMI (UKB regenie, full stratum)
#   - Outcome: Protein levels (DECODE Iceland SomaScan pQTLs)
# =============================================================================

.libPaths("/n/groups/patel/sivateja/.R/library")
library(data.table)

cat("================================================================================\n")
cat("  Bidirectional MR: BMI -> LEP, LEPR, SHBG\n")
cat("  Exposure: UKB BMI GWAS (full stratum)\n")
cat("  Outcome:  DECODE Iceland pQTLs (no sample overlap)\n")
cat("================================================================================\n\n")

setwd("/n/groups/patel/sivateja/regenie_pipeline")

# =============================================================================
# 1. Extract genome-wide significant BMI instruments from UKB GWAS
# =============================================================================
cat("Step 1: Loading UKB BMI GWAS and extracting instruments...\n")

bmi_file <- "results/GWAS/BMI/BMI_full_all_chr.regenie.gz"

# Extract genome-wide significant hits
cmd <- sprintf("zcat %s | awk 'NR==1 || $12 > 7.30103'", bmi_file)  # LOG10P > 7.3 = P < 5e-8
bmi_sig <- fread(cmd = cmd)
cat(sprintf("  GW-significant BMI SNPs (P < 5e-8): %d\n", nrow(bmi_sig)))

# Filter to SNPs with rsIDs
bmi_sig <- bmi_sig[grepl("^rs", ID), ]
cat(sprintf("  With rsIDs: %d\n", nrow(bmi_sig)))

# Distance-based pruning: keep top SNP per 500kb window
bmi_sig <- bmi_sig[order(-LOG10P)]
bmi_sig$pval <- 10^(-as.numeric(bmi_sig$LOG10P))

keep <- rep(FALSE, nrow(bmi_sig))
keep[1] <- TRUE
kept_chr <- bmi_sig$CHROM[1]
kept_pos <- bmi_sig$GENPOS[1]

for (i in 2:nrow(bmi_sig)) {
  cur_chr <- bmi_sig$CHROM[i]
  cur_pos <- bmi_sig$GENPOS[i]
  
  # Check if far enough from all kept SNPs on same chromosome
  same_chr_idx <- which(keep & bmi_sig$CHROM == cur_chr)
  if (length(same_chr_idx) == 0 || 
      all(abs(cur_pos - bmi_sig$GENPOS[same_chr_idx]) > 500000)) {
    keep[i] <- TRUE
  }
}

bmi_instruments <- bmi_sig[keep, ]
cat(sprintf("  Independent instruments after 500kb pruning: %d\n\n", nrow(bmi_instruments)))

# Show top 10
cat("  Top 10 BMI instruments:\n")
print(bmi_instruments[1:min(10, nrow(bmi_instruments)), 
  .(ID, CHROM, GENPOS, ALLELE0, ALLELE1, A1FREQ, BETA, SE, LOG10P)])
cat("\n")

# Compute F-statistics
bmi_instruments$F_stat <- (as.numeric(bmi_instruments$BETA) / as.numeric(bmi_instruments$SE))^2
cat(sprintf("  F-statistic: median = %.1f, min = %.1f, max = %.1f\n",
    median(bmi_instruments$F_stat), min(bmi_instruments$F_stat), max(bmi_instruments$F_stat)))
cat(sprintf("  All F > 10: %s\n\n", ifelse(all(bmi_instruments$F_stat > 10), "YES", "NO")))

# =============================================================================
# 2. Define DECODE pQTL outcome files
# =============================================================================

proteins <- list(
  LEP = list(
    file = "/n/groups/patel/IGLOO/DECODE/pQTL/final_somascan_smp/Proteomics_SMP_PC0_8484_24_LEP_Leptin_10032022.txt.gz",
    name = "Leptin"
  ),
  LEPR = list(
    file = "/n/groups/patel/IGLOO/DECODE/pQTL/final_somascan_smp/Proteomics_SMP_PC0_5400_52_LEPR_sLeptin_R_10032022.txt.gz",
    name = "Leptin Receptor"
  ),
  SHBG = list(
    file = "/n/groups/patel/IGLOO/DECODE/pQTL/final_somascan_smp/Proteomics_SMP_PC0_4929_55_SHBG_SHBG_10032022.txt.gz",
    name = "Sex Hormone-Binding Globulin"
  )
)

# =============================================================================
# 3. For each protein, look up BMI instruments in DECODE pQTL and run MR
# =============================================================================

all_results <- data.frame()
snp_list <- bmi_instruments$ID

for (prot_name in names(proteins)) {
  prot <- proteins[[prot_name]]
  
  cat("================================================================================\n")
  cat(sprintf("  BMI -> %s (%s)\n", prot_name, prot$name))
  cat("================================================================================\n\n")
  
  # Extract matching SNPs from DECODE pQTL
  # Build grep pattern (do in chunks to avoid argument too long)
  snp_pattern <- paste(snp_list, collapse = "|")
  
  # Write SNP list to temp file for grep
  tmp_snp_file <- tempfile()
  writeLines(snp_list, tmp_snp_file)
  
  cmd <- sprintf("zcat %s | head -1 && zcat %s | grep -wFf %s", 
                 prot$file, prot$file, tmp_snp_file)
  
  pqtl <- tryCatch({
    fread(cmd = cmd)
  }, error = function(e) {
    cat(sprintf("  grep failed, trying full scan...\n"))
    full <- fread(cmd = sprintf("zcat %s", prot$file))
    rbind(full[1,][0,], full[rsids %in% snp_list, ])  # keep header structure
  })
  
  unlink(tmp_snp_file)
  
  cat(sprintf("  DECODE pQTL SNPs matching BMI instruments: %d\n", nrow(pqtl)))
  
  if (nrow(pqtl) == 0) {
    cat("  WARNING: No matching SNPs found. Trying rsID column matching...\n")
    # Some DECODE files have multiple rsIDs per row
    full_pqtl <- fread(cmd = sprintf("zcat %s", prot$file))
    # Check if any rsIDs match
    matched_idx <- which(sapply(full_pqtl$rsids, function(rs) {
      any(unlist(strsplit(rs, ",")) %in% snp_list)
    }))
    pqtl <- full_pqtl[matched_idx, ]
    cat(sprintf("  After flexible matching: %d\n", nrow(pqtl)))
  }
  
  if (nrow(pqtl) < 3) {
    cat("  Too few matching SNPs for robust MR. Skipping.\n\n")
    next
  }
  
  # Extract first rsID for matching
  pqtl$SNP <- sapply(strsplit(pqtl$rsids, ","), function(x) x[1])
  
  # Merge BMI instruments with pQTL outcomes
  merged <- merge(
    bmi_instruments[, .(SNP = ID, chr_bmi = CHROM, pos_bmi = GENPOS,
                        ea_bmi = ALLELE1, oa_bmi = ALLELE0, 
                        eaf_bmi = as.numeric(A1FREQ),
                        beta_bmi = as.numeric(BETA), se_bmi = as.numeric(SE),
                        pval_bmi = pval, F_stat)],
    pqtl[, .(SNP, ea_pqtl = effectAllele, oa_pqtl = otherAllele,
             beta_pqtl = Beta, se_pqtl = SE, pval_pqtl = Pval)],
    by = "SNP"
  )
  
  cat(sprintf("  Merged SNPs: %d\n", nrow(merged)))
  
  # Harmonize alleles
  merged$harmonized <- TRUE
  for (i in 1:nrow(merged)) {
    if (toupper(merged$ea_bmi[i]) == toupper(merged$ea_pqtl[i]) &&
        toupper(merged$oa_bmi[i]) == toupper(merged$oa_pqtl[i])) {
      # Same orientation - keep as is
    } else if (toupper(merged$ea_bmi[i]) == toupper(merged$oa_pqtl[i]) &&
               toupper(merged$oa_bmi[i]) == toupper(merged$ea_pqtl[i])) {
      # Flip pQTL effect
      merged$beta_pqtl[i] <- -merged$beta_pqtl[i]
    } else {
      # Can't harmonize (strand issues)
      merged$harmonized[i] <- FALSE
    }
  }
  
  merged <- merged[merged$harmonized == TRUE, ]
  cat(sprintf("  After allele harmonization: %d\n", nrow(merged)))
  
  if (nrow(merged) < 3) {
    cat("  Too few harmonized SNPs. Skipping.\n\n")
    next
  }
  
  # Wald ratios
  merged$wald_beta <- merged$beta_pqtl / merged$beta_bmi
  merged$wald_se <- abs(merged$se_pqtl / merged$beta_bmi)
  
  # IVW (inverse-variance weighted)
  weights <- 1 / merged$wald_se^2
  ivw_beta <- sum(weights * merged$wald_beta) / sum(weights)
  ivw_se <- sqrt(1 / sum(weights))
  ivw_z <- ivw_beta / ivw_se
  ivw_p <- 2 * pnorm(-abs(ivw_z))
  
  # Weighted median (simple approximation)
  # Sort by absolute weight
  ord <- order(abs(merged$wald_beta))
  cum_weight <- cumsum(weights[ord]) / sum(weights)
  median_idx <- which(cum_weight >= 0.5)[1]
  wm_beta <- merged$wald_beta[ord[median_idx]]
  # Bootstrap SE approximation (using 1.4826 * MAD)
  wm_se <- 1.4826 * median(abs(merged$wald_beta - median(merged$wald_beta))) / 
            sqrt(length(merged$wald_beta))
  if (wm_se == 0) wm_se <- ivw_se  # fallback
  wm_z <- wm_beta / wm_se
  wm_p <- 2 * pnorm(-abs(wm_z))
  
  # MR-Egger
  if (nrow(merged) >= 3) {
    egger_fit <- lm(merged$beta_pqtl ~ merged$beta_bmi, 
                    weights = 1/merged$se_pqtl^2)
    egger_coefs <- summary(egger_fit)$coefficients
    egger_beta <- egger_coefs[2, "Estimate"]
    egger_se <- egger_coefs[2, "Std. Error"]
    egger_p <- egger_coefs[2, "Pr(>|t|)"]
    egger_intercept <- egger_coefs[1, "Estimate"]
    egger_intercept_p <- egger_coefs[1, "Pr(>|t|)"]
  }
  
  # Cochran's Q for heterogeneity
  Q <- sum(weights * (merged$wald_beta - ivw_beta)^2)
  Q_df <- nrow(merged) - 1
  Q_p <- pchisq(Q, Q_df, lower.tail = FALSE)
  
  cat(sprintf("\n  --- IVW ---\n"))
  cat(sprintf("    Beta = %.4f, SE = %.4f, P = %.4e (N SNPs = %d)\n", 
              ivw_beta, ivw_se, ivw_p, nrow(merged)))
  cat(sprintf("    Interpretation: 1 SD increase in BMI -> %.4f SD change in %s\n",
              ivw_beta, prot_name))
  
  cat(sprintf("\n  --- MR-Egger ---\n"))
  cat(sprintf("    Slope:     Beta = %.4f, SE = %.4f, P = %.4e\n", egger_beta, egger_se, egger_p))
  cat(sprintf("    Intercept: %.4f, P = %.4f %s\n", egger_intercept, egger_intercept_p,
              ifelse(egger_intercept_p < 0.05, "(PLEIOTROPY DETECTED)", "(no evidence of pleiotropy)")))
  
  cat(sprintf("\n  --- Heterogeneity ---\n"))
  cat(sprintf("    Cochran's Q = %.2f, df = %d, P = %.4e %s\n", Q, Q_df, Q_p,
              ifelse(Q_p < 0.05, "(significant heterogeneity)", "")))
  
  cat("\n")
  
  # Store results
  all_results <- rbind(all_results, data.frame(
    Direction = "BMI -> Protein",
    Protein = prot_name,
    Method = c("IVW", "MR-Egger", "Egger intercept"),
    N_SNPs = nrow(merged),
    Beta = c(ivw_beta, egger_beta, egger_intercept),
    SE = c(ivw_se, egger_se, egger_coefs[1, "Std. Error"]),
    P = c(ivw_p, egger_p, egger_intercept_p),
    stringsAsFactors = FALSE
  ))
}

# =============================================================================
# 4. Load forward MR results for comparison
# =============================================================================
cat("\n================================================================================\n")
cat("  Loading forward MR (Protein -> BMI) results for comparison\n")
cat("================================================================================\n\n")

# From cis-pQTL MR table
mr_cis <- read.csv("results/twosampleMR/supplemental_tables/Table_MR_Cis_pQTL_BMI_Full.csv")
forward <- mr_cis[mr_cis$Protein %in% c("LEP", "LEPR", "SHBG"), ]

forward_results <- data.frame(
  Direction = "Protein -> BMI",
  Protein = forward$Protein,
  Method = "Cis-pQTL (Wald/IVW)",
  N_SNPs = forward$N.SNPs,
  Beta = as.numeric(forward$Beta),
  SE = as.numeric(forward$SE),
  P = as.numeric(forward$P.value),
  stringsAsFactors = FALSE
)

# Also load DECODE forward MR results
for (prot in c("LEP", "LEPR")) {
  fname <- sprintf("results/twosampleMR/decode_%s_bmi/MR_decode_%s_BMI_full.csv", 
                   tolower(prot), prot)
  if (file.exists(fname)) {
    decode_mr <- read.csv(fname)
    primary <- decode_mr[decode_mr$method %in% c("Inverse variance weighted", "Wald ratio"), ]
    if (nrow(primary) > 0) {
      forward_results <- rbind(forward_results, data.frame(
        Direction = "Protein -> BMI",
        Protein = prot,
        Method = paste0("DECODE cis (", primary$method[1], ")"),
        N_SNPs = primary$nsnp[1],
        Beta = primary$b[1],
        SE = primary$se[1],
        P = primary$pval[1],
        stringsAsFactors = FALSE
      ))
    }
  }
}

# =============================================================================
# 5. Print combined summary table
# =============================================================================
cat("\n================================================================================\n")
cat("         BIDIRECTIONAL MR SUMMARY: BMI <-> LEP, LEPR, SHBG\n")
cat("================================================================================\n\n")

combined <- rbind(forward_results, all_results)
combined$Beta <- round(combined$Beta, 4)
combined$SE <- round(combined$SE, 4)
combined$P_fmt <- formatC(combined$P, format = "e", digits = 3)

# Print by protein
for (prot in c("LEP", "LEPR", "SHBG")) {
  cat(sprintf("--- %s ---\n", prot))
  sub <- combined[combined$Protein == prot, ]
  cat(sprintf("%-25s %-25s %6s %8s %8s %12s\n", 
              "Direction", "Method", "N_SNP", "Beta", "SE", "P"))
  cat(paste(rep("-", 90), collapse=""), "\n")
  for (i in 1:nrow(sub)) {
    cat(sprintf("%-25s %-25s %6d %8.4f %8.4f %12s\n",
                sub$Direction[i], sub$Method[i], sub$N_SNPs[i],
                sub$Beta[i], sub$SE[i], sub$P_fmt[i]))
  }
  cat("\n")
}

# Save results
output_dir <- "results/twosampleMR/bidirectional_bmi_proteins"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(combined, file.path(output_dir, "bidirectional_MR_BMI_proteins_summary.csv"),
          row.names = FALSE)
cat(sprintf("Results saved to: %s/\n", output_dir))

cat("\n================================================================================\n")
cat("                          DONE\n")
cat("================================================================================\n")
