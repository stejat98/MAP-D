#!/usr/bin/env Rscript
# ============================================================================
# Compare LEPR MR estimates: BMI Normoglycemic vs BMI Full vs Observational
# ============================================================================

cat("================================================================\n")
cat("LEPR -> BMI: MR (Normoglycemic vs Full) vs Observational (UKB)\n")
cat("================================================================\n\n")

# ----------------------------------------------------------
# 1. LEPR pQTL instruments (exposure: protein levels)
#    From Sun et al. 2023 (UKB-PPP), via ST10
# ----------------------------------------------------------
lepr_snps <- data.frame(
  SNP = c("rs1260326", "rs13232120", "rs1801689", "rs186021206", "rs2049865",
          "rs2298475", "rs2393775", "rs28929474", "rs445810", "rs61779759", "rs7287617"),
  effect_allele = c("C", "T", "C", "A", "A", "C", "A", "T", "C", "A", "A"),
  other_allele  = c("T", "A", "A", "G", "C", "T", "G", "C", "T", "G", "G"),
  beta_exp = c(-0.1425720, -0.0751232, 0.2658970, 0.4218220, 0.0424502,
               -0.0841103, -0.0543263, 0.3104240, 0.0459359, -0.2702120, 0.0468152),
  se_exp   = c(0.00616048, 0.00913617, 0.01747750, 0.04227510, 0.00615188,
               0.01139590, 0.00613681, 0.02158090, 0.00677751, 0.00798028, 0.00639830),
  pval_exp = c(1.71e-118, 1.99e-16, 2.87e-52, 1.90e-23, 5.19e-12,
               1.57e-13, 8.56e-19, 6.50e-47, 1.22e-11, 1e-200, 2.54e-13),
  chr = c(2, 7, 17, 17, 8, 11, 12, 14, 5, 1, 22),
  pos = c(27730940, 72983310, 64210580, 7069412, 116588546,
          126278203, 121424574, 94844947, 81731077, 66161965, 39897163),
  stringsAsFactors = FALSE
)

# LEPR gene: chr1:65,420,652-65,641,559 (GRCh37)
lepr_start <- 65420652; lepr_end <- 65641559; cis_window <- 1000000
lepr_snps$is_cis <- lepr_snps$chr == 1 &
  lepr_snps$pos >= (lepr_start - cis_window) &
  lepr_snps$pos <= (lepr_end + cis_window)

cat(sprintf("LEPR instruments: %d total (%d cis, %d trans)\n\n",
            nrow(lepr_snps), sum(lepr_snps$is_cis), sum(!lepr_snps$is_cis)))

# ----------------------------------------------------------
# 2. BMI outcome data: NORMOGLYCEMIC stratum
#    GWAS columns: CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE ...
#    BETA is for ALLELE1 (the effect allele), NOT ALLELE0
# ----------------------------------------------------------
normo <- data.frame(
  SNP           = c("rs61779759","rs1260326","rs445810","rs13232120","rs2049865",
                     "rs2298475","rs2393775","rs28929474","rs186021206","rs1801689","rs7287617"),
  ea_out        = c("G","T","T","A","C","T","G","C","G","A","G"),   # ALLELE1 (effect)
  oa_out        = c("A","C","C","T","A","C","A","T","A","C","A"),   # ALLELE0 (other)
  beta_out_norm = c(0.0248193,-0.0196274,0.0215351,-0.0415159,0.0284892,
                    -0.0109624,-0.0172886,0.110826,0.000295637,-0.0137909,-0.0285189),
  se_out_norm   = c(0.0154103,0.0118601,0.0132876,0.0176687,0.0119119,
                    0.0220195,0.0119901,0.0413688,0.08906,0.0337728,0.0124125),
  n_norm        = c(275715,276279,271242,276278,269436,276279,276250,276279,274011,276279,267367),
  stringsAsFactors = FALSE
)

# ----------------------------------------------------------
# 3. BMI outcome data: FULL stratum
# ----------------------------------------------------------
full <- data.frame(
  SNP           = c("rs61779759","rs1260326","rs445810","rs13232120","rs2049865",
                     "rs2298475","rs2393775","rs28929474","rs186021206","rs1801689","rs7287617"),
  ea_out_f      = c("G","T","T","A","C","T","G","C","G","A","G"),   # ALLELE1 (effect)
  oa_out_f      = c("A","C","C","T","A","C","A","T","A","C","A"),   # ALLELE0 (other)
  beta_out_full = c(0.0316348,-0.060529,0.0281001,-0.0647144,0.0371286,
                    -0.00391242,-0.000300556,0.124772,-0.0243298,-0.00875199,-0.0288421),
  se_out_full   = c(0.0148411,0.0114305,0.0128055,0.0169739,0.0114733,
                    0.0212296,0.0115396,0.0400906,0.0857953,0.0324863,0.0119672),
  n_full        = c(349136,349849,343508,349847,341093,349849,349813,349849,347003,349849,338495),
  stringsAsFactors = FALSE
)

# ----------------------------------------------------------
# 4. Merge & harmonize alleles
# ----------------------------------------------------------
dat <- merge(lepr_snps, normo, by = "SNP")
dat <- merge(dat, full, by = "SNP")

# Harmonize: flip outcome beta if alleles are swapped
for (i in 1:nrow(dat)) {
  # Normoglycemic
  if (dat$effect_allele[i] != dat$ea_out[i]) {
    if (dat$effect_allele[i] == dat$oa_out[i]) {
      dat$beta_out_norm[i] <- -dat$beta_out_norm[i]
    }
  }
  # Full
  if (dat$effect_allele[i] != dat$ea_out_f[i]) {
    if (dat$effect_allele[i] == dat$oa_out_f[i]) {
      dat$beta_out_full[i] <- -dat$beta_out_full[i]
    }
  }
}

# ----------------------------------------------------------
# 5. MR helper functions
# ----------------------------------------------------------
ivw_mr <- function(bx, sx, by, sy) {
  w <- 1 / sy^2
  b <- sum(w * by / bx) / sum(w / bx^2)
  se <- sqrt(1 / sum(w / bx^2))
  p <- 2 * pnorm(-abs(b / se))
  bsnp <- by / bx
  Q <- sum(w * (bsnp - b)^2)
  Q_df <- length(bx) - 1
  Q_p <- pchisq(Q, Q_df, lower.tail = FALSE)
  list(beta = b, se = se, pval = p, Q = Q, Q_df = Q_df, Q_pval = Q_p)
}

wald_ratio <- function(bx, sx, by, sy) {
  b <- by / bx
  se <- abs(sy / bx)
  p <- 2 * pnorm(-abs(b / se))
  list(beta = b, se = se, pval = p)
}

weighted_median_mr <- function(bx, sx, by, sy, nboot = 1000) {
  bsnp <- by / bx
  se_snp <- abs(sy / bx)
  w <- 1 / se_snp^2
  w <- w / sum(w)
  o <- order(bsnp)
  bsnp_o <- bsnp[o]; w_o <- w[o]
  b <- bsnp_o[which(cumsum(w_o) >= 0.5)[1]]
  set.seed(42)
  boots <- replicate(nboot, {
    bb <- rnorm(length(bsnp), bsnp, se_snp)
    oo <- order(bb); bb <- bb[oo]; ww <- w[oo]
    bb[which(cumsum(ww) >= 0.5)[1]]
  })
  se <- sd(boots)
  p <- 2 * pnorm(-abs(b / se))
  list(beta = b, se = se, pval = p)
}

fmt_p <- function(p) ifelse(p < 0.001, sprintf("%.2e", p), sprintf("%.4f", p))
fmt_sig <- function(p) ifelse(p < 0.05, " *", "")

# ----------------------------------------------------------
# 6. Run MR: ALL pQTLs (IVW)
# ----------------------------------------------------------
cat("================================================================\n")
cat("A. ALL pQTLs IVW (11 SNPs)\n")
cat("================================================================\n\n")

ivw_norm <- ivw_mr(dat$beta_exp, dat$se_exp, dat$beta_out_norm, dat$se_out_norm)
ivw_full <- ivw_mr(dat$beta_exp, dat$se_exp, dat$beta_out_full, dat$se_out_full)
wm_norm  <- weighted_median_mr(dat$beta_exp, dat$se_exp, dat$beta_out_norm, dat$se_out_norm)
wm_full  <- weighted_median_mr(dat$beta_exp, dat$se_exp, dat$beta_out_full, dat$se_out_full)

cat(sprintf("%-25s | %-14s | %5s | %8s | %8s | %12s\n",
            "Method", "Stratum", "nSNP", "Beta", "SE", "P-value"))
cat(paste(rep("-", 80), collapse=""), "\n")

cat(sprintf("%-25s | %-14s | %5d | %8.4f | %8.4f | %12s%s\n",
            "IVW", "Normoglycemic", 11, ivw_norm$beta, ivw_norm$se, fmt_p(ivw_norm$pval), fmt_sig(ivw_norm$pval)))
cat(sprintf("%-25s | %-14s | %5d | %8.4f | %8.4f | %12s%s\n",
            "IVW", "Full", 11, ivw_full$beta, ivw_full$se, fmt_p(ivw_full$pval), fmt_sig(ivw_full$pval)))
cat(sprintf("%-25s | %-14s | %5d | %8.4f | %8.4f | %12s%s\n",
            "Weighted median", "Normoglycemic", 11, wm_norm$beta, wm_norm$se, fmt_p(wm_norm$pval), fmt_sig(wm_norm$pval)))
cat(sprintf("%-25s | %-14s | %5d | %8.4f | %8.4f | %12s%s\n",
            "Weighted median", "Full", 11, wm_full$beta, wm_full$se, fmt_p(wm_full$pval), fmt_sig(wm_full$pval)))

cat(sprintf("\nHeterogeneity (Normoglycemic): Q=%.2f, df=%d, P=%.4f\n",
            ivw_norm$Q, ivw_norm$Q_df, ivw_norm$Q_pval))
cat(sprintf("Heterogeneity (Full):          Q=%.2f, df=%d, P=%.4f\n\n",
            ivw_full$Q, ivw_full$Q_df, ivw_full$Q_pval))

# ----------------------------------------------------------
# 7. Run MR: CIS-only (Wald ratio, 1 SNP: rs61779759)
# ----------------------------------------------------------
cat("================================================================\n")
cat("B. CIS-ONLY pQTL (Wald ratio: rs61779759)\n")
cat("================================================================\n\n")

dc <- dat[dat$is_cis, ]
cis_norm <- wald_ratio(dc$beta_exp, dc$se_exp, dc$beta_out_norm, dc$se_out_norm)
cis_full <- wald_ratio(dc$beta_exp, dc$se_exp, dc$beta_out_full, dc$se_out_full)

cat(sprintf("%-14s | %8s | %8s | %12s | %8s | %18s\n",
            "Stratum", "Beta", "SE", "P-value", "OR", "95% CI"))
cat(paste(rep("-", 80), collapse=""), "\n")

for (res_list in list(
  list("Normoglycemic", cis_norm, dc$n_norm[1]),
  list("Full", cis_full, dc$n_full[1])
)) {
  nm <- res_list[[1]]; r <- res_list[[2]]; n <- res_list[[3]]
  or <- exp(r$beta); lo <- exp(r$beta - 1.96*r$se); hi <- exp(r$beta + 1.96*r$se)
  cat(sprintf("%-14s | %8.4f | %8.4f | %12s | %8.3f | (%.3f, %.3f)%s  [N=%d]\n",
              nm, r$beta, r$se, fmt_p(r$pval), or, lo, hi, fmt_sig(r$pval), n))
}

# ----------------------------------------------------------
# 8. UKB Observational estimates (from merged_step1_2_ukb_T2D_coloc.csv)
# ----------------------------------------------------------
cat("\n\n================================================================\n")
cat("C. UKB OBSERVATIONAL: LEPR ~ BMI (protein-phenotype association)\n")
cat("================================================================\n\n")

obs <- data.frame(
  Stratum = c("Normoglycemic (non-T2D)", "Prediabetes", "T2D"),
  Estimate_UKB = c(-0.2817, -0.2919, -0.2299),
  P_UKB = c(0, 1.30e-102, 5.24e-13),
  stringsAsFactors = FALSE
)

cat(sprintf("%-25s | %10s | %12s\n", "Stratum", "Beta (obs)", "P-value"))
cat(paste(rep("-", 55), collapse=""), "\n")
for (i in 1:nrow(obs)) {
  pstr <- if (obs$P_UKB[i] == 0) "< 1e-300" else sprintf("%.2e", obs$P_UKB[i])
  cat(sprintf("%-25s | %10.4f | %12s\n", obs$Stratum[i], obs$Estimate_UKB[i], pstr))
}

cat("\nNote: Observational estimates are standardized regression coefficients\n")
cat("      of LEPR protein levels on BMI from UKB (Estimate_UKB column).\n")
cat("      Negative = higher LEPR levels associated with LOWER BMI.\n")

# ----------------------------------------------------------
# 9. GRAND COMPARISON TABLE
# ----------------------------------------------------------
cat("\n\n================================================================\n")
cat("D. GRAND COMPARISON: LEPR -> BMI\n")
cat("================================================================\n\n")

cat(sprintf("%-30s | %-15s | %8s | %8s | %12s | %5s\n",
            "Analysis", "Stratum", "Beta", "SE", "P-value", "Sig"))
cat(paste(rep("-", 90), collapse=""), "\n")

# MR: All pQTLs IVW
cat(sprintf("%-30s | %-15s | %8.4f | %8.4f | %12s | %5s\n",
            "MR (IVW, 11 SNPs)", "Normoglycemic", ivw_norm$beta, ivw_norm$se,
            fmt_p(ivw_norm$pval), fmt_sig(ivw_norm$pval)))
cat(sprintf("%-30s | %-15s | %8.4f | %8.4f | %12s | %5s\n",
            "MR (IVW, 11 SNPs)", "Full", ivw_full$beta, ivw_full$se,
            fmt_p(ivw_full$pval), fmt_sig(ivw_full$pval)))

cat(paste(rep("-", 90), collapse=""), "\n")

# MR: Cis-only Wald
cat(sprintf("%-30s | %-15s | %8.4f | %8.4f | %12s | %5s\n",
            "MR (cis Wald, 1 SNP)", "Normoglycemic", cis_norm$beta, cis_norm$se,
            fmt_p(cis_norm$pval), fmt_sig(cis_norm$pval)))
cat(sprintf("%-30s | %-15s | %8.4f | %8.4f | %12s | %5s\n",
            "MR (cis Wald, 1 SNP)", "Full", cis_full$beta, cis_full$se,
            fmt_p(cis_full$pval), fmt_sig(cis_full$pval)))

cat(paste(rep("-", 90), collapse=""), "\n")

# MR: Weighted Median
cat(sprintf("%-30s | %-15s | %8.4f | %8.4f | %12s | %5s\n",
            "MR (Wt. Median, 11 SNPs)", "Normoglycemic", wm_norm$beta, wm_norm$se,
            fmt_p(wm_norm$pval), fmt_sig(wm_norm$pval)))
cat(sprintf("%-30s | %-15s | %8.4f | %8.4f | %12s | %5s\n",
            "MR (Wt. Median, 11 SNPs)", "Full", wm_full$beta, wm_full$se,
            fmt_p(wm_full$pval), fmt_sig(wm_full$pval)))

cat(paste(rep("-", 90), collapse=""), "\n")

# Observational
cat(sprintf("%-30s | %-15s | %8.4f | %8s | %12s | %5s\n",
            "Observational (UKB)", "Non-T2D", obs$Estimate_UKB[1], "    --",
            "< 1e-300", "  ***"))
cat(sprintf("%-30s | %-15s | %8.4f | %8s | %12s | %5s\n",
            "Observational (UKB)", "Prediabetes", obs$Estimate_UKB[2], "    --",
            "1.30e-102", "  ***"))
cat(sprintf("%-30s | %-15s | %8.4f | %8s | %12s | %5s\n",
            "Observational (UKB)", "T2D", obs$Estimate_UKB[3], "    --",
            "5.24e-13", "  ***"))

# ----------------------------------------------------------
# 10. Key interpretation
# ----------------------------------------------------------
cat("\n\n================================================================\n")
cat("INTERPRETATION\n")
cat("================================================================\n\n")

cat("1. OBSERVATIONAL (UKB): LEPR is strongly NEGATIVELY associated with BMI\n")
cat("   across all strata (beta ~ -0.23 to -0.29, P < 1e-13).\n")
cat("   Higher LEPR protein -> lower BMI (consistent with leptin signaling).\n\n")

cat("2. MR (Full strata):\n")
cat(sprintf("   - IVW (11 SNPs):  beta = %.4f, P = %s\n",
            ivw_full$beta, fmt_p(ivw_full$pval)))
cat(sprintf("   - Cis Wald:       beta = %.4f, P = %s%s\n",
            cis_full$beta, fmt_p(cis_full$pval), fmt_sig(cis_full$pval)))

cat("\n3. MR (Normoglycemic):\n")
cat(sprintf("   - IVW (11 SNPs):  beta = %.4f, P = %s\n",
            ivw_norm$beta, fmt_p(ivw_norm$pval)))
cat(sprintf("   - Cis Wald:       beta = %.4f, P = %s\n",
            cis_norm$beta, fmt_p(cis_norm$pval)))

# Direction comparison
cat("\n4. DIRECTION COMPARISON:\n")
cat(sprintf("   Observational (non-T2D):   beta = %.4f  (NEGATIVE, highly sig)\n", obs$Estimate_UKB[1]))
cat(sprintf("   MR cis (Normoglycemic):    beta = %.4f  (%s)\n",
            cis_norm$beta, ifelse(cis_norm$beta > 0, "POSITIVE", "NEGATIVE")))
cat(sprintf("   MR cis (Full):             beta = %.4f  (%s, P=%s)\n",
            cis_full$beta, ifelse(cis_full$beta > 0, "POSITIVE", "NEGATIVE"), fmt_p(cis_full$pval)))

# Concordance with observational
cat("\n5. CONCORDANCE WITH OBSERVATIONAL:\n")
conc_norm <- sign(cis_norm$beta) == sign(obs$Estimate_UKB[1])
conc_full <- sign(cis_full$beta) == sign(obs$Estimate_UKB[1])
cat(sprintf("   Cis-MR Normoglycemic vs Obs: %s\n", ifelse(conc_norm, "CONCORDANT", "DISCORDANT")))
cat(sprintf("   Cis-MR Full vs Obs:          %s\n", ifelse(conc_full, "CONCORDANT", "DISCORDANT")))

if (!conc_full) {
  cat("\n   => DISCORDANCE: Observational = NEGATIVE, cis-MR Full = POSITIVE.\n")
  cat("      The observational LEPR-BMI association is likely confounded by\n")
  cat("      reverse causation (adiposity drives LEPR levels upward via leptin).\n")
  cat("      The MR cis estimate captures the true causal direction.\n")
}

cat(sprintf("\n6. NORMOGLYCEMIC vs FULL cis-MR:\n"))
if (sign(cis_norm$beta) == sign(cis_full$beta)) {
  cat(sprintf("   Same direction: Normo (%.4f) and Full (%.4f)\n",
              cis_norm$beta, cis_full$beta))
} else {
  cat(sprintf("   DIRECTION CHANGE: Normo (%.4f) vs Full (%.4f)\n",
              cis_norm$beta, cis_full$beta))
}

# Implication for quadrant plot
cat("\n7. QUADRANT PLOT IMPLICATION:\n")
cat("   The concordant plot (validation_step1_points_BMI_concordant.pdf)\n")
cat("   shows cis-MR hits as triangles only if sign(Beta_cisMR) == sign(Estimate_UKB).\n")
cat(sprintf("   LEPR: Beta_cisMR = %.4f (POSITIVE), Estimate_UKB = %.4f (NEGATIVE)\n",
            cis_full$beta, obs$Estimate_UKB[1]))
cat(sprintf("   => LEPR is %s with observational estimate.\n",
            ifelse(conc_full, "CONCORDANT (shown as triangle)", "DISCORDANT (not shown as triangle)")))
cat("   This is CORRECT: the MR causal estimate is opposite to the observational.\n")

cat("\n================================================================\n")
cat("Done.\n")
cat("================================================================\n")
