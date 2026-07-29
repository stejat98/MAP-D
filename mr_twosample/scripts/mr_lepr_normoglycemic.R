#!/usr/bin/env Rscript
# LEPR → BMI MR for Normoglycemic Stratum (Base R only)

cat("================================================================\n")
cat("LEPR -> BMI MR ANALYSIS: NORMOGLYCEMIC STRATUM\n")
cat("================================================================\n\n")

# LEPR instruments from existing analysis
lepr_snps <- data.frame(
  SNP = c("rs1260326", "rs13232120", "rs1801689", "rs186021206", "rs2049865",
          "rs2298475", "rs2393775", "rs28929474", "rs445810", "rs61779759", "rs7287617"),
  effect_allele = c("C", "T", "C", "A", "A", "C", "A", "T", "C", "A", "A"),
  other_allele = c("T", "A", "A", "G", "C", "T", "G", "C", "T", "G", "G"),
  beta = c(-0.1425720, -0.0751232, 0.2658970, 0.4218220, 0.0424502,
           -0.0841103, -0.0543263, 0.3104240, 0.0459359, -0.2702120, 0.0468152),
  se = c(0.00616048, 0.00913617, 0.01747750, 0.04227510, 0.00615188,
         0.01139590, 0.00613681, 0.02158090, 0.00677751, 0.00798028, 0.00639830),
  pval = c(1.71e-118, 1.99e-16, 2.87e-52, 1.90e-23, 5.19e-12,
           1.57e-13, 8.56e-19, 6.50e-47, 1.22e-11, 1e-200, 2.54e-13),
  chr = c(2, 7, 17, 17, 8, 11, 12, 14, 5, 1, 22),
  pos = c(27730940, 72983310, 64210580, 7069412, 116588546,
          126278203, 121424574, 94844947, 81731077, 66161965, 39897163),
  stringsAsFactors = FALSE
)

# LEPR gene location: chr1:65,420,652-65,641,559 (GRCh37)
lepr_start <- 65420652
lepr_end <- 65641559
cis_window <- 1000000

lepr_snps$is_cis <- lepr_snps$chr == 1 & 
                    lepr_snps$pos >= (lepr_start - cis_window) & 
                    lepr_snps$pos <= (lepr_end + cis_window)

cat("LEPR instruments:\n")
cat(sprintf("  Total: %d SNPs\n", nrow(lepr_snps)))
cat(sprintf("  Cis (chr1, within 1Mb of LEPR): %d SNP(s)\n", sum(lepr_snps$is_cis)))
cat(sprintf("  Trans: %d SNPs\n\n", sum(!lepr_snps$is_cis)))

# Read normoglycemic BMI GWAS
cat("Reading normoglycemic BMI GWAS results...\n")
gwas_dir <- "/n/scratch/users/s/st320/regenie/BMI/normoglycemic/step2"

outcome_list <- list()
for (i in 1:nrow(lepr_snps)) {
  snp <- lepr_snps$SNP[i]
  chr <- lepr_snps$chr[i]
  gwas_file <- file.path(gwas_dir, sprintf("regenie_step2_BMI_chr%d_BMI.regenie", chr))
  
  if (file.exists(gwas_file)) {
    gwas <- read.table(gwas_file, header = TRUE, stringsAsFactors = FALSE)
    match_row <- gwas[gwas$ID == snp, ]
    if (nrow(match_row) > 0) {
      outcome_list[[snp]] <- data.frame(
        SNP = snp,
        beta_out = match_row$BETA,
        se_out = match_row$SE,
        pval_out = 10^(-match_row$LOG10P),
        ea_out = match_row$ALLELE1,
        oa_out = match_row$ALLELE0,
        n_out = match_row$N,
        stringsAsFactors = FALSE
      )
    } else {
      cat(sprintf("  Warning: %s not found\n", snp))
    }
  }
}

outcome_df <- do.call(rbind, outcome_list)
cat(sprintf("  Found %d/%d SNPs\n\n", nrow(outcome_df), nrow(lepr_snps)))

# Merge exposure and outcome
dat <- merge(lepr_snps, outcome_df, by = "SNP")
cat(sprintf("Sample size (normoglycemic): N = %d\n\n", dat$n_out[1]))

# Harmonize alleles
for (i in 1:nrow(dat)) {
  if (dat$effect_allele[i] != dat$ea_out[i]) {
    if (dat$effect_allele[i] == dat$oa_out[i]) {
      dat$beta_out[i] <- -dat$beta_out[i]
    }
  }
}

# MR Functions
wald_ratio <- function(beta_exp, se_exp, beta_out, se_out) {
  beta_mr <- beta_out / beta_exp
  se_mr <- abs(se_out / beta_exp)
  pval <- 2 * pnorm(-abs(beta_mr / se_mr))
  list(beta = beta_mr, se = se_mr, pval = pval)
}

ivw_mr <- function(beta_exp, se_exp, beta_out, se_out) {
  weights <- 1 / se_out^2
  beta_mr <- sum(weights * beta_out / beta_exp) / sum(weights / beta_exp^2)
  se_mr <- sqrt(1 / sum(weights / beta_exp^2))
  pval <- 2 * pnorm(-abs(beta_mr / se_mr))
  
  # Cochran Q
  beta_snp <- beta_out / beta_exp
  Q <- sum(weights * (beta_snp - beta_mr)^2)
  Q_df <- length(beta_exp) - 1
  Q_pval <- pchisq(Q, Q_df, lower.tail = FALSE)
  
  list(beta = beta_mr, se = se_mr, pval = pval, Q = Q, Q_df = Q_df, Q_pval = Q_pval)
}

weighted_median_mr <- function(beta_exp, se_exp, beta_out, se_out, n_boot = 1000) {
  beta_snp <- beta_out / beta_exp
  se_snp <- abs(se_out / beta_exp)
  weights <- 1 / se_snp^2
  weights <- weights / sum(weights)
  
  ord <- order(beta_snp)
  beta_snp <- beta_snp[ord]
  weights <- weights[ord]
  cum_weights <- cumsum(weights)
  
  idx <- which(cum_weights >= 0.5)[1]
  beta_mr <- beta_snp[idx]
  
  # Bootstrap SE
  boot_betas <- numeric(n_boot)
  for (b in 1:n_boot) {
    boot_beta_snp <- rnorm(length(beta_snp), beta_snp, se_snp)
    ord_b <- order(boot_beta_snp)
    boot_beta_snp <- boot_beta_snp[ord_b]
    weights_b <- weights[ord_b]
    cum_b <- cumsum(weights_b)
    idx_b <- which(cum_b >= 0.5)[1]
    boot_betas[b] <- boot_beta_snp[idx_b]
  }
  se_mr <- sd(boot_betas)
  pval <- 2 * pnorm(-abs(beta_mr / se_mr))
  
  list(beta = beta_mr, se = se_mr, pval = pval)
}

# ============================================
# 1. ALL pQTLs MR
# ============================================
cat("================================================================\n")
cat("1. ALL pQTLs MR (LEPR -> BMI Normoglycemic)\n")
cat("================================================================\n\n")

ivw_res <- ivw_mr(dat$beta, dat$se, dat$beta_out, dat$se_out)
wm_res <- weighted_median_mr(dat$beta, dat$se, dat$beta_out, dat$se_out)

cat(sprintf("%-28s | %5s | %8s | %8s | %12s | %8s\n",
            "Method", "nSNP", "Beta", "SE", "P-value", "OR"))
cat(paste(rep("-", 85), collapse=""), "\n")

# IVW
pval_fmt <- ifelse(ivw_res$pval < 0.001, sprintf("%.2e", ivw_res$pval), sprintf("%.4f", ivw_res$pval))
sig <- ifelse(ivw_res$pval < 0.05, " *", "")
cat(sprintf("%-28s | %5d | %8.4f | %8.4f | %12s | %8.3f%s\n",
            "Inverse variance weighted", nrow(dat), ivw_res$beta, ivw_res$se, pval_fmt, exp(ivw_res$beta), sig))

# Weighted median
pval_fmt <- ifelse(wm_res$pval < 0.001, sprintf("%.2e", wm_res$pval), sprintf("%.4f", wm_res$pval))
sig <- ifelse(wm_res$pval < 0.05, " *", "")
cat(sprintf("%-28s | %5d | %8.4f | %8.4f | %12s | %8.3f%s\n",
            "Weighted median", nrow(dat), wm_res$beta, wm_res$se, pval_fmt, exp(wm_res$beta), sig))

cat(sprintf("\nHeterogeneity: Q = %.2f, df = %d, P = %.4f\n", 
            ivw_res$Q, ivw_res$Q_df, ivw_res$Q_pval))

# ============================================
# 2. CIS-ONLY MR
# ============================================
cat("\n\n================================================================\n")
cat("2. CIS-ONLY pQTL MR (LEPR -> BMI Normoglycemic)\n")
cat("================================================================\n\n")

dat_cis <- dat[dat$is_cis, ]
cat(sprintf("Cis instruments: %d SNP(s)\n\n", nrow(dat_cis)))

if (nrow(dat_cis) == 1) {
  cis_res <- wald_ratio(dat_cis$beta, dat_cis$se, dat_cis$beta_out, dat_cis$se_out)
  or_cis <- exp(cis_res$beta)
  or_lci <- exp(cis_res$beta - 1.96 * cis_res$se)
  or_uci <- exp(cis_res$beta + 1.96 * cis_res$se)
  
  cat("Method: Wald ratio (single SNP)\n\n")
  cat(sprintf("SNP: %s (chr%d:%d)\n", dat_cis$SNP, dat_cis$chr, dat_cis$pos))
  cat(sprintf("Exposure beta: %.4f (SE: %.4f, P: %.2e)\n", dat_cis$beta, dat_cis$se, dat_cis$pval))
  cat(sprintf("Outcome beta:  %.4f (SE: %.4f)\n\n", dat_cis$beta_out, dat_cis$se_out))
  
  cat(sprintf("%-20s | %5s | %8s | %8s | %12s | %8s | %s\n",
              "Method", "nSNP", "Beta", "SE", "P-value", "OR", "95% CI"))
  cat(paste(rep("-", 90), collapse=""), "\n")
  
  pval_fmt <- ifelse(cis_res$pval < 0.001, sprintf("%.2e", cis_res$pval), sprintf("%.4f", cis_res$pval))
  sig <- ifelse(cis_res$pval < 0.05, " *", "")
  bonf <- ifelse(cis_res$pval < 0.05/349, " **Bonf", "")
  cat(sprintf("%-20s | %5d | %8.4f | %8.4f | %12s | %8.3f | (%.3f, %.3f)%s%s\n",
              "Wald ratio", 1, cis_res$beta, cis_res$se, pval_fmt, or_cis, or_lci, or_uci, sig, bonf))
}

# ============================================
# 3. TRANS-ONLY MR
# ============================================
cat("\n\n================================================================\n")
cat("3. TRANS-ONLY pQTL MR (LEPR -> BMI Normoglycemic)\n")
cat("================================================================\n\n")

dat_trans <- dat[!dat$is_cis, ]
cat(sprintf("Trans instruments: %d SNPs\n\n", nrow(dat_trans)))

ivw_trans <- ivw_mr(dat_trans$beta, dat_trans$se, dat_trans$beta_out, dat_trans$se_out)
wm_trans <- weighted_median_mr(dat_trans$beta, dat_trans$se, dat_trans$beta_out, dat_trans$se_out)

cat(sprintf("%-28s | %5s | %8s | %8s | %12s | %8s\n",
            "Method", "nSNP", "Beta", "SE", "P-value", "OR"))
cat(paste(rep("-", 85), collapse=""), "\n")

pval_fmt <- ifelse(ivw_trans$pval < 0.001, sprintf("%.2e", ivw_trans$pval), sprintf("%.4f", ivw_trans$pval))
sig <- ifelse(ivw_trans$pval < 0.05, " *", "")
cat(sprintf("%-28s | %5d | %8.4f | %8.4f | %12s | %8.3f%s\n",
            "Inverse variance weighted", nrow(dat_trans), ivw_trans$beta, ivw_trans$se, pval_fmt, exp(ivw_trans$beta), sig))

pval_fmt <- ifelse(wm_trans$pval < 0.001, sprintf("%.2e", wm_trans$pval), sprintf("%.4f", wm_trans$pval))
sig <- ifelse(wm_trans$pval < 0.05, " *", "")
cat(sprintf("%-28s | %5d | %8.4f | %8.4f | %12s | %8.3f%s\n",
            "Weighted median", nrow(dat_trans), wm_trans$beta, wm_trans$se, pval_fmt, exp(wm_trans$beta), sig))

# ============================================
# 4. COMPARISON ACROSS STRATA
# ============================================
cat("\n\n================================================================\n")
cat("4. COMPARISON: LEPR CIS-MR -> BMI ACROSS ALL STRATA\n")
cat("================================================================\n\n")

cat(sprintf("%-16s | %10s | %8s | %8s | %12s | %8s | %s\n",
            "Stratum", "N", "Beta", "SE", "P-value", "OR", "95% CI"))
cat(paste(rep("-", 95), collapse=""), "\n")

# Normoglycemic
pval_fmt <- ifelse(cis_res$pval < 0.001, sprintf("%.2e", cis_res$pval), sprintf("%.4f", cis_res$pval))
sig <- ifelse(cis_res$pval < 0.05, " *", "")
cat(sprintf("%-16s | %10d | %8.4f | %8.4f | %12s | %8.3f | (%.3f, %.3f)%s\n",
            "NORMOGLYCEMIC", dat_cis$n_out[1], cis_res$beta, cis_res$se, pval_fmt, 
            or_cis, or_lci, or_uci, sig))

# Other strata
strata_data <- list(
  list("Full cohort", 500000, 0.1171, 0.0549, 0.0330),
  list("Prediabetes", 50000, 0.0991, 0.1832, 0.5886),
  list("T2D", 30000, 0.7034, 0.4160, 0.0908)
)

for (s in strata_data) {
  b <- s[[3]]; se_v <- s[[4]]; p <- s[[5]]
  or_v <- exp(b); or_l <- exp(b - 1.96 * se_v); or_u <- exp(b + 1.96 * se_v)
  pval_fmt <- ifelse(p < 0.001, sprintf("%.2e", p), sprintf("%.4f", p))
  sig <- ifelse(p < 0.05, " *", "")
  cat(sprintf("%-16s | %10d | %8.4f | %8.4f | %12s | %8.3f | (%.3f, %.3f)%s\n",
              s[[1]], s[[2]], b, se_v, pval_fmt, or_v, or_l, or_u, sig))
}

cat("\n================================================================\n")
cat("KEY FINDING\n")
cat("================================================================\n")
cat(sprintf("\nLEPR cis-MR in NORMOGLYCEMIC individuals:\n"))
cat(sprintf("  Beta = %.4f, SE = %.4f, P = %.4f\n", cis_res$beta, cis_res$se, cis_res$pval))
cat(sprintf("  OR = %.3f (95%% CI: %.3f - %.3f)\n", or_cis, or_lci, or_uci))
cat(sprintf("  Sample size: N = %d\n", dat_cis$n_out[1]))
cat("================================================================\n")
