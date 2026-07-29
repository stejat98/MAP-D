#!/usr/bin/env Rscript
# =============================================================================
# Reverse MR: BMI -> All Forward MR Proteins (using DECODE pQTLs as outcome)
#
# Exposure: UKB BMI GWAS (full stratum, held-out non-Olink sample)
# Outcome:  DECODE Iceland SomaScan pQTLs (no sample overlap)
# =============================================================================

.libPaths("/path/to/project/.R/library")
library(data.table)

cat("================================================================================\n")
cat("  Reverse MR: BMI -> All Proteins (DECODE pQTL outcomes)\n")
cat("================================================================================\n\n")

setwd("/path/to/project/regenie_pipeline")
output_dir <- "results/twosampleMR/bidirectional_bmi_proteins"

# =============================================================================
# 1. Load BMI instruments (pre-computed from pilot run)
# =============================================================================
cat("Step 1: Loading UKB BMI GWAS and extracting instruments...\n")

bmi_file <- "results/GWAS/BMI/BMI_full_all_chr.regenie.gz"

cmd <- sprintf("zcat %s | awk 'NR==1 || $12 > 7.30103'", bmi_file)
bmi_sig <- fread(cmd = cmd)
bmi_sig <- bmi_sig[grepl("^rs", ID), ]
bmi_sig <- bmi_sig[order(-LOG10P)]
bmi_sig$pval <- 10^(-as.numeric(bmi_sig$LOG10P))

# Distance-based pruning
keep <- rep(FALSE, nrow(bmi_sig))
keep[1] <- TRUE
for (i in 2:nrow(bmi_sig)) {
  cur_chr <- bmi_sig$CHROM[i]
  cur_pos <- bmi_sig$GENPOS[i]
  same_chr_idx <- which(keep & bmi_sig$CHROM == cur_chr)
  if (length(same_chr_idx) == 0 || 
      all(abs(cur_pos - bmi_sig$GENPOS[same_chr_idx]) > 500000)) {
    keep[i] <- TRUE
  }
}
bmi_instruments <- bmi_sig[keep, ]
cat(sprintf("  Independent BMI instruments: %d\n", nrow(bmi_instruments)))
cat(sprintf("  F-stat: median = %.1f, min = %.1f\n\n", 
    median((as.numeric(bmi_instruments$BETA)/as.numeric(bmi_instruments$SE))^2),
    min((as.numeric(bmi_instruments$BETA)/as.numeric(bmi_instruments$SE))^2)))

snp_list <- bmi_instruments$ID

# Write SNP list to temp file for efficient grep
tmp_snp_file <- tempfile()
writeLines(snp_list, tmp_snp_file)

# =============================================================================
# 2. Load protein-DECODE mapping
# =============================================================================
mapping <- read.csv(file.path(output_dir, "protein_decode_mapping.csv"))
cat(sprintf("Proteins to test: %d\n\n", nrow(mapping)))

decode_dir <- "/path/to/shared_data/DECODE/pQTL/final_somascan_smp"

# =============================================================================
# 3. Run reverse MR for each protein
# =============================================================================
all_results <- data.frame()

for (idx in 1:nrow(mapping)) {
  prot_name <- mapping$protein[idx]
  decode_file <- file.path(decode_dir, mapping$decode_file[idx])
  
  if (idx %% 25 == 1) {
    cat(sprintf("Processing protein %d/%d: %s...\n", idx, nrow(mapping), prot_name))
  }
  
  # Extract matching SNPs from DECODE pQTL
  cmd <- sprintf("zcat %s | head -1 && zcat %s | grep -wFf %s", 
                 decode_file, decode_file, tmp_snp_file)
  
  pqtl <- tryCatch({
    fread(cmd = cmd, showProgress = FALSE)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(pqtl) || nrow(pqtl) < 3) {
    all_results <- rbind(all_results, data.frame(
      Protein = prot_name, Method = "IVW", N_SNPs = ifelse(is.null(pqtl), 0, nrow(pqtl)),
      Beta = NA, SE = NA, P = NA, Q = NA, Q_P = NA,
      Egger_Beta = NA, Egger_SE = NA, Egger_P = NA,
      Egger_Intercept = NA, Egger_Intercept_P = NA,
      stringsAsFactors = FALSE))
    next
  }
  
  # Extract first rsID for matching
  pqtl$SNP <- sapply(strsplit(pqtl$rsids, ","), function(x) x[1])
  
  # Merge
  merged <- merge(
    bmi_instruments[, .(SNP = ID, ea_bmi = ALLELE1, oa_bmi = ALLELE0,
                        beta_bmi = as.numeric(BETA), se_bmi = as.numeric(SE))],
    pqtl[, .(SNP, ea_pqtl = effectAllele, oa_pqtl = otherAllele,
             beta_pqtl = Beta, se_pqtl = SE, pval_pqtl = Pval)],
    by = "SNP"
  )
  
  if (nrow(merged) < 3) {
    all_results <- rbind(all_results, data.frame(
      Protein = prot_name, Method = "IVW", N_SNPs = nrow(merged),
      Beta = NA, SE = NA, P = NA, Q = NA, Q_P = NA,
      Egger_Beta = NA, Egger_SE = NA, Egger_P = NA,
      Egger_Intercept = NA, Egger_Intercept_P = NA,
      stringsAsFactors = FALSE))
    next
  }
  
  # Harmonize alleles
  for (i in 1:nrow(merged)) {
    if (toupper(merged$ea_bmi[i]) == toupper(merged$oa_pqtl[i]) &&
        toupper(merged$oa_bmi[i]) == toupper(merged$ea_pqtl[i])) {
      merged$beta_pqtl[i] <- -merged$beta_pqtl[i]
    } else if (!(toupper(merged$ea_bmi[i]) == toupper(merged$ea_pqtl[i]) &&
                 toupper(merged$oa_bmi[i]) == toupper(merged$oa_pqtl[i]))) {
      merged$beta_pqtl[i] <- NA  # can't harmonize
    }
  }
  merged <- merged[!is.na(merged$beta_pqtl), ]
  
  if (nrow(merged) < 3) {
    all_results <- rbind(all_results, data.frame(
      Protein = prot_name, Method = "IVW", N_SNPs = nrow(merged),
      Beta = NA, SE = NA, P = NA, Q = NA, Q_P = NA,
      Egger_Beta = NA, Egger_SE = NA, Egger_P = NA,
      Egger_Intercept = NA, Egger_Intercept_P = NA,
      stringsAsFactors = FALSE))
    next
  }
  
  # Wald ratios
  merged$wald_beta <- merged$beta_pqtl / merged$beta_bmi
  merged$wald_se <- abs(merged$se_pqtl / merged$beta_bmi)
  
  # IVW
  weights <- 1 / merged$wald_se^2
  ivw_beta <- sum(weights * merged$wald_beta) / sum(weights)
  ivw_se <- sqrt(1 / sum(weights))
  ivw_z <- ivw_beta / ivw_se
  ivw_p <- 2 * pnorm(-abs(ivw_z))
  
  # Cochran's Q
  Q <- sum(weights * (merged$wald_beta - ivw_beta)^2)
  Q_df <- nrow(merged) - 1
  Q_p <- pchisq(Q, Q_df, lower.tail = FALSE)
  
  # MR-Egger
  egger_beta <- egger_se <- egger_p <- egger_int <- egger_int_p <- NA
  if (nrow(merged) >= 3) {
    egger_fit <- tryCatch({
      lm(merged$beta_pqtl ~ merged$beta_bmi, weights = 1/merged$se_pqtl^2)
    }, error = function(e) NULL)
    
    if (!is.null(egger_fit)) {
      ec <- summary(egger_fit)$coefficients
      egger_beta <- ec[2, "Estimate"]
      egger_se <- ec[2, "Std. Error"]
      egger_p <- ec[2, "Pr(>|t|)"]
      egger_int <- ec[1, "Estimate"]
      egger_int_p <- ec[1, "Pr(>|t|)"]
    }
  }
  
  all_results <- rbind(all_results, data.frame(
    Protein = prot_name, Method = "IVW", N_SNPs = nrow(merged),
    Beta = ivw_beta, SE = ivw_se, P = ivw_p, Q = Q, Q_P = Q_p,
    Egger_Beta = egger_beta, Egger_SE = egger_se, Egger_P = egger_p,
    Egger_Intercept = egger_int, Egger_Intercept_P = egger_int_p,
    stringsAsFactors = FALSE))
}

unlink(tmp_snp_file)

# =============================================================================
# 4. Merge with forward MR results
# =============================================================================
cat("\n\nStep 4: Merging with forward MR results...\n")

forward <- read.csv("results/twosampleMR/supplemental_tables/Table_MR_Cis_pQTL_BMI_Full.csv")
forward <- forward[, c("Protein", "N.SNPs", "Beta", "SE", "P.value", "Bonferroni.Significant")]
names(forward) <- c("Protein", "Fwd_N_SNPs", "Fwd_Beta", "Fwd_SE", "Fwd_P", "Fwd_Bonf")

all_results$Rev_Beta <- all_results$Beta
all_results$Rev_SE <- all_results$SE
all_results$Rev_P <- all_results$P

combined <- merge(forward, 
                  all_results[, c("Protein", "N_SNPs", "Rev_Beta", "Rev_SE", "Rev_P",
                                  "Q", "Q_P", "Egger_Beta", "Egger_Intercept", "Egger_Intercept_P")],
                  by = "Protein", all.x = TRUE)

names(combined)[names(combined) == "N_SNPs"] <- "Rev_N_SNPs"

# Multiple testing correction for reverse MR
combined$Rev_P_BH <- p.adjust(combined$Rev_P, method = "BH")
combined$Rev_Bonf <- ifelse(combined$Rev_P < 0.05 / sum(!is.na(combined$Rev_P)), "Yes", "No")

# Direction concordance
combined$Fwd_Dir <- ifelse(as.numeric(combined$Fwd_Beta) > 0, "+", "-")
combined$Rev_Dir <- ifelse(combined$Rev_Beta > 0, "+", "-")
combined$Concordant <- combined$Fwd_Dir == combined$Rev_Dir

# Sort by reverse MR p-value
combined <- combined[order(combined$Rev_P), ]

# =============================================================================
# 5. Save and print results
# =============================================================================
write.csv(combined, file.path(output_dir, "bidirectional_MR_BMI_all_proteins.csv"),
          row.names = FALSE)

cat("\n================================================================================\n")
cat("  TOP 30 REVERSE MR RESULTS: BMI -> Protein (by P-value)\n")
cat("================================================================================\n\n")

top <- head(combined[!is.na(combined$Rev_P), ], 30)
cat(sprintf("%-12s %6s %8s %8s %12s %6s | %6s %8s %12s %6s | %5s\n",
            "Protein", "N_SNP", "Beta", "SE", "P", "Bonf",
            "N_SNP", "Beta", "P", "Bonf", "Conc"))
cat(sprintf("%-12s %6s %8s %8s %12s %6s | %6s %8s %12s %6s | %5s\n",
            "", "(Rev)", "(Rev)", "(Rev)", "(Rev)", "(Rev)",
            "(Fwd)", "(Fwd)", "(Fwd)", "(Fwd)", ""))
cat(paste(rep("-", 110), collapse=""), "\n")
for (i in 1:nrow(top)) {
  cat(sprintf("%-12s %6d %8.4f %8.4f %12s %6s | %6d %8.4f %12s %6s | %5s\n",
              top$Protein[i], top$Rev_N_SNPs[i], top$Rev_Beta[i], top$Rev_SE[i],
              formatC(top$Rev_P[i], format="e", digits=2), top$Rev_Bonf[i],
              top$Fwd_N_SNPs[i], as.numeric(top$Fwd_Beta[i]),
              formatC(as.numeric(top$Fwd_P[i]), format="e", digits=2), top$Fwd_Bonf[i],
              ifelse(top$Concordant[i], "Yes", "NO")))
}

cat("\n\n================================================================================\n")
cat("  SUMMARY STATISTICS\n")
cat("================================================================================\n\n")

n_tested <- sum(!is.na(combined$Rev_P))
n_sig_bonf <- sum(combined$Rev_Bonf == "Yes", na.rm=TRUE)
n_sig_fdr <- sum(combined$Rev_P_BH < 0.05, na.rm=TRUE)
n_sig_nom <- sum(combined$Rev_P < 0.05, na.rm=TRUE)

cat(sprintf("  Proteins tested:                     %d\n", n_tested))
cat(sprintf("  Significant (Bonferroni P < %.2e): %d\n", 0.05/n_tested, n_sig_bonf))
cat(sprintf("  Significant (FDR < 0.05):            %d\n", n_sig_fdr))
cat(sprintf("  Nominally significant (P < 0.05):    %d\n\n", n_sig_nom))

# Direction concordance among proteins significant in both directions
both_sig <- combined[!is.na(combined$Rev_P) & combined$Rev_P < 0.05 & 
                     as.numeric(combined$Fwd_P) < 0.05, ]
if (nrow(both_sig) > 0) {
  cat(sprintf("  Proteins nominally significant in BOTH directions: %d\n", nrow(both_sig)))
  cat(sprintf("  Direction concordant: %d (%.1f%%)\n", 
              sum(both_sig$Concordant), 100*mean(both_sig$Concordant)))
  cat(sprintf("  Direction discordant: %d (%.1f%%)\n\n", 
              sum(!both_sig$Concordant), 100*mean(!both_sig$Concordant)))
  
  cat("  Proteins with DISCORDANT directions (nominally sig both ways):\n")
  disc <- both_sig[!both_sig$Concordant, ]
  if (nrow(disc) > 0) {
    for (i in 1:nrow(disc)) {
      cat(sprintf("    %s: Fwd=%s%.4f (P=%s), Rev=%s%.4f (P=%s)\n",
                  disc$Protein[i], disc$Fwd_Dir[i], abs(as.numeric(disc$Fwd_Beta[i])),
                  formatC(as.numeric(disc$Fwd_P[i]), format="e", digits=2),
                  disc$Rev_Dir[i], abs(disc$Rev_Beta[i]),
                  formatC(disc$Rev_P[i], format="e", digits=2)))
    }
  }
}

# Forward-only significant (suggests protein -> BMI causal)
fwd_only <- combined[!is.na(combined$Rev_P) & 
                     as.numeric(combined$Fwd_P) < 0.05/nrow(forward) & 
                     combined$Rev_P > 0.05, ]
cat(sprintf("\n  Proteins with Bonferroni-sig FORWARD but non-sig REVERSE: %d\n", nrow(fwd_only)))
if (nrow(fwd_only) > 0) {
  cat("  (Suggests protein -> BMI causal direction):\n")
  for (i in 1:nrow(fwd_only)) {
    cat(sprintf("    %s: Fwd Beta=%.4f P=%s\n", fwd_only$Protein[i], 
                as.numeric(fwd_only$Fwd_Beta[i]),
                formatC(as.numeric(fwd_only$Fwd_P[i]), format="e", digits=2)))
  }
}

cat(sprintf("\nResults saved to: %s/bidirectional_MR_BMI_all_proteins.csv\n", output_dir))
cat("\n================================================================================\n")
cat("                          DONE\n")
cat("================================================================================\n")
