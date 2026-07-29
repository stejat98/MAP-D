#!/usr/bin/env Rscript
# =============================================================================
# SENSITIVITY ANALYSIS: Reverse MR (Phenotype -> Protein) with r2 < 0.01 LD clumping
#
# Identical to scripts/reverse_mr_phenotype_to_all_proteins.R EXCEPT the
# instrument-independence step: instead of 500 kb distance pruning, instruments
# are LD-clumped at r2 < 0.01 (10 Mb window) using TwoSampleMR/ieugwasr local
# clumping against the 1000G EUR reference panel.
#
# Purpose: test whether reverse-MR Bonferroni hits are robust to stricter LD
# filtering (collaborator request: "LD less than 0.01 or similar").
#
# Usage: Rscript reverse_mr_r2clump_sensitivity.R <PHENOTYPE> <CHUNK_ID> <N_CHUNKS>
#   e.g.: Rscript reverse_mr_r2clump_sensitivity.R BMI 1 16
#
# Writes per-chunk raw stats to results/twosampleMR/sensitivity_r2clump/<PHENO>/
# (does NOT touch the existing bidirectional_* results). Bonferroni/FDR are
# computed in the separate combine step over the full protein set.
# =============================================================================

.libPaths("/n/groups/patel/sivateja/.R/library")
suppressMessages(library(data.table))
suppressMessages(library(ieugwasr))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript reverse_mr_r2clump_sensitivity.R <PHENOTYPE> <CHUNK_ID> <N_CHUNKS>")
PHENO    <- args[1]
CHUNK_ID <- as.integer(args[2])
N_CHUNKS <- as.integer(args[3])

CLUMP_R2 <- 0.01
CLUMP_KB <- 10000
BFILE    <- "/n/groups/patel/IGLOO/LDref/EUR"
PLINK    <- Sys.getenv("PLINK_BIN", "/n/app/plink/1.90b7.7_20241022/bin/plink")

setwd("/n/groups/patel/sivateja/regenie_pipeline")
output_dir <- sprintf("results/twosampleMR/sensitivity_r2clump/%s", PHENO)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("================================================================================\n")
cat(sprintf("  SENSITIVITY Reverse MR (r2<%.3f clump): %s -> Proteins  [chunk %d/%d]\n",
            CLUMP_R2, PHENO, CHUNK_ID, N_CHUNKS))
cat("================================================================================\n\n")

# =============================================================================
# 1. Load phenotype GWAS and extract instruments (P<5e-8, rsID) then r2<0.01 clump
# =============================================================================
cat("Step 1: Loading GWAS and extracting instruments...\n")
gwas_file <- sprintf("results/GWAS/%s/%s_full_all_chr.regenie.gz", PHENO, PHENO)
if (!file.exists(gwas_file)) stop(sprintf("GWAS file not found: %s", gwas_file))

cmd <- sprintf("zcat %s | awk 'NR==1 || $12 > 7.30103'", gwas_file)  # P < 5e-8
pheno_sig <- fread(cmd = cmd)
pheno_sig <- pheno_sig[grepl("^rs", ID), ]
pheno_sig <- pheno_sig[order(-LOG10P)]
pheno_sig$pval <- 10^(-as.numeric(pheno_sig$LOG10P))
cat(sprintf("  GW-significant %s SNPs with rsIDs (P<5e-8): %d\n", PHENO, nrow(pheno_sig)))

# ---- LD clumping at r2 < 0.01 (replaces 500kb distance pruning) -------------
clump_in <- data.frame(rsid = pheno_sig$ID, pval = pheno_sig$pval,
                       id = PHENO, stringsAsFactors = FALSE)
clumped <- ld_clump(clump_in, clump_r2 = CLUMP_R2, clump_kb = CLUMP_KB,
                    plink_bin = PLINK, bfile = BFILE)
pheno_instruments <- pheno_sig[pheno_sig$ID %in% clumped$rsid, ]
cat(sprintf("  Independent instruments after r2<%.3f clump (%dkb, EUR): %d\n",
            CLUMP_R2, CLUMP_KB, nrow(pheno_instruments)))

fstats <- (as.numeric(pheno_instruments$BETA)/as.numeric(pheno_instruments$SE))^2
cat(sprintf("  F-stat: median = %.1f, min = %.1f, max = %.1f  (all F>10: %s)\n\n",
    median(fstats), min(fstats), max(fstats), ifelse(all(fstats > 10), "YES", "NO")))

snp_list <- pheno_instruments$ID
tmp_snp_file <- tempfile()
writeLines(snp_list, tmp_snp_file)

# =============================================================================
# 2. Protein-DECODE mapping (reuse existing mapping; restrict to forward MR proteins)
# =============================================================================
mapping_file <- "results/twosampleMR/bidirectional_bmi_proteins/protein_decode_mapping.csv"
mapping <- read.csv(mapping_file)

fwd_file <- sprintf("results/twosampleMR/supplemental_tables/Table_MR_Cis_pQTL_%s_Full.csv", PHENO)
has_forward <- file.exists(fwd_file)
if (has_forward) {
  fwd_mr <- read.csv(fwd_file)
  mapping <- mapping[mapping$protein %in% fwd_mr$Protein, ]
}
mapping <- mapping[order(mapping$protein), ]

# ---- chunk the protein list ----
n_prot <- nrow(mapping)
idx_all <- seq_len(n_prot)
chunk_assign <- ((idx_all - 1) %% N_CHUNKS) + 1
my_idx <- idx_all[chunk_assign == CHUNK_ID]
cat(sprintf("Total proteins: %d | this chunk (%d/%d): %d proteins\n\n",
            n_prot, CHUNK_ID, N_CHUNKS, length(my_idx)))

decode_dir <- "/n/groups/patel/IGLOO/DECODE/pQTL/final_somascan_smp"

# =============================================================================
# 3. Reverse MR per protein (identical estimation to the original script)
# =============================================================================
na_row <- function(prot, n) data.frame(
  Protein = prot, N_SNPs = n,
  IVW_Beta = NA, IVW_SE = NA, IVW_P = NA, Q = NA, Q_P = NA,
  Egger_Beta = NA, Egger_SE = NA, Egger_P = NA,
  Egger_Intercept = NA, Egger_Intercept_P = NA, stringsAsFactors = FALSE)

all_results <- data.frame()
for (k in seq_along(my_idx)) {
  idx <- my_idx[k]
  prot_name <- mapping$protein[idx]
  decode_file <- file.path(decode_dir, mapping$decode_file[idx])
  cat(sprintf("[%d/%d] %s\n", k, length(my_idx), prot_name))

  cmd <- sprintf("zcat %s | head -1 && zcat %s | grep -wFf %s",
                 decode_file, decode_file, tmp_snp_file)
  pqtl <- tryCatch(fread(cmd = cmd, showProgress = FALSE), error = function(e) NULL)

  if (is.null(pqtl) || nrow(pqtl) < 3) {
    all_results <- rbind(all_results, na_row(prot_name, ifelse(is.null(pqtl), 0, nrow(pqtl)))); next
  }

  pqtl$SNP <- sapply(strsplit(pqtl$rsids, ","), function(x) x[1])
  merged <- merge(
    pheno_instruments[, .(SNP = ID, ea_exp = ALLELE1, oa_exp = ALLELE0,
                          beta_exp = as.numeric(BETA), se_exp = as.numeric(SE))],
    pqtl[, .(SNP, ea_pqtl = effectAllele, oa_pqtl = otherAllele,
             beta_pqtl = Beta, se_pqtl = SE, pval_pqtl = Pval)],
    by = "SNP")

  if (nrow(merged) < 3) { all_results <- rbind(all_results, na_row(prot_name, nrow(merged))); next }

  for (i in 1:nrow(merged)) {
    if (toupper(merged$ea_exp[i]) == toupper(merged$oa_pqtl[i]) &&
        toupper(merged$oa_exp[i]) == toupper(merged$ea_pqtl[i])) {
      merged$beta_pqtl[i] <- -merged$beta_pqtl[i]
    } else if (!(toupper(merged$ea_exp[i]) == toupper(merged$ea_pqtl[i]) &&
                 toupper(merged$oa_exp[i]) == toupper(merged$oa_pqtl[i]))) {
      merged$beta_pqtl[i] <- NA
    }
  }
  merged <- merged[!is.na(merged$beta_pqtl), ]
  if (nrow(merged) < 3) { all_results <- rbind(all_results, na_row(prot_name, nrow(merged))); next }

  merged$wald_beta <- merged$beta_pqtl / merged$beta_exp
  merged$wald_se <- abs(merged$se_pqtl / merged$beta_exp)
  weights <- 1 / merged$wald_se^2
  ivw_beta <- sum(weights * merged$wald_beta) / sum(weights)
  ivw_se <- sqrt(1 / sum(weights))
  ivw_p <- 2 * pnorm(-abs(ivw_beta / ivw_se))
  Q <- sum(weights * (merged$wald_beta - ivw_beta)^2)
  Q_p <- pchisq(Q, nrow(merged) - 1, lower.tail = FALSE)

  egger_beta <- egger_se <- egger_p <- egger_int <- egger_int_p <- NA
  egger_fit <- tryCatch(lm(merged$beta_pqtl ~ merged$beta_exp, weights = 1/merged$se_pqtl^2),
                        error = function(e) NULL)
  if (!is.null(egger_fit)) {
    ec <- summary(egger_fit)$coefficients
    egger_beta <- ec[2, "Estimate"]; egger_se <- ec[2, "Std. Error"]; egger_p <- ec[2, "Pr(>|t|)"]
    egger_int <- ec[1, "Estimate"]; egger_int_p <- ec[1, "Pr(>|t|)"]
  }

  all_results <- rbind(all_results, data.frame(
    Protein = prot_name, N_SNPs = nrow(merged),
    IVW_Beta = ivw_beta, IVW_SE = ivw_se, IVW_P = ivw_p, Q = Q, Q_P = Q_p,
    Egger_Beta = egger_beta, Egger_SE = egger_se, Egger_P = egger_p,
    Egger_Intercept = egger_int, Egger_Intercept_P = egger_int_p, stringsAsFactors = FALSE))
}
unlink(tmp_snp_file)

# record instrument count used (constant across chunks) for the combine step
out_file <- file.path(output_dir, sprintf("chunk_%03d_of_%03d.csv", CHUNK_ID, N_CHUNKS))
all_results$n_instruments <- nrow(pheno_instruments)
write.csv(all_results, out_file, row.names = FALSE)
cat(sprintf("\nWrote %d protein results to %s\n", nrow(all_results), out_file))
cat("DONE chunk\n")
