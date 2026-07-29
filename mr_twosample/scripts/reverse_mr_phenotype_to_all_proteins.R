#!/usr/bin/env Rscript
# =============================================================================
# Reverse MR: Phenotype -> All Proteins (using DECODE pQTLs as outcome)
#
# Usage: Rscript reverse_mr_phenotype_to_all_proteins.R <PHENOTYPE>
#   e.g.: Rscript reverse_mr_phenotype_to_all_proteins.R HbA1c
#         Rscript reverse_mr_phenotype_to_all_proteins.R TRIG_HDL_RATIO
#
# Exposure: UKB hallmark GWAS (full stratum, held-out non-Olink sample)
# Outcome:  DECODE Iceland SomaScan pQTLs (no sample overlap)
# =============================================================================

.libPaths("/path/to/project/.R/library")
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript reverse_mr_phenotype_to_all_proteins.R <PHENOTYPE>")
PHENO <- args[1]

cat("================================================================================\n")
cat(sprintf("  Reverse MR: %s -> All Proteins (DECODE pQTL outcomes)\n", PHENO))
cat("================================================================================\n\n")

setwd("/path/to/project/regenie_pipeline")
output_dir <- sprintf("results/twosampleMR/bidirectional_%s_proteins", tolower(PHENO))
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. Load phenotype GWAS and extract instruments
# =============================================================================
cat("Step 1: Loading GWAS and extracting instruments...\n")

gwas_file <- sprintf("results/GWAS/%s/%s_full_all_chr.regenie.gz", PHENO, PHENO)
if (!file.exists(gwas_file)) stop(sprintf("GWAS file not found: %s", gwas_file))

cmd <- sprintf("zcat %s | awk 'NR==1 || $12 > 7.30103'", gwas_file)  # P < 5e-8
pheno_sig <- fread(cmd = cmd)
cat(sprintf("  GW-significant %s SNPs (P < 5e-8): %d\n", PHENO, nrow(pheno_sig)))

pheno_sig <- pheno_sig[grepl("^rs", ID), ]
cat(sprintf("  With rsIDs: %d\n", nrow(pheno_sig)))

pheno_sig <- pheno_sig[order(-LOG10P)]
pheno_sig$pval <- 10^(-as.numeric(pheno_sig$LOG10P))

# Distance-based pruning: keep top SNP per 500kb window
keep <- rep(FALSE, nrow(pheno_sig))
if (nrow(pheno_sig) > 0) {
  keep[1] <- TRUE
  for (i in 2:nrow(pheno_sig)) {
    cur_chr <- pheno_sig$CHROM[i]
    cur_pos <- pheno_sig$GENPOS[i]
    same_chr_idx <- which(keep & pheno_sig$CHROM == cur_chr)
    if (length(same_chr_idx) == 0 || 
        all(abs(cur_pos - pheno_sig$GENPOS[same_chr_idx]) > 500000)) {
      keep[i] <- TRUE
    }
  }
}
pheno_instruments <- pheno_sig[keep, ]
cat(sprintf("  Independent instruments after 500kb pruning: %d\n", nrow(pheno_instruments)))

fstats <- (as.numeric(pheno_instruments$BETA)/as.numeric(pheno_instruments$SE))^2
cat(sprintf("  F-stat: median = %.1f, min = %.1f, max = %.1f\n", 
    median(fstats), min(fstats), max(fstats)))
cat(sprintf("  All F > 10: %s\n\n", ifelse(all(fstats > 10), "YES", "NO")))

snp_list <- pheno_instruments$ID
tmp_snp_file <- tempfile()
writeLines(snp_list, tmp_snp_file)

# =============================================================================
# 2. Load protein-DECODE mapping
# =============================================================================
# Use same mapping as BMI (same set of forward MR proteins)
mapping_file <- "results/twosampleMR/bidirectional_bmi_proteins/protein_decode_mapping.csv"
if (!file.exists(mapping_file)) {
  # Generate mapping on the fly
  cat("Generating protein-DECODE mapping...\n")
  decode_dir <- "/path/to/shared_data/DECODE/pQTL/final_somascan_smp"
  decode_files <- list.files(decode_dir)
  
  get_gene <- function(f) {
    parts <- strsplit(f, "_")[[1]]
    if (length(parts) >= 7) return(parts[6])
    return(NA)
  }
  decode_genes <- sapply(decode_files, get_gene)
  
  # Try to load forward MR table for this phenotype; fall back to BMI
  fwd_file <- sprintf("results/twosampleMR/supplemental_tables/Table_MR_Cis_pQTL_%s_Full.csv", PHENO)
  if (file.exists(fwd_file)) {
    mr <- read.csv(fwd_file)
  } else {
    mr <- read.csv("results/twosampleMR/supplemental_tables/Table_MR_Cis_pQTL_BMI_Full.csv")
    cat(sprintf("  No forward MR table for %s; using BMI protein list\n", PHENO))
  }
  
  mapping <- data.frame(protein = mr$Protein, decode_file = NA, stringsAsFactors = FALSE)
  for (i in 1:nrow(mapping)) {
    idx <- which(decode_genes == mapping$protein[i])
    if (length(idx) > 0) mapping$decode_file[i] <- names(decode_genes)[idx[1]]
  }
  mapping <- mapping[!is.na(mapping$decode_file), ]
  write.csv(mapping, file.path(output_dir, "protein_decode_mapping.csv"), row.names = FALSE)
} else {
  mapping <- read.csv(mapping_file)
}

# If forward MR exists for this phenotype, also get its protein list
fwd_file <- sprintf("results/twosampleMR/supplemental_tables/Table_MR_Cis_pQTL_%s_Full.csv", PHENO)
has_forward <- file.exists(fwd_file)
if (has_forward) {
  fwd_mr <- read.csv(fwd_file)
  # Restrict to proteins that were in the forward MR for this phenotype
  mapping <- mapping[mapping$protein %in% fwd_mr$Protein, ]
}

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
  
  cmd <- sprintf("zcat %s | head -1 && zcat %s | grep -wFf %s", 
                 decode_file, decode_file, tmp_snp_file)
  
  pqtl <- tryCatch(fread(cmd = cmd, showProgress = FALSE), error = function(e) NULL)
  
  if (is.null(pqtl) || nrow(pqtl) < 3) {
    all_results <- rbind(all_results, data.frame(
      Protein = prot_name, N_SNPs = ifelse(is.null(pqtl), 0, nrow(pqtl)),
      IVW_Beta = NA, IVW_SE = NA, IVW_P = NA, Q = NA, Q_P = NA,
      Egger_Beta = NA, Egger_SE = NA, Egger_P = NA,
      Egger_Intercept = NA, Egger_Intercept_P = NA,
      stringsAsFactors = FALSE))
    next
  }
  
  pqtl$SNP <- sapply(strsplit(pqtl$rsids, ","), function(x) x[1])
  
  merged <- merge(
    pheno_instruments[, .(SNP = ID, ea_exp = ALLELE1, oa_exp = ALLELE0,
                          beta_exp = as.numeric(BETA), se_exp = as.numeric(SE))],
    pqtl[, .(SNP, ea_pqtl = effectAllele, oa_pqtl = otherAllele,
             beta_pqtl = Beta, se_pqtl = SE, pval_pqtl = Pval)],
    by = "SNP"
  )
  
  if (nrow(merged) < 3) {
    all_results <- rbind(all_results, data.frame(
      Protein = prot_name, N_SNPs = nrow(merged),
      IVW_Beta = NA, IVW_SE = NA, IVW_P = NA, Q = NA, Q_P = NA,
      Egger_Beta = NA, Egger_SE = NA, Egger_P = NA,
      Egger_Intercept = NA, Egger_Intercept_P = NA,
      stringsAsFactors = FALSE))
    next
  }
  
  # Harmonize alleles
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
  
  if (nrow(merged) < 3) {
    all_results <- rbind(all_results, data.frame(
      Protein = prot_name, N_SNPs = nrow(merged),
      IVW_Beta = NA, IVW_SE = NA, IVW_P = NA, Q = NA, Q_P = NA,
      Egger_Beta = NA, Egger_SE = NA, Egger_P = NA,
      Egger_Intercept = NA, Egger_Intercept_P = NA,
      stringsAsFactors = FALSE))
    next
  }
  
  # Wald ratios & IVW
  merged$wald_beta <- merged$beta_pqtl / merged$beta_exp
  merged$wald_se <- abs(merged$se_pqtl / merged$beta_exp)
  
  weights <- 1 / merged$wald_se^2
  ivw_beta <- sum(weights * merged$wald_beta) / sum(weights)
  ivw_se <- sqrt(1 / sum(weights))
  ivw_p <- 2 * pnorm(-abs(ivw_beta / ivw_se))
  
  Q <- sum(weights * (merged$wald_beta - ivw_beta)^2)
  Q_p <- pchisq(Q, nrow(merged) - 1, lower.tail = FALSE)
  
  # MR-Egger
  egger_beta <- egger_se <- egger_p <- egger_int <- egger_int_p <- NA
  if (nrow(merged) >= 3) {
    egger_fit <- tryCatch(
      lm(merged$beta_pqtl ~ merged$beta_exp, weights = 1/merged$se_pqtl^2),
      error = function(e) NULL)
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
    Protein = prot_name, N_SNPs = nrow(merged),
    IVW_Beta = ivw_beta, IVW_SE = ivw_se, IVW_P = ivw_p, Q = Q, Q_P = Q_p,
    Egger_Beta = egger_beta, Egger_SE = egger_se, Egger_P = egger_p,
    Egger_Intercept = egger_int, Egger_Intercept_P = egger_int_p,
    stringsAsFactors = FALSE))
}

unlink(tmp_snp_file)

# =============================================================================
# 4. Merge with forward MR if available
# =============================================================================
if (has_forward) {
  cat("\nMerging with forward MR results...\n")
  forward <- fwd_mr[, c("Protein", "N.SNPs", "Beta", "SE", "P.value", "Bonferroni.Significant")]
  names(forward) <- c("Protein", "Fwd_N_SNPs", "Fwd_Beta", "Fwd_SE", "Fwd_P", "Fwd_Bonf")
  
  combined <- merge(forward, all_results, by = "Protein", all.x = TRUE)
} else {
  cat(sprintf("\nNo forward MR table for %s; reporting reverse-only.\n", PHENO))
  combined <- all_results
  combined$Fwd_Beta <- NA
  combined$Fwd_P <- NA
  combined$Fwd_Bonf <- NA
}

# Multiple testing correction
combined$IVW_P_BH <- p.adjust(combined$IVW_P, method = "BH")
n_tested <- sum(!is.na(combined$IVW_P))
combined$IVW_Bonf <- ifelse(!is.na(combined$IVW_P) & combined$IVW_P < 0.05 / n_tested, "Yes", "No")

# Direction concordance (if forward available)
if (has_forward) {
  combined$Fwd_Dir <- ifelse(as.numeric(combined$Fwd_Beta) > 0, "+", "-")
  combined$Rev_Dir <- ifelse(combined$IVW_Beta > 0, "+", "-")
  combined$Concordant <- combined$Fwd_Dir == combined$Rev_Dir
}

combined <- combined[order(combined$IVW_P), ]

# Save
write.csv(combined, file.path(output_dir, sprintf("bidirectional_MR_%s_all_proteins.csv", PHENO)),
          row.names = FALSE)

# =============================================================================
# 5. Print results
# =============================================================================
cat(sprintf("\n================================================================================\n"))
cat(sprintf("  TOP 30 REVERSE MR: %s -> Protein (by P-value)\n", PHENO))
cat(sprintf("================================================================================\n\n"))

top <- head(combined[!is.na(combined$IVW_P), ], 30)
if (has_forward) {
  cat(sprintf("%-12s %6s %8s %8s %12s %6s | %8s %12s %6s | %5s\n",
              "Protein", "N_SNP", "RevBeta", "RevSE", "RevP", "Bonf",
              "FwdBeta", "FwdP", "FBonf", "Conc"))
  cat(paste(rep("-", 105), collapse=""), "\n")
  for (i in 1:nrow(top)) {
    cat(sprintf("%-12s %6d %8.4f %8.4f %12s %6s | %8.4f %12s %6s | %5s\n",
                top$Protein[i], top$N_SNPs[i], top$IVW_Beta[i], top$IVW_SE[i],
                formatC(top$IVW_P[i], format="e", digits=2), top$IVW_Bonf[i],
                as.numeric(top$Fwd_Beta[i]),
                formatC(as.numeric(top$Fwd_P[i]), format="e", digits=2), top$Fwd_Bonf[i],
                ifelse(is.na(top$Concordant[i]), "NA", ifelse(top$Concordant[i], "Yes", "NO"))))
  }
} else {
  cat(sprintf("%-12s %6s %8s %8s %12s %6s\n",
              "Protein", "N_SNP", "RevBeta", "RevSE", "RevP", "Bonf"))
  cat(paste(rep("-", 60), collapse=""), "\n")
  for (i in 1:nrow(top)) {
    cat(sprintf("%-12s %6d %8.4f %8.4f %12s %6s\n",
                top$Protein[i], top$N_SNPs[i], top$IVW_Beta[i], top$IVW_SE[i],
                formatC(top$IVW_P[i], format="e", digits=2), top$IVW_Bonf[i]))
  }
}

cat(sprintf("\n\n================================================================================\n"))
cat(sprintf("  SUMMARY: Reverse MR %s -> Protein\n", PHENO))
cat(sprintf("================================================================================\n\n"))

cat(sprintf("  Instruments: %d independent %s SNPs (P < 5e-8, 500kb pruned)\n", 
    nrow(pheno_instruments), PHENO))
cat(sprintf("  Proteins tested:                     %d\n", n_tested))
cat(sprintf("  Significant (Bonferroni P < %.2e): %d\n", 0.05/n_tested,
    sum(combined$IVW_Bonf == "Yes", na.rm=TRUE)))
cat(sprintf("  Significant (FDR < 0.05):            %d\n", 
    sum(combined$IVW_P_BH < 0.05, na.rm=TRUE)))
cat(sprintf("  Nominally significant (P < 0.05):    %d\n\n", 
    sum(combined$IVW_P < 0.05, na.rm=TRUE)))

if (has_forward) {
  both_sig <- combined[!is.na(combined$IVW_P) & combined$IVW_P < 0.05 & 
                       as.numeric(combined$Fwd_P) < 0.05, ]
  if (nrow(both_sig) > 0) {
    cat(sprintf("  Nom. sig both directions: %d | Concordant: %d | Discordant: %d\n",
                nrow(both_sig), sum(both_sig$Concordant, na.rm=TRUE), 
                sum(!both_sig$Concordant, na.rm=TRUE)))
  }
}

cat(sprintf("\nResults: %s/bidirectional_MR_%s_all_proteins.csv\n", output_dir, PHENO))
cat("\n================================================================================\n")
cat("                          DONE\n")
cat("================================================================================\n")
