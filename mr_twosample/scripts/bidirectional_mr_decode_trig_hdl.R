#!/usr/bin/env Rscript
# =============================================================================
# Bidirectional MR: DECODE cis-pQTLs <-> UKB held-out TRIG_HDL_RATIO GWAS
#
# Forward MR: DECODE cis-pQTL instruments -> UKB TRIG_HDL_RATIO outcome
# Reverse MR: Load existing results (UKB TRIG_HDL_RATIO instruments -> DECODE)
#
# No sample overlap: DECODE (Iceland) vs UKB (UK held-out non-Olink)
# =============================================================================

.libPaths("/path/to/project/.R/library")
library(data.table)
library(readxl)

cat("================================================================================\n")
cat("  Bidirectional MR: DECODE cis-pQTLs <-> UKB TRIG_HDL_RATIO\n")
cat("================================================================================\n\n")

setwd("/path/to/project/regenie_pipeline")
output_dir <- "results/twosampleMR/bidirectional_decode_trig_hdl"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# STEP 1: Define cis regions per protein using UKB-PPP ST10 annotations
#   ST10 has cis/trans labels; we use the rsIDs of cis-pQTLs as "anchors"
#   to locate the gene region in DECODE coordinates (avoids build issues)
# =============================================================================
cat("Step 1: Defining cis regions from UKB-PPP ST10 annotations...\n")

st10 <- read_excel("41586_2023_6592_MOESM3_ESM.xlsx", sheet = "ST10", skip = 3)
colnames(st10) <- gsub("[ ():/]", "_", colnames(st10))
st10$Gene_Symbol <- sapply(strsplit(st10$UKBPPP_ProteinID, ":"), function(x) x[1])

# Get cis-pQTL entries only
cis_entries <- st10[st10$cis_trans == "cis", ]
cat(sprintf("  Total cis-pQTL entries in ST10: %d\n", nrow(cis_entries)))
cat(sprintf("  Proteins with cis-pQTLs: %d\n", length(unique(cis_entries$Gene_Symbol))))

# Store cis-pQTL rsIDs per protein (for anchor-based region definition in DECODE)
cis_rsids_by_protein <- split(cis_entries$rsID, cis_entries$Gene_Symbol)

# =============================================================================
# STEP 2: Load protein-DECODE mapping
# =============================================================================
cat("\nStep 2: Loading protein-DECODE mapping...\n")

mapping_file <- "results/twosampleMR/bidirectional_bmi_proteins/protein_decode_mapping.csv"
mapping <- read.csv(mapping_file)

# Only keep proteins that have cis-pQTLs defined in ST10
mapping <- mapping[mapping$protein %in% names(cis_rsids_by_protein), ]
cat(sprintf("  Proteins with both DECODE data and ST10 cis-pQTLs: %d\n", nrow(mapping)))

decode_dir <- "/path/to/shared_data/DECODE/pQTL/final_somascan_smp"

# =============================================================================
# STEP 3: Load UKB TRIG_HDL_RATIO GWAS (outcome for forward MR)
# =============================================================================
cat("\nStep 3: Loading UKB TRIG_HDL_RATIO GWAS...\n")

gwas_file <- "results/GWAS/TRIG_HDL_RATIO/TRIG_HDL_RATIO_full_all_chr.regenie.gz"
outcome_gwas <- fread(cmd = sprintf("zcat %s", gwas_file))
cat(sprintf("  Total SNPs in outcome GWAS: %d\n", nrow(outcome_gwas)))

# Index by SNP ID for fast lookup
setkey(outcome_gwas, ID)

# =============================================================================
# STEP 4: Forward MR — DECODE cis-pQTL instruments -> UKB TRIG_HDL_RATIO
#
# For each protein:
#   a) Load DECODE full GWAS for that protein
#   b) Find ST10 cis-pQTL rsIDs in DECODE -> get positions -> define cis window
#   c) Extract genome-wide significant SNPs in cis window from DECODE
#   d) Distance-prune (500kb)
#   e) Harmonize with UKB TRIG_HDL_RATIO GWAS
#   f) Run MR (Wald ratio for 1 SNP, IVW for 2+, +Egger for 3+)
# =============================================================================
cat("\nStep 4: Running forward MR (DECODE cis-pQTLs -> TRIG_HDL_RATIO)...\n\n")

fwd_results <- data.frame()
skipped <- list(no_decode = 0, no_anchor = 0, no_cis_sig = 0,
                no_overlap = 0, no_harmonize = 0)

for (idx in 1:nrow(mapping)) {
  prot <- mapping$protein[idx]
  decode_file <- file.path(decode_dir, mapping$decode_file[idx])

  if (idx %% 25 == 1) {
    cat(sprintf("  Processing protein %d/%d: %s...\n", idx, nrow(mapping), prot))
  }

  # a) Load DECODE full GWAS for this protein
  pqtl <- tryCatch(
    fread(cmd = sprintf("zcat %s", decode_file), showProgress = FALSE),
    error = function(e) NULL
  )
  if (is.null(pqtl) || nrow(pqtl) == 0) {
    skipped$no_decode <- skipped$no_decode + 1
    next
  }

  # Parse first rsID from comma-separated list
  pqtl$SNP <- sapply(strsplit(pqtl$rsids, ","), `[`, 1)

  # b) Find anchor SNPs: ST10 cis-pQTL rsIDs present in DECODE
  anchor_rsids <- cis_rsids_by_protein[[prot]]
  anchors <- pqtl[SNP %in% anchor_rsids, ]

  if (nrow(anchors) == 0) {
    skipped$no_anchor <- skipped$no_anchor + 1
    next
  }

  # Define cis region from anchor positions (±1Mb)
  cis_chr <- anchors$Chrom[which.min(anchors$Pval)]
  cis_start <- max(1, min(anchors[Chrom == cis_chr, ]$Pos) - 1e6)
  cis_end <- max(anchors[Chrom == cis_chr, ]$Pos) + 1e6

  # c) Extract DECODE cis-significant SNPs (P < 5e-8) in cis window
  cis_sig <- pqtl[Chrom == cis_chr & Pos >= cis_start & Pos <= cis_end & Pval < 5e-8, ]

  if (nrow(cis_sig) == 0) {
    skipped$no_cis_sig <- skipped$no_cis_sig + 1
    next
  }

  # d) Distance-based pruning (500kb, keep most significant first)
  cis_sig <- cis_sig[order(Pval), ]
  keep <- rep(FALSE, nrow(cis_sig))
  keep[1] <- TRUE
  if (nrow(cis_sig) > 1) {
    for (i in 2:nrow(cis_sig)) {
      kept_pos <- cis_sig$Pos[which(keep)]
      if (all(abs(cis_sig$Pos[i] - kept_pos) > 500000)) {
        keep[i] <- TRUE
      }
    }
  }
  instruments <- cis_sig[keep, ]
  instruments <- instruments[grepl("^rs", SNP), ]

  if (nrow(instruments) == 0) {
    skipped$no_cis_sig <- skipped$no_cis_sig + 1
    next
  }

  # F-statistics for instruments
  fstats <- (instruments$Beta / instruments$SE)^2

  # e) Look up instruments in UKB TRIG_HDL_RATIO GWAS & harmonize
  merged <- merge(
    instruments[, .(SNP, ea_pqtl = effectAllele, oa_pqtl = otherAllele,
                    beta_pqtl = Beta, se_pqtl = SE)],
    outcome_gwas[, .(SNP = ID, ea_gwas = ALLELE1, oa_gwas = ALLELE0,
                     beta_gwas = as.numeric(BETA), se_gwas = as.numeric(SE))],
    by = "SNP"
  )

  if (nrow(merged) == 0) {
    skipped$no_overlap <- skipped$no_overlap + 1
    next
  }

  # Harmonize alleles
  for (i in 1:nrow(merged)) {
    ea_p <- toupper(merged$ea_pqtl[i])
    oa_p <- toupper(merged$oa_pqtl[i])
    ea_g <- toupper(merged$ea_gwas[i])
    oa_g <- toupper(merged$oa_gwas[i])

    if (ea_p == ea_g && oa_p == oa_g) {
      # Already aligned — do nothing
    } else if (ea_p == oa_g && oa_p == ea_g) {
      # Flipped — flip outcome beta
      merged$beta_gwas[i] <- -merged$beta_gwas[i]
    } else {
      # Can't harmonize (strand ambiguity etc.)
      merged$beta_gwas[i] <- NA
    }
  }
  merged <- merged[!is.na(merged$beta_gwas), ]

  if (nrow(merged) == 0) {
    skipped$no_harmonize <- skipped$no_harmonize + 1
    next
  }

  # f) MR estimation
  n_snps <- nrow(merged)
  mean_f <- mean(fstats[1:min(length(fstats), n_snps)])

  if (n_snps == 1) {
    # Wald ratio
    beta_mr <- merged$beta_gwas[1] / merged$beta_pqtl[1]
    se_mr <- abs(merged$se_gwas[1] / merged$beta_pqtl[1])
    p_mr <- 2 * pnorm(-abs(beta_mr / se_mr))

    fwd_results <- rbind(fwd_results, data.frame(
      Protein = prot, Fwd_Method = "Wald ratio", Fwd_N_SNPs = 1,
      Fwd_Beta = beta_mr, Fwd_SE = se_mr, Fwd_P = p_mr,
      Fwd_Q = NA, Fwd_Q_P = NA,
      Fwd_Egger_Beta = NA, Fwd_Egger_P = NA,
      Fwd_Egger_Intercept = NA, Fwd_Egger_Intercept_P = NA,
      Fwd_MeanF = mean_f,
      stringsAsFactors = FALSE))

  } else {
    # IVW (2+ SNPs)
    merged$wald_beta <- merged$beta_gwas / merged$beta_pqtl
    merged$wald_se <- abs(merged$se_gwas / merged$beta_pqtl)
    weights <- 1 / merged$wald_se^2
    ivw_beta <- sum(weights * merged$wald_beta) / sum(weights)
    ivw_se <- sqrt(1 / sum(weights))
    ivw_p <- 2 * pnorm(-abs(ivw_beta / ivw_se))
    Q <- sum(weights * (merged$wald_beta - ivw_beta)^2)
    Q_p <- pchisq(Q, n_snps - 1, lower.tail = FALSE)

    # MR-Egger (3+ SNPs)
    egger_beta <- egger_p <- egger_int <- egger_int_p <- NA
    if (n_snps >= 3) {
      egger_fit <- tryCatch(
        lm(merged$beta_gwas ~ merged$beta_pqtl, weights = 1 / merged$se_gwas^2),
        error = function(e) NULL)
      if (!is.null(egger_fit)) {
        ec <- summary(egger_fit)$coefficients
        egger_beta <- ec[2, "Estimate"]
        egger_p <- ec[2, "Pr(>|t|)"]
        egger_int <- ec[1, "Estimate"]
        egger_int_p <- ec[1, "Pr(>|t|)"]
      }
    }

    fwd_results <- rbind(fwd_results, data.frame(
      Protein = prot, Fwd_Method = "IVW", Fwd_N_SNPs = n_snps,
      Fwd_Beta = ivw_beta, Fwd_SE = ivw_se, Fwd_P = ivw_p,
      Fwd_Q = Q, Fwd_Q_P = Q_p,
      Fwd_Egger_Beta = egger_beta, Fwd_Egger_P = egger_p,
      Fwd_Egger_Intercept = egger_int, Fwd_Egger_Intercept_P = egger_int_p,
      Fwd_MeanF = mean_f,
      stringsAsFactors = FALSE))
  }
}

cat(sprintf("\n  Forward MR completed: %d proteins with results\n", nrow(fwd_results)))
cat(sprintf("  Skipped: no DECODE=%d, no anchor=%d, no cis-sig=%d, no overlap=%d, no harmonize=%d\n",
    skipped$no_decode, skipped$no_anchor, skipped$no_cis_sig,
    skipped$no_overlap, skipped$no_harmonize))

# Instrument summary
cat(sprintf("  Wald ratio (1 SNP): %d proteins\n", sum(fwd_results$Fwd_Method == "Wald ratio")))
cat(sprintf("  IVW (2+ SNPs):      %d proteins\n", sum(fwd_results$Fwd_Method == "IVW")))
cat(sprintf("  Mean F-stat across proteins: %.1f\n",
    mean(fwd_results$Fwd_MeanF, na.rm = TRUE)))

# Multiple testing correction
fwd_results$Fwd_P_BH <- p.adjust(fwd_results$Fwd_P, method = "BH")
n_fwd <- sum(!is.na(fwd_results$Fwd_P))
fwd_results$Fwd_Bonf <- ifelse(!is.na(fwd_results$Fwd_P) &
                                  fwd_results$Fwd_P < 0.05 / n_fwd, "Yes", "No")

# Save forward-only results
write.csv(fwd_results, file.path(output_dir, "forward_MR_DECODE_cis_TRIG_HDL.csv"),
          row.names = FALSE)
cat(sprintf("  Saved: %s\n", file.path(output_dir, "forward_MR_DECODE_cis_TRIG_HDL.csv")))

# =============================================================================
# STEP 5: Load existing reverse MR results
# =============================================================================
cat("\nStep 5: Loading existing reverse MR results...\n")

rev_file <- "results/twosampleMR/bidirectional_trig_hdl_ratio_proteins/bidirectional_MR_TRIG_HDL_RATIO_all_proteins.csv"
rev_raw <- read.csv(rev_file)

# Rename to Rev_ prefix for clarity
rev_results <- data.frame(
  Protein = rev_raw$Protein,
  Rev_N_SNPs = rev_raw$N_SNPs,
  Rev_Beta = rev_raw$IVW_Beta,
  Rev_SE = rev_raw$IVW_SE,
  Rev_P = rev_raw$IVW_P,
  Rev_Q = rev_raw$Q,
  Rev_Q_P = rev_raw$Q_P,
  Rev_Egger_Beta = rev_raw$Egger_Beta,
  Rev_Egger_P = rev_raw$Egger_P,
  Rev_Egger_Intercept = rev_raw$Egger_Intercept,
  Rev_Egger_Intercept_P = rev_raw$Egger_Intercept_P,
  Rev_P_BH = rev_raw$IVW_P_BH,
  Rev_Bonf = rev_raw$IVW_Bonf,
  stringsAsFactors = FALSE
)

cat(sprintf("  Reverse MR results loaded: %d proteins\n", nrow(rev_results)))

# =============================================================================
# STEP 6: Merge forward + reverse
# =============================================================================
cat("\nStep 6: Merging bidirectional results...\n")

combined <- merge(fwd_results, rev_results, by = "Protein", all = TRUE)

# Direction concordance
combined$Fwd_Dir <- ifelse(combined$Fwd_Beta > 0, "+", "-")
combined$Rev_Dir <- ifelse(combined$Rev_Beta > 0, "+", "-")
combined$Concordant <- ifelse(!is.na(combined$Fwd_Dir) & !is.na(combined$Rev_Dir),
                               combined$Fwd_Dir == combined$Rev_Dir, NA)

# Sort by forward MR p-value
combined <- combined[order(combined$Fwd_P), ]

# Save combined results
write.csv(combined, file.path(output_dir, "bidirectional_MR_DECODE_TRIG_HDL.csv"),
          row.names = FALSE)
cat(sprintf("  Saved: %s\n", file.path(output_dir, "bidirectional_MR_DECODE_TRIG_HDL.csv")))

# =============================================================================
# STEP 7: Print results
# =============================================================================
cat("\n================================================================================\n")
cat("  TOP 30 FORWARD MR: Protein (DECODE cis-pQTL) -> TRIG_HDL_RATIO\n")
cat("================================================================================\n\n")

top_fwd <- head(combined[!is.na(combined$Fwd_P), ], 30)

cat(sprintf("%-12s %5s %5s %8s %12s %5s | %5s %8s %12s %5s | %5s\n",
            "Protein", "Meth", "N_SNP", "FwdBeta", "FwdP", "FBonf",
            "N_SNP", "RevBeta", "RevP", "RBonf", "Conc"))
cat(paste(rep("-", 100), collapse = ""), "\n")

for (i in 1:nrow(top_fwd)) {
  meth <- ifelse(top_fwd$Fwd_Method[i] == "Wald ratio", "Wald", "IVW")
  rev_n <- ifelse(is.na(top_fwd$Rev_N_SNPs[i]), 0, top_fwd$Rev_N_SNPs[i])
  rev_b <- ifelse(is.na(top_fwd$Rev_Beta[i]), "    NA", sprintf("%8.4f", top_fwd$Rev_Beta[i]))
  rev_p <- ifelse(is.na(top_fwd$Rev_P[i]), "          NA", formatC(top_fwd$Rev_P[i], format = "e", digits = 2))
  rev_bonf <- ifelse(is.na(top_fwd$Rev_Bonf[i]), "  NA", sprintf("%5s", top_fwd$Rev_Bonf[i]))
  conc <- ifelse(is.na(top_fwd$Concordant[i]), "  NA",
                 ifelse(top_fwd$Concordant[i], "  Yes", "   NO"))

  cat(sprintf("%-12s %5s %5d %8.4f %12s %5s | %5d %s %s %s | %s\n",
              top_fwd$Protein[i], meth, top_fwd$Fwd_N_SNPs[i],
              top_fwd$Fwd_Beta[i],
              formatC(top_fwd$Fwd_P[i], format = "e", digits = 2),
              top_fwd$Fwd_Bonf[i],
              rev_n, rev_b, rev_p, rev_bonf, conc))
}

# =============================================================================
# STEP 8: Summary statistics
# =============================================================================
cat("\n\n================================================================================\n")
cat("  SUMMARY\n")
cat("================================================================================\n\n")

cat("  --- Forward MR: Protein (DECODE cis-pQTL) -> TRIG_HDL_RATIO ---\n")
cat(sprintf("  Proteins tested:              %d\n", n_fwd))
cat(sprintf("  Bonferroni significant:       %d (P < %.2e)\n",
    sum(fwd_results$Fwd_Bonf == "Yes", na.rm = TRUE), 0.05 / n_fwd))
cat(sprintf("  FDR < 0.05:                   %d\n",
    sum(fwd_results$Fwd_P_BH < 0.05, na.rm = TRUE)))
cat(sprintf("  Nominally significant:        %d\n\n",
    sum(fwd_results$Fwd_P < 0.05, na.rm = TRUE)))

n_rev <- sum(!is.na(rev_results$Rev_P))
cat("  --- Reverse MR: TRIG_HDL_RATIO -> Protein (DECODE outcome) ---\n")
cat(sprintf("  Proteins tested:              %d\n", n_rev))
cat(sprintf("  Bonferroni significant:       %d\n",
    sum(rev_results$Rev_Bonf == "Yes", na.rm = TRUE)))
cat(sprintf("  FDR < 0.05:                   %d\n",
    sum(rev_results$Rev_P_BH < 0.05, na.rm = TRUE)))
cat(sprintf("  Nominally significant:        %d\n\n",
    sum(rev_results$Rev_P < 0.05, na.rm = TRUE)))

# Bidirectional summary
both_tested <- combined[!is.na(combined$Fwd_P) & !is.na(combined$Rev_P), ]
both_sig_nom <- both_tested[both_tested$Fwd_P < 0.05 & both_tested$Rev_P < 0.05, ]
fwd_only_bonf <- both_tested[!is.na(both_tested$Fwd_Bonf) & both_tested$Fwd_Bonf == "Yes" &
                               (is.na(both_tested$Rev_P) | both_tested$Rev_P > 0.05), ]
rev_only_bonf <- both_tested[!is.na(both_tested$Rev_Bonf) & both_tested$Rev_Bonf == "Yes" &
                               (is.na(both_tested$Fwd_P) | both_tested$Fwd_P > 0.05), ]

cat("  --- Bidirectional ---\n")
cat(sprintf("  Proteins tested both directions:   %d\n", nrow(both_tested)))

if (nrow(both_sig_nom) > 0) {
  cat(sprintf("  Nominally sig both directions:     %d\n", nrow(both_sig_nom)))
  n_conc <- sum(both_sig_nom$Concordant, na.rm = TRUE)
  n_disc <- sum(!both_sig_nom$Concordant, na.rm = TRUE)
  cat(sprintf("    Concordant: %d | Discordant: %d\n", n_conc, n_disc))

  if (n_disc > 0) {
    disc <- both_sig_nom[!is.na(both_sig_nom$Concordant) & !both_sig_nom$Concordant, ]
    cat("\n  Discordant proteins (nom. sig both ways, opposite directions):\n")
    for (i in 1:nrow(disc)) {
      cat(sprintf("    %s: Fwd=%s%.4f (P=%s), Rev=%s%.4f (P=%s)\n",
                  disc$Protein[i], disc$Fwd_Dir[i], abs(disc$Fwd_Beta[i]),
                  formatC(disc$Fwd_P[i], format = "e", digits = 2),
                  disc$Rev_Dir[i], abs(disc$Rev_Beta[i]),
                  formatC(disc$Rev_P[i], format = "e", digits = 2)))
    }
  }
}

cat(sprintf("\n  Bonf-sig FORWARD only (Protein -> TRIG_HDL): %d\n", nrow(fwd_only_bonf)))
if (nrow(fwd_only_bonf) > 0) {
  cat("    (Evidence for protein -> TRIG_HDL causal direction):\n")
  for (i in 1:min(nrow(fwd_only_bonf), 20)) {
    cat(sprintf("      %s: Fwd Beta=%.4f P=%s (F=%.1f, %d SNPs)\n",
                fwd_only_bonf$Protein[i], fwd_only_bonf$Fwd_Beta[i],
                formatC(fwd_only_bonf$Fwd_P[i], format = "e", digits = 2),
                fwd_only_bonf$Fwd_MeanF[i], fwd_only_bonf$Fwd_N_SNPs[i]))
  }
}

cat(sprintf("\n  Bonf-sig REVERSE only (TRIG_HDL -> Protein): %d\n", nrow(rev_only_bonf)))
if (nrow(rev_only_bonf) > 0) {
  cat("    (Evidence for TRIG_HDL -> protein causal direction):\n")
  for (i in 1:min(nrow(rev_only_bonf), 20)) {
    cat(sprintf("      %s: Rev Beta=%.4f P=%s\n",
                rev_only_bonf$Protein[i], rev_only_bonf$Rev_Beta[i],
                formatC(rev_only_bonf$Rev_P[i], format = "e", digits = 2)))
  }
}

cat(sprintf("\nResults: %s/bidirectional_MR_DECODE_TRIG_HDL.csv\n", output_dir))
cat("\n================================================================================\n")
cat("                          DONE\n")
cat("================================================================================\n")
