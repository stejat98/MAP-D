#!/usr/bin/env Rscript
# =============================================================================
# Merge Forward (UKB cis-pQTL -> TRIG_HDL_RATIO) + Reverse (TRIG_HDL_RATIO -> protein)
# into a single bidirectional file matching the BMI/HbA1c format
# =============================================================================

setwd("/n/groups/patel/sivateja/regenie_pipeline")

cat("================================================================================\n")
cat("  Merging bidirectional MR for TRIG_HDL_RATIO\n")
cat("  Forward: UKB cis-pQTLs -> TRIG_HDL_RATIO (from cis_trans_chunked)\n")
cat("  Reverse: TRIG_HDL_RATIO -> Protein (from reverse MR)\n")
cat("================================================================================\n\n")

# ---- 1. Load forward cis-pQTL results (protein -> TRIG_HDL_RATIO) ----
cat("Step 1: Loading forward cis-pQTL MR results...\n")

cis_dir <- "results/twosampleMR/cis_trans_chunked"
cis_files <- list.files(cis_dir, pattern = "MR_cis_chunk.*\\.csv", full.names = TRUE)
cis_all <- do.call(rbind, lapply(cis_files, read.csv))

fwd <- cis_all[cis_all$phenotype == "TRIG_HDL_RATIO" & cis_all$stratum == "full", ]
cat(sprintf("  Forward results: %d rows for %d proteins\n", nrow(fwd), length(unique(fwd$protein))))

n_fwd_proteins <- length(unique(fwd$protein))
bonf_fwd <- 0.05 / n_fwd_proteins

fwd_summary <- data.frame(
  Protein = fwd$protein,
  Fwd_N_SNPs = fwd$nsnp,
  Fwd_Beta = fwd$b,
  Fwd_SE = fwd$se,
  Fwd_P = fwd$pval,
  Fwd_Bonf = ifelse(fwd$pval < bonf_fwd, "Yes", "No"),
  stringsAsFactors = FALSE
)

cat(sprintf("  Bonferroni threshold (forward): %.2e\n", bonf_fwd))
cat(sprintf("  Forward Bonferroni significant: %d\n", sum(fwd_summary$Fwd_Bonf == "Yes")))

# ---- 2. Load reverse results (TRIG_HDL_RATIO -> protein) ----
cat("\nStep 2: Loading reverse MR results...\n")

rev_file <- "results/twosampleMR/bidirectional_trig_hdl_ratio_proteins/bidirectional_MR_TRIG_HDL_RATIO_all_proteins.csv"
rev <- read.csv(rev_file)
# Drop old Fwd_ columns from reverse file (they were all NA)
rev <- rev[, !grepl("^Fwd_", names(rev))]
cat(sprintf("  Reverse results: %d proteins\n", nrow(rev)))

# ---- 3. Merge ----
cat("\nStep 3: Merging forward and reverse...\n")

merged <- merge(fwd_summary, rev, by = "Protein", all = TRUE)

# Add direction columns
merged$Fwd_Dir <- ifelse(is.na(merged$Fwd_Beta), NA,
                         ifelse(merged$Fwd_Beta > 0, "+", "-"))
merged$Rev_Dir <- ifelse(is.na(merged$IVW_Beta), NA,
                         ifelse(merged$IVW_Beta > 0, "+", "-"))
merged$Concordant <- ifelse(is.na(merged$Fwd_Dir) | is.na(merged$Rev_Dir), NA,
                            merged$Fwd_Dir == merged$Rev_Dir)

# ---- 4. Reformat to match BMI bidirectional structure ----
out <- data.frame(
  Protein = merged$Protein,
  Fwd_N_SNPs = merged$Fwd_N_SNPs,
  Fwd_Beta = merged$Fwd_Beta,
  Fwd_SE = merged$Fwd_SE,
  Fwd_P = merged$Fwd_P,
  Fwd_Bonf = merged$Fwd_Bonf,
  Rev_N_SNPs = merged$N_SNPs,
  Rev_Beta = merged$IVW_Beta,
  Rev_SE = merged$IVW_SE,
  Rev_P = merged$IVW_P,
  Q = merged$Q,
  Q_P = merged$Q_P,
  Egger_Beta = merged$Egger_Beta,
  Egger_Intercept = merged$Egger_Intercept,
  Egger_Intercept_P = merged$Egger_Intercept_P,
  Rev_P_BH = merged$IVW_P_BH,
  Rev_Bonf = merged$IVW_Bonf,
  Fwd_Dir = merged$Fwd_Dir,
  Rev_Dir = merged$Rev_Dir,
  Concordant = merged$Concordant,
  stringsAsFactors = FALSE
)

# Sort by forward p-value
out <- out[order(out$Fwd_P), ]

# ---- 5. Save ----
outfile <- "results/twosampleMR/bidirectional_trig_hdl_ratio_proteins/bidirectional_MR_TRIG_HDL_RATIO_all_proteins.csv"
write.csv(out, outfile, row.names = FALSE)

cat(sprintf("\nSaved: %s\n", outfile))
cat(sprintf("  Total proteins: %d\n", nrow(out)))
cat(sprintf("  Forward Bonf sig: %d\n", sum(out$Fwd_Bonf == "Yes", na.rm = TRUE)))
cat(sprintf("  Reverse Bonf sig: %d\n", sum(out$Rev_Bonf == "Yes", na.rm = TRUE)))
cat(sprintf("  Both Bonf sig: %d\n", sum(out$Fwd_Bonf == "Yes" & out$Rev_Bonf == "Yes", na.rm = TRUE)))

both_sig <- out[!is.na(out$Fwd_Bonf) & out$Fwd_Bonf == "Yes" &
                !is.na(out$Rev_Bonf) & out$Rev_Bonf == "Yes", ]
if (nrow(both_sig) > 0) {
  cat("\n  Both directions Bonferroni significant:\n")
  for (i in 1:nrow(both_sig)) {
    cat(sprintf("    %s: Fwd=%.4f (P=%.2e), Rev=%.4f (P=%.2e) | %s\n",
                both_sig$Protein[i],
                both_sig$Fwd_Beta[i], both_sig$Fwd_P[i],
                both_sig$Rev_Beta[i], both_sig$Rev_P[i],
                ifelse(both_sig$Concordant[i], "CONCORDANT", "DISCORDANT")))
  }
}

fwd_only <- out[!is.na(out$Fwd_Bonf) & out$Fwd_Bonf == "Yes" &
                (is.na(out$Rev_Bonf) | out$Rev_Bonf != "Yes"), ]
if (nrow(fwd_only) > 0) {
  cat(sprintf("\n  Forward-only Bonf sig (protein -> TRIG_HDL causal): %d\n", nrow(fwd_only)))
  for (i in 1:nrow(fwd_only)) {
    cat(sprintf("    %s: Fwd Beta=%.4f P=%.2e\n",
                fwd_only$Protein[i], fwd_only$Fwd_Beta[i], fwd_only$Fwd_P[i]))
  }
}

cat("\n================================================================================\n")
cat("                          DONE\n")
cat("================================================================================\n")
