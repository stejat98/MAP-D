#!/usr/bin/env Rscript
# Combine sensitivity (r2<0.01) reverse-MR chunks and compare hit robustness vs
# the original 500kb-distance-pruned results.
#
# Usage: Rscript combine_reverse_mr_r2clump.R

.libPaths("/path/to/project/.R/library")
suppressMessages(library(data.table))
setwd("/path/to/project/regenie_pipeline")

orig_files <- c(
  BMI            = "results/twosampleMR/bidirectional_bmi_proteins/bidirectional_MR_BMI_all_proteins.csv",
  TRIG_HDL_RATIO = "results/twosampleMR/bidirectional_trig_hdl_ratio_proteins/bidirectional_MR_TRIG_HDL_RATIO_all_proteins.csv"
)

summarize <- function(PHENO) {
  cdir <- sprintf("results/twosampleMR/sensitivity_r2clump/%s", PHENO)
  chunks <- list.files(cdir, pattern = "^chunk_.*\\.csv$", full.names = TRUE)
  if (length(chunks) == 0) { cat(sprintf("[%s] no chunk files found in %s\n", PHENO, cdir)); return(invisible()) }
  res <- rbindlist(lapply(chunks, fread), fill = TRUE)
  n_instr <- unique(res$n_instruments[!is.na(res$n_instruments)])

  n_tested <- sum(!is.na(res$IVW_P))
  res$IVW_P_BH <- p.adjust(res$IVW_P, method = "BH")
  bonf_thr <- 0.05 / n_tested
  res$IVW_Bonf <- ifelse(!is.na(res$IVW_P) & res$IVW_P < bonf_thr, "Yes", "No")
  res <- res[order(res$IVW_P), ]

  out <- sprintf("%s/COMBINED_%s_r2clump.csv", cdir, PHENO)
  fwrite(res, out)

  new_hits <- res$Protein[res$IVW_Bonf == "Yes"]

  # original hits
  orig <- fread(orig_files[[PHENO]])
  obonf_col <- if ("Rev_Bonf" %in% names(orig)) "Rev_Bonf" else "IVW_Bonf"
  oprot_col <- "Protein"
  orig_hits <- orig[[oprot_col]][orig[[obonf_col]] == "Yes"]
  orig_hits <- orig_hits[!is.na(orig_hits)]

  retained <- intersect(orig_hits, new_hits)
  lost     <- setdiff(orig_hits, new_hits)
  gained   <- setdiff(new_hits, orig_hits)

  cat(sprintf("\n================= %s =================\n", PHENO))
  cat(sprintf("  Instruments: original=(see manuscript)  r2<0.01-clumped=%s\n", paste(n_instr, collapse=",")))
  cat(sprintf("  Proteins tested (non-NA IVW): %d   Bonferroni thr: %.2e\n", n_tested, bonf_thr))
  cat(sprintf("  Original Bonferroni hits:     %d\n", length(orig_hits)))
  cat(sprintf("  Sensitivity Bonferroni hits:  %d\n", length(new_hits)))
  cat(sprintf("  RETAINED (robust):            %d\n", length(retained)))
  cat(sprintf("  LOST (no longer sig):         %d  %s\n", length(lost),
              if (length(lost)) paste0("[", paste(head(lost, 30), collapse=", "), "]") else ""))
  cat(sprintf("  GAINED (newly sig):           %d  %s\n", length(gained),
              if (length(gained)) paste0("[", paste(head(gained, 30), collapse=", "), "]") else ""))
  cat(sprintf("  -> robustness: %.0f%% of original hits retained\n",
              100 * length(retained) / max(1, length(orig_hits))))
  cat(sprintf("  Combined results: %s\n", out))

  data.frame(PHENO = PHENO, orig = length(orig_hits), new = length(new_hits),
             retained = length(retained), lost = length(lost), gained = length(gained))
}

cat("================================================================================\n")
cat("  Reverse-MR LD sensitivity (r2<0.01) vs original (500kb distance prune)\n")
cat("================================================================================\n")
tab <- rbindlist(lapply(c("BMI", "TRIG_HDL_RATIO"), summarize), fill = TRUE)
cat("\n\n==== OVERALL ====\n")
print(tab)
fwrite(tab, "results/twosampleMR/sensitivity_r2clump/SENSITIVITY_SUMMARY.csv")
cat("\nSaved: results/twosampleMR/sensitivity_r2clump/SENSITIVITY_SUMMARY.csv\n")
