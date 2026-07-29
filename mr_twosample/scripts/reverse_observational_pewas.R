#!/usr/bin/env Rscript
# =============================================================================
# Reverse Observational PEWAS: Protein ~ I(scale(Phenotype)) + Covariates
# 
# Generates reverse beta estimates (Phenotype -> Protein direction) for
# integration with bidirectional MR results.
#
# Forward (existing): Phenotype ~ I(scale(Protein)) + Covariates
# Reverse (this script): Protein ~ I(scale(Phenotype)) + Covariates
#
# Note: P-values are identical in both directions (partial correlation is
# symmetric). This script is run to obtain reverse beta estimates which are
# more interpretable when compared to reverse MR (Phenotype -> Protein).
# =============================================================================

source("/path/to/home/Baseline_PEWAS_Linear_Functions_script.R")

cat("================================================================================\n")
cat("  Reverse PEWAS: Phenotype -> Protein Observational Regressions\n")
cat("================================================================================\n\n")

# ---- Load data ----
cat("Loading data...\n")
data <- readRDS("/path/to/project/UKB/PEWAS_results/data_plus_GLP_complications_glycemic_status_HbA1c_adjusted.RDS")
cat(sprintf("Data dimensions: %d x %d\n\n", nrow(data), ncol(data)))

# ---- Adjustments ----
load("/path/to/home/UKB_PEWAS/adjustments_survival_analysis.Rdata")
adjustments <- c(adjustments, "f.74.0.0")
cat(sprintf("Adjustments (%d vars): %s ...\n\n", length(adjustments),
            paste(head(adjustments, 5), collapse = ", ")))

# ---- Protein variables ----
protein_vars <- colnames(data)[grep("protein_x.", colnames(data))][-1]
cat(sprintf("Number of proteins: %d\n\n", length(protein_vars)))

# ---- Stratify by glycemic status ----
data_non_T2D    <- data[!is.na(data$GlycemicStatus) & data$GlycemicStatus == "Normoglycemic", ]
data_prediabetes <- data[!is.na(data$GlycemicStatus) & data$GlycemicStatus == "Prediabetes", ]
data_T2D        <- data[!is.na(data$GlycemicStatus) & data$GlycemicStatus == "Diabetes", ]

cat(sprintf("Stratum sizes: Normoglycemic=%d, Prediabetes=%d, Diabetes=%d\n\n",
            nrow(data_non_T2D), nrow(data_prediabetes), nrow(data_T2D)))

rm(data); gc()

# ---- Phenotypes ----
phenotypes <- c("BMI", "HbA1c", "TRIG_HDL_RATIO")

# ---- Output directory ----
outdir <- "/path/to/project/UKB/PEWAS_results"

# =============================================================================
# REVERSE EWAS FUNCTION
# =============================================================================
run_reverse_ewas <- function(data_stratum, phenotype, adjustments, protein_vars, stratum_name) {

  cat(sprintf("\n--- Reverse EWAS: %s -> Proteins [%s] (%d individuals) ---\n",
              phenotype, stratum_name, nrow(data_stratum)))

  adj_str <- paste(adjustments, collapse = " + ")
  results_list <- vector("list", length(protein_vars))

  base_cols <- c(phenotype, adjustments)
  base_cols <- base_cols[base_cols %in% colnames(data_stratum)]
  base_complete <- complete.cases(data_stratum[, base_cols])

  n_success <- 0
  n_fail <- 0
  t0 <- Sys.time()

  for (i in seq_along(protein_vars)) {
    protein <- protein_vars[i]
    if (i %% 500 == 0) {
      elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
      cat(sprintf("  Protein %d/%d (%.1f min elapsed)...\n", i, length(protein_vars), elapsed))
    }

    tryCatch({
      protein_complete <- base_complete & !is.na(data_stratum[[protein]])
      df <- data_stratum[protein_complete, c(protein, phenotype, adjustments), drop = FALSE]

      if (nrow(df) < 50) {
        n_fail <- n_fail + 1
        next
      }

      formula_str <- sprintf("`%s` ~ I(scale(%s)) + %s", protein, phenotype, adj_str)

      mod <- biglm(as.formula(formula_str), df)
      res <- broom::tidy(mod)

      pheno_pattern <- sprintf("I\\(scale\\(%s\\)\\)", phenotype)
      pheno_row <- res[grep(pheno_pattern, res$term), ]

      if (nrow(pheno_row) > 0) {
        pheno_row$Protein   <- protein
        pheno_row$Phenotype <- phenotype
        pheno_row$Stratum   <- stratum_name
        pheno_row$SampleSize <- nrow(df)
        results_list[[i]] <- pheno_row
        n_success <- n_success + 1
      }
    }, error = function(e) {
      n_fail <<- n_fail + 1
    })
  }

  results_df <- dplyr::bind_rows(results_list)

  # Add FDR correction across all proteins within this phenotype-stratum
  if (nrow(results_df) > 0) {
    results_df$FDR <- p.adjust(results_df$p.value, method = "BH")
  }

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("  Done: %d success, %d failed (%.1f min)\n", n_success, n_fail, elapsed))

  return(results_df)
}

# =============================================================================
# RUN ALL COMBINATIONS
# =============================================================================

cat("\n================================================================================\n")
cat("  Running Reverse Regressions\n")
cat("================================================================================\n")

strata_list <- list(
  non_T2D     = data_non_T2D,
  prediabetes = data_prediabetes,
  T2D         = data_T2D
)

all_results <- list()
total_t0 <- Sys.time()

for (pheno in phenotypes) {
  for (stratum_name in names(strata_list)) {

    results <- run_reverse_ewas(
      strata_list[[stratum_name]], pheno, adjustments, protein_vars, stratum_name
    )

    # Save individual file (matching forward naming convention)
    outfile <- file.path(outdir,
      sprintf("Reverse_%s_Linear_regression_proteomic_lm_results_%s_adj_fasting_time.RDS",
              pheno, stratum_name))
    saveRDS(results, outfile)
    cat(sprintf("  Saved: %s (%d proteins)\n", outfile, nrow(results)))

    all_results[[paste(pheno, stratum_name, sep = "_")]] <- results
  }
}

# =============================================================================
# COMBINED SUMMARY
# =============================================================================

cat("\n================================================================================\n")
cat("  Summary\n")
cat("================================================================================\n\n")

combined <- dplyr::bind_rows(all_results)
combined_file <- file.path(outdir, "Reverse_PEWAS_all_phenotypes_all_strata.RDS")
saveRDS(combined, combined_file)
cat(sprintf("Combined file: %s (%d total rows)\n", combined_file, nrow(combined)))

# Also save as CSV for easy viewing
combined_csv <- file.path(outdir, "Reverse_PEWAS_all_phenotypes_all_strata.csv")
write.csv(combined, combined_csv, row.names = FALSE)
cat(sprintf("Combined CSV: %s\n", combined_csv))

# Print summary statistics
cat("\nResults per phenotype-stratum:\n")
for (key in names(all_results)) {
  r <- all_results[[key]]
  n_nom <- sum(r$p.value < 0.05, na.rm = TRUE)
  n_bonf <- sum(r$p.value < 0.05 / nrow(r), na.rm = TRUE)
  n_fdr <- sum(r$FDR < 0.05, na.rm = TRUE)
  cat(sprintf("  %-30s: %4d proteins (Nominal: %d, Bonf: %d, FDR: %d)\n",
              key, nrow(r), n_nom, n_bonf, n_fdr))
}

total_elapsed <- as.numeric(difftime(Sys.time(), total_t0, units = "mins"))
cat(sprintf("\nTotal runtime: %.1f minutes\n", total_elapsed))

cat("\n================================================================================\n")
cat("  DONE\n")
cat("================================================================================\n")
