# Generate REGENIE input files for Normoglycemic BMI GWAS (held-out sample)

.libPaths("/n/groups/patel/sivateja/.R/library")

library(data.table)
library(dplyr)

cat("=== Generating Normoglycemic BMI GWAS Inputs ===\n\n")

setwd("/n/groups/patel/sivateja/regenie_pipeline")

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------
cat("Step 1: Loading data...\n")

main_df <- as.data.table(readRDS("/n/groups/patel/sivateja/UKB/PEWAS_results/data_plus_GLP_complications_glycemic_status_HbA1c_adjusted.RDS"))
cat(sprintf("  Main data: %d individuals\n", nrow(main_df)))

load("/n/groups/patel/sivateja/olink_eids_for_proteins_gwas.RDS")
olink_eids <- olink_eids_for_proteins_gwas
cat(sprintf("  Proteomics EIDs: %d\n", length(olink_eids)))

# -----------------------------------------------------------------------------
# 2. Filter for normoglycemic held-out sample
# -----------------------------------------------------------------------------
cat("\nStep 2: Filtering for normoglycemic held-out sample...\n")

df_normo_heldout <- main_df[
  GlycemicStatus == "Normoglycemic" & 
  !(eid %in% olink_eids)
]

cat(sprintf("  Normoglycemic held-out: %d individuals\n", nrow(df_normo_heldout)))
cat(sprintf("  HbA1c range: %.1f - %.1f (all < 39 confirmed)\n", 
    min(df_normo_heldout$HbA1c, na.rm=TRUE),
    max(df_normo_heldout$HbA1c, na.rm=TRUE)))

# -----------------------------------------------------------------------------
# 3. Define covariates (same as other strata)
# -----------------------------------------------------------------------------
cat("\nStep 3: Preparing covariates...\n")

covar_cols <- c(
  "eid",
  "x.age", 
  "x.sex",
  "f.54.0.0",
  paste0("f.22009.0.", 1:40),
  "x.738"
)

available_covars <- covar_cols[covar_cols %in% names(df_normo_heldout)]
cat(sprintf("  Using %d covariates\n", length(available_covars) - 1))

# -----------------------------------------------------------------------------
# 4. Create output directory
# -----------------------------------------------------------------------------
output_dir <- "inputs/normoglycemic/hallmarks_heldout/BMI"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("\nOutput directory: %s\n", output_dir))

# -----------------------------------------------------------------------------
# 5. Filter for complete cases
# -----------------------------------------------------------------------------
cat("\nStep 4: Filtering for complete cases...\n")

required_cols <- c("eid", "BMI", available_covars[-1])
df_complete <- df_normo_heldout[complete.cases(df_normo_heldout[, ..required_cols])]
cat(sprintf("  Complete cases (BMI + all covariates): %d\n", nrow(df_complete)))

# -----------------------------------------------------------------------------
# 6. Create phenotype file
# -----------------------------------------------------------------------------
cat("\nStep 5: Creating phenotype file...\n")

pheno_df <- data.frame(
  FID = df_complete$eid,
  IID = df_complete$eid,
  BMI = df_complete$BMI
)

pheno_file <- file.path(output_dir, "pheno.txt")
fwrite(pheno_df, pheno_file, sep = "\t", na = "NA", quote = FALSE)
cat(sprintf("  Saved: %s (%d individuals)\n", pheno_file, nrow(pheno_df)))

# -----------------------------------------------------------------------------
# 7. Create covariate file
# -----------------------------------------------------------------------------
cat("\nStep 6: Creating covariate file...\n")

covar_df <- data.frame(FID = df_complete$eid, IID = df_complete$eid)
for (col in available_covars[-1]) {
  covar_df[[col]] <- df_complete[[col]]
}

covar_file <- file.path(output_dir, "covar.txt")
fwrite(covar_df, covar_file, sep = "\t", na = "NA", quote = FALSE)
cat(sprintf("  Saved: %s (%d covariates)\n", covar_file, ncol(covar_df) - 2))

# -----------------------------------------------------------------------------
# 8. Create pheno list file for per-chromosome jobs
# -----------------------------------------------------------------------------
cat("\nStep 7: Creating pheno list file for per-chromosome jobs...\n")

pheno_list_file <- "scripts/pheno_list_normoglycemic_per_chr.txt"
pheno_lines <- paste0("normoglycemic|hallmarks_heldout|BMI|", 1:22)
writeLines(pheno_lines, pheno_list_file)
cat(sprintf("  Saved: %s (22 chromosome jobs)\n", pheno_list_file))

# -----------------------------------------------------------------------------
# 9. Summary
# -----------------------------------------------------------------------------
cat("\n=== Summary ===\n")
cat(sprintf("Stratum: Normoglycemic (held-out)\n"))
cat(sprintf("Phenotype: BMI\n"))
cat(sprintf("Sample size: %d\n", nrow(df_complete)))
cat(sprintf("HbA1c range: %.1f - %.1f mmol/mol (all < 39)\n",
    min(df_complete$HbA1c, na.rm=TRUE),
    max(df_complete$HbA1c, na.rm=TRUE)))
cat(sprintf("BMI range: %.1f - %.1f\n",
    min(df_complete$BMI), max(df_complete$BMI)))
cat(sprintf("Jobs to submit: 22 (one per chromosome)\n"))

cat("\n=== Done ===\n")
