#!/usr/bin/env Rscript
# ============================================================================
# LEPR ~ BMI: Observational Estimate Adjusted for LEP (Leptin)
#
# The strong negative observational LEPR-BMI association (beta ~ -0.28) is
# likely driven by reverse causation: BMI → LEP↑ → LEPR↑
# Adjusting for LEP removes this pathway and should yield an estimate
# closer to the cis-MR causal estimate (+0.12).
# ============================================================================

.libPaths(c("/n/groups/patel/sivateja/.R/library", "/n/groups/patel/sivateja/R_libs", .libPaths()))
suppressPackageStartupMessages(library(data.table))

cat("================================================================\n")
cat("LEPR ~ BMI: Unadjusted vs Adjusted for LEP\n")
cat("================================================================\n\n")

# Load data
cat("Loading main dataframe...\n")
main_df <- readRDS("/n/groups/patel/sivateja/UKB/PEWAS_results/data_plus_GLP_complications_glycemic_status_HbA1c_adjusted.RDS")
cat(sprintf("  Loaded: %d rows, %d columns\n", nrow(main_df), ncol(main_df)))

# Check columns for LEPR, LEP, BMI
lepr_col <- grep("^protein_x\\.1573$", colnames(main_df), value = TRUE)
lep_cols <- grep("^protein_x\\.1572$", colnames(main_df), value = TRUE)
bmi_col <- "BMI"

cat(sprintf("  LEPR column: %s\n", ifelse(length(lepr_col) > 0, lepr_col, "NOT FOUND")))
cat(sprintf("  LEP column: %s\n", ifelse(length(lep_cols) > 0, paste(lep_cols, collapse = ", "), "NOT FOUND")))
cat(sprintf("  BMI column: %s\n", ifelse(bmi_col %in% colnames(main_df), "found", "NOT FOUND")))

if (length(lepr_col) == 0 || length(lep_cols) == 0 || !(bmi_col %in% colnames(main_df))) {
  cat("\nSearching for alternative column names...\n")
  cat("Columns matching 'LEP':\n")
  print(grep("LEP|lep|leptin", colnames(main_df), value = TRUE, ignore.case = TRUE))
  cat("Columns matching 'BMI':\n")
  print(grep("BMI|bmi", colnames(main_df), value = TRUE, ignore.case = TRUE))
  cat("Columns matching 'protein_x.157':\n")
  print(grep("protein_x\\.157", colnames(main_df), value = TRUE))
}

# Also check for glycemic status for stratification
glyc_col <- grep("glycemic_status|glyc|diabetes_status", colnames(main_df), value = TRUE, ignore.case = TRUE)
cat(sprintf("  Glycemic status column(s): %s\n", paste(glyc_col, collapse = ", ")))

# Standard covariates
age_col <- grep("^x\\.age$|^age$", colnames(main_df), value = TRUE)
sex_col <- grep("^x\\.sex$|^sex$", colnames(main_df), value = TRUE)
cat(sprintf("  Age column: %s\n", paste(age_col, collapse = ", ")))
cat(sprintf("  Sex column: %s\n", paste(sex_col, collapse = ", ")))

# ============================================================================
# Build analysis dataset
# ============================================================================
cat("\n--- Building analysis dataset ---\n")

# Use first matching columns
lepr_var <- lepr_col[1]
lep_var <- lep_cols[1]

# Build dataset with non-missing values
cols_needed <- c("eid", bmi_col, lepr_var, lep_var, age_col[1], sex_col[1], glyc_col[1])
cols_needed <- cols_needed[cols_needed %in% colnames(main_df)]

df <- as.data.table(main_df[, cols_needed])
colnames(df) <- c("eid", "BMI", "LEPR", "LEP", "age", "sex", "glycemic_status")[1:length(cols_needed)]

# Complete cases for LEPR and BMI
df_full <- df[complete.cases(df[, c("BMI", "LEPR"), with = FALSE]), ]
cat(sprintf("  N with BMI + LEPR: %d\n", nrow(df_full)))

df_both <- df[complete.cases(df[, c("BMI", "LEPR", "LEP"), with = FALSE]), ]
cat(sprintf("  N with BMI + LEPR + LEP: %d\n", nrow(df_both)))

# Standardize all variables for comparable betas
standardize <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
df_both$BMI_z <- standardize(df_both$BMI)
df_both$LEPR_z <- standardize(df_both$LEPR)
df_both$LEP_z <- standardize(df_both$LEP)

# ============================================================================
# 1. UNADJUSTED: BMI ~ LEPR (full sample with both)
# ============================================================================
cat("\n================================================================\n")
cat("1. UNADJUSTED: BMI_z ~ LEPR_z\n")
cat("================================================================\n")

m1 <- lm(BMI_z ~ LEPR_z + age + sex, data = df_both)
cat(sprintf("   N = %d\n", nobs(m1)))
cat(sprintf("   LEPR beta = %.4f (SE = %.4f, P = %.2e)\n",
            coef(m1)["LEPR_z"], summary(m1)$coefficients["LEPR_z", "Std. Error"],
            summary(m1)$coefficients["LEPR_z", "Pr(>|t|)"]))

# ============================================================================
# 2. ADJUSTED FOR LEP: BMI ~ LEPR + LEP
# ============================================================================
cat("\n================================================================\n")
cat("2. ADJUSTED FOR LEP: BMI_z ~ LEPR_z + LEP_z\n")
cat("================================================================\n")

m2 <- lm(BMI_z ~ LEPR_z + LEP_z + age + sex, data = df_both)
cat(sprintf("   N = %d\n", nobs(m2)))
cat(sprintf("   LEPR beta = %.4f (SE = %.4f, P = %.2e)\n",
            coef(m2)["LEPR_z"], summary(m2)$coefficients["LEPR_z", "Std. Error"],
            summary(m2)$coefficients["LEPR_z", "Pr(>|t|)"]))
cat(sprintf("   LEP  beta = %.4f (SE = %.4f, P = %.2e)\n",
            coef(m2)["LEP_z"], summary(m2)$coefficients["LEP_z", "Std. Error"],
            summary(m2)$coefficients["LEP_z", "Pr(>|t|)"]))

# ============================================================================
# 3. RESIDUALIZED: Regress out LEP from LEPR, then BMI ~ LEPR_resid
# ============================================================================
cat("\n================================================================\n")
cat("3. RESIDUALIZED: LEPR_resid (= LEPR | LEP), then BMI_z ~ LEPR_resid\n")
cat("================================================================\n")

lep_model <- lm(LEPR_z ~ LEP_z, data = df_both)
df_both$LEPR_resid <- residuals(lep_model)
df_both$LEPR_resid_z <- standardize(df_both$LEPR_resid)

m3 <- lm(BMI_z ~ LEPR_resid_z + age + sex, data = df_both)
cat(sprintf("   N = %d\n", nobs(m3)))
cat(sprintf("   LEPR_resid beta = %.4f (SE = %.4f, P = %.2e)\n",
            coef(m3)["LEPR_resid_z"], summary(m3)$coefficients["LEPR_resid_z", "Std. Error"],
            summary(m3)$coefficients["LEPR_resid_z", "Pr(>|t|)"]))

r2_lep <- summary(lep_model)$r.squared
cat(sprintf("\n   LEP explains %.1f%% of LEPR variance\n", r2_lep * 100))

# ============================================================================
# 4. STRATIFIED ANALYSIS (by glycemic status)
# ============================================================================
if ("glycemic_status" %in% colnames(df_both)) {
  cat("\n================================================================\n")
  cat("4. STRATIFIED BY GLYCEMIC STATUS\n")
  cat("================================================================\n\n")
  
  cat(sprintf("%-20s | %5s | %10s | %10s | %10s | %10s\n",
              "Stratum", "N", "Unadj LEPR", "Adj LEPR", "LEP beta", "Resid LEPR"))
  cat(paste(rep("-", 80), collapse=""), "\n")
  
  strata <- unique(df_both$glycemic_status)
  strata <- strata[!is.na(strata)]
  
  for (s in sort(strata)) {
    ds <- df_both[glycemic_status == s, ]
    if (nrow(ds) < 50) next
    
    ds$BMI_z <- standardize(ds$BMI)
    ds$LEPR_z <- standardize(ds$LEPR)
    ds$LEP_z <- standardize(ds$LEP)
    
    # Unadjusted
    m_unadj <- lm(BMI_z ~ LEPR_z + age + sex, data = ds)
    b_unadj <- coef(m_unadj)["LEPR_z"]
    
    # Adjusted for LEP
    m_adj <- lm(BMI_z ~ LEPR_z + LEP_z + age + sex, data = ds)
    b_adj <- coef(m_adj)["LEPR_z"]
    b_lep <- coef(m_adj)["LEP_z"]
    
    # Residualized
    lm_resid <- lm(LEPR_z ~ LEP_z, data = ds)
    ds$LEPR_resid_z <- standardize(residuals(lm_resid))
    m_resid <- lm(BMI_z ~ LEPR_resid_z + age + sex, data = ds)
    b_resid <- coef(m_resid)["LEPR_resid_z"]
    
    cat(sprintf("%-20s | %5d | %10.4f | %10.4f | %10.4f | %10.4f\n",
                s, nrow(ds), b_unadj, b_adj, b_lep, b_resid))
  }
}

# ============================================================================
# 5. COMPARISON WITH MR
# ============================================================================
cat("\n\n================================================================\n")
cat("5. GRAND COMPARISON TABLE\n")
cat("================================================================\n\n")

cat(sprintf("%-40s | %10s | %10s | %10s\n", "Estimate", "Beta", "SE", "P"))
cat(paste(rep("-", 80), collapse=""), "\n")

# Unadjusted
cat(sprintf("%-40s | %10.4f | %10.4f | %10.2e\n",
            "Observational (unadjusted)",
            coef(m1)["LEPR_z"], summary(m1)$coefficients["LEPR_z", 2],
            summary(m1)$coefficients["LEPR_z", 4]))

# Adjusted for LEP
cat(sprintf("%-40s | %10.4f | %10.4f | %10.2e\n",
            "Observational (adj. for LEP)",
            coef(m2)["LEPR_z"], summary(m2)$coefficients["LEPR_z", 2],
            summary(m2)$coefficients["LEPR_z", 4]))

# Residualized
cat(sprintf("%-40s | %10.4f | %10.4f | %10.2e\n",
            "Observational (LEPR resid. on LEP)",
            coef(m3)["LEPR_resid_z"], summary(m3)$coefficients["LEPR_resid_z", 2],
            summary(m3)$coefficients["LEPR_resid_z", 4]))

cat(paste(rep("-", 80), collapse=""), "\n")

# MR estimates (from TwoSampleMR)
cat(sprintf("%-40s | %10.4f | %10.4f | %10.4f\n",
            "MR cis-Wald (Full, rs61779759)", 0.1171, 0.0549, 0.0330))
cat(sprintf("%-40s | %10.4f | %10.4f | %10.4f\n",
            "MR cis-Wald (Normoglycemic)", 0.0919, 0.0570, 0.1073))

cat("\n================================================================\n")
cat("INTERPRETATION\n")
cat("================================================================\n\n")

unadj_sign <- sign(coef(m1)["LEPR_z"])
adj_sign <- sign(coef(m2)["LEPR_z"])
mr_sign <- 1  # positive

cat(sprintf("Unadjusted LEPR-BMI:   %s (beta = %.4f)\n",
            ifelse(unadj_sign > 0, "POSITIVE", "NEGATIVE"), coef(m1)["LEPR_z"]))
cat(sprintf("LEP-adjusted LEPR-BMI: %s (beta = %.4f)\n",
            ifelse(adj_sign > 0, "POSITIVE", "NEGATIVE"), coef(m2)["LEPR_z"]))
cat(sprintf("MR cis-Wald:           POSITIVE (beta = 0.1171)\n"))

if (unadj_sign != adj_sign) {
  cat("\n=> SIGN FLIP! Adjusting for LEP reverses the LEPR-BMI association.\n")
  cat("   This confirms reverse causation via the BMI->LEP->LEPR pathway.\n")
  if (adj_sign == mr_sign) {
    cat("   The LEP-adjusted observational estimate is now CONCORDANT with MR.\n")
  }
} else if (abs(coef(m2)["LEPR_z"]) < abs(coef(m1)["LEPR_z"])) {
  cat("\n=> Adjusting for LEP ATTENUATES the LEPR-BMI association.\n")
  cat("   Partial confounding via the BMI->LEP->LEPR pathway confirmed.\n")
}

cat("\n================================================================\n")
cat("Done.\n")
cat("================================================================\n")
