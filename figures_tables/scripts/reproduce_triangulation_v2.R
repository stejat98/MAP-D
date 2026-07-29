.libPaths(c("/path/to/home/R/x86_64-pc-linux-gnu-library/4.4",
            "/path/to/project/R_libs", .libPaths()))
library(dplyr)

cat("================================================================\n")
cat("TRIANGULATION v2: BMI->STEP1, HbA1c->STEP2, TG/HDL->both\n")
cat("================================================================\n\n")

step1 <- read.csv("/path/to/project/STEP1_merged_results.csv")
step2 <- read.csv("/path/to/project/STEP2_merged_results.csv")
rev_obs <- read.csv("/path/to/project/UKB/PEWAS_results/Reverse_PEWAS_all_phenotypes_all_strata.csv")
prot_map <- step1 %>% dplyr::select(Exposure, code) %>% distinct()

mr_dir <- "/path/to/project/regenie_pipeline/results/twosampleMR/supplemental_tables"

# trial_mode: "step1" = STEP1 only, "step2" = STEP2 only, "both" = average
run_triangulation <- function(pheno, pheno_label, stratum, mr_file, step_subgroup, trial_mode) {

  cat("================================================================\n")
  cat(sprintf("  %s  (trial_mode = %s)\n", pheno_label, trial_mode))
  cat("================================================================\n\n")

  # A: Reverse observational
  rev_sub <- rev_obs %>%
    filter(Phenotype == pheno, Stratum == stratum) %>%
    left_join(prot_map, by = c("Protein" = "Exposure")) %>%
    filter(!is.na(code))
  cat(sprintf("STEP A: %d reverse obs proteins (Phenotype='%s', Stratum='%s')\n", nrow(rev_sub), pheno, stratum))

  # B: Bidirectional MR
  bidir <- read.csv(mr_file, check.names = FALSE)
  rev_beta_col <- grep("Rev.*IVW.*Beta", names(bidir), value = TRUE)[1]
  rev_p_col <- grep("Rev.*IVW.*P", names(bidir), value = TRUE)[1]
  bidir$Rev_Beta <- as.numeric(bidir[[rev_beta_col]])
  bidir$Rev_P <- as.numeric(bidir[[rev_p_col]])
  rev_mr_sig <- bidir %>% filter(Rev_P < 0.05) %>%
    dplyr::select(Protein, Rev_Beta) %>%
    rename(code = Protein, Beta_revMR = Rev_Beta)
  cat(sprintf("STEP B: %d proteins with nominally sig reverse MR (P<0.05)\n", nrow(rev_mr_sig)))

  # C: STEP trial effects
  step1_sub <- step1 %>%
    filter(Phenotype == pheno, Subgroup == step_subgroup) %>%
    dplyr::select(code, effect_size, pvalue) %>%
    group_by(code) %>% arrange(pvalue) %>% slice_head(n = 1) %>% ungroup() %>%
    rename(Estimate_STEP1 = effect_size, Pval_STEP1 = pvalue)

  step2_sub <- step2 %>%
    filter(Phenotype == pheno, Subgroup == step_subgroup) %>%
    dplyr::select(code, effect_size, pvalue) %>%
    group_by(code) %>% arrange(pvalue) %>% slice_head(n = 1) %>% ungroup() %>%
    rename(Estimate_STEP2 = effect_size, Pval_STEP2 = pvalue)

  step_merged <- full_join(step1_sub, step2_sub, by = "code")

  if (trial_mode == "step1") {
    step_merged$ESTIMATE_TRIAL <- step_merged$Estimate_STEP1
    step_merged$Pval_TRIAL <- step_merged$Pval_STEP1
    trial_label <- "STEP 1"
  } else if (trial_mode == "step2") {
    step_merged$ESTIMATE_TRIAL <- step_merged$Estimate_STEP2
    step_merged$Pval_TRIAL <- step_merged$Pval_STEP2
    trial_label <- "STEP 2"
  } else {
    step_merged$ESTIMATE_TRIAL <- rowMeans(step_merged[, c("Estimate_STEP1", "Estimate_STEP2")], na.rm = TRUE)
    step_merged$Pval_TRIAL <- pmin(step_merged$Pval_STEP1, step_merged$Pval_STEP2, na.rm = TRUE)
    trial_label <- "STEP 1/2 (min P, mean effect)"
  }
  cat(sprintf("STEP C: Using %s for trial effect + filtering\n", trial_label))

  # D: Merge + filter
  rev_for_merge <- rev_sub %>%
    dplyr::select(code, Protein, Beta_obs_rev = estimate, Pval_obs_rev = p.value)
  n_rev <- nrow(rev_for_merge)
  bonf <- 0.05 / n_rev
  cat(sprintf("STEP D: Bonferroni = 0.05/%d = %.2e\n", n_rev, bonf))

  merged <- step_merged %>% inner_join(rev_for_merge, by = "code")
  cat(sprintf("  inner_join:                  %d proteins\n", nrow(merged)))

  after_bonf <- merged %>% filter(Pval_obs_rev < bonf)
  cat(sprintf("  Pval_obs_rev < Bonf:         %d proteins\n", nrow(after_bonf)))

  after_trial <- after_bonf %>% filter(Pval_TRIAL < 0.05)
  cat(sprintf("  Pval_%s < 0.05:        %d proteins  *** REPORTED N ***\n", trial_label, nrow(after_trial)))

  plot_df <- after_trial %>%
    left_join(rev_mr_sig, by = "code") %>%
    group_by(code) %>% slice_head(n = 1) %>% ungroup()

  plot_df$rev_mr_sig <- !is.na(plot_df$Beta_revMR)
  plot_df$concordant <- plot_df$rev_mr_sig &
    (sign(plot_df$Beta_revMR) == sign(plot_df$Beta_obs_rev))

  cat(sprintf("  Nominally sig reverse MR:    %d proteins\n", sum(plot_df$rev_mr_sig)))
  cat(sprintf("  CONCORDANT (triple triang):  %d proteins  ***\n\n", sum(plot_df$concordant)))

  # Anti-correlation check
  n_anti <- sum(sign(plot_df$Beta_obs_rev) != sign(plot_df$ESTIMATE_TRIAL), na.rm = TRUE)
  n_same <- sum(sign(plot_df$Beta_obs_rev) == sign(plot_df$ESTIMATE_TRIAL), na.rm = TRUE)
  cat(sprintf("  Direction: %d anti-correlated (obs vs trial), %d same-direction\n\n", n_anti, n_same))

  # Top concordant
  concordant_df <- plot_df %>%
    filter(concordant) %>%
    mutate(combined_effect = abs(Beta_obs_rev) + abs(ESTIMATE_TRIAL) + abs(Beta_revMR)) %>%
    arrange(desc(combined_effect))

  cat(sprintf("TOP %s CONCORDANT PROTEINS:\n", pheno_label))
  cat("---------------------------------------------------------------\n")
  top <- concordant_df %>%
    select(code, Beta_obs_rev, ESTIMATE_TRIAL, Pval_TRIAL, Beta_revMR, combined_effect) %>%
    head(10)
  print(as.data.frame(top), row.names = FALSE)

  # Quadrant breakdown
  cat(sprintf("\nQuadrant breakdown (all %d):\n", nrow(plot_df)))
  plot_df$quadrant <- paste0(
    ifelse(plot_df$Beta_obs_rev > 0, "Inc_UKB", "Dec_UKB"), " / ",
    ifelse(plot_df$ESTIMATE_TRIAL > 0, "Inc_TRIAL", "Dec_TRIAL"))
  quad_tbl <- plot_df %>%
    group_by(quadrant) %>%
    summarise(n = n(), n_concordant = sum(concordant), .groups = "drop")
  print(as.data.frame(quad_tbl), row.names = FALSE)
  cat("\n\n")

  return(list(plot_df = plot_df, concordant_df = concordant_df))
}

# ---- BMI: STEP 1 only ----
bmi_res <- run_triangulation(
  "BMI", "BMI", "non_T2D",
  file.path(mr_dir, "Table_Bidirectional_MR_BMI_Full.csv"),
  "non T2D", trial_mode = "step1")

# ---- HbA1c: STEP 2 only ----
hba1c_res <- run_triangulation(
  "HbA1c", "HbA1c", "non_T2D",
  file.path(mr_dir, "Table_Bidirectional_MR_HbA1c_Full.csv"),
  "non T2D", trial_mode = "step2")

# ---- TG/HDL: both STEP 1 and STEP 2 ----
trig_res_both <- run_triangulation(
  "TRIG_HDL_RATIO", "TG/HDL-C Ratio (STEP 1+2 combined)", "non_T2D",
  file.path(mr_dir, "Table_Bidirectional_MR_TRIG_HDL_RATIO_Full.csv"),
  "non T2D", trial_mode = "both")

# Also show TG/HDL with each trial separately for comparison
trig_res_s1 <- run_triangulation(
  "TRIG_HDL_RATIO", "TG/HDL-C Ratio (STEP 1 only)", "non_T2D",
  file.path(mr_dir, "Table_Bidirectional_MR_TRIG_HDL_RATIO_Full.csv"),
  "non T2D", trial_mode = "step1")

trig_res_s2 <- run_triangulation(
  "TRIG_HDL_RATIO", "TG/HDL-C Ratio (STEP 2 only)", "non_T2D",
  file.path(mr_dir, "Table_Bidirectional_MR_TRIG_HDL_RATIO_Full.csv"),
  "non T2D", trial_mode = "step2")

cat("================================================================\n")
cat("FINAL SUMMARY\n")
cat("================================================================\n")
cat(sprintf("BMI (STEP 1):            %d Bonf+nom, %d concordant\n",
            nrow(bmi_res$plot_df), nrow(bmi_res$concordant_df)))
cat(sprintf("HbA1c (STEP 2):          %d Bonf+nom, %d concordant\n",
            nrow(hba1c_res$plot_df), nrow(hba1c_res$concordant_df)))
cat(sprintf("TG/HDL (STEP 1+2 comb):  %d Bonf+nom, %d concordant\n",
            nrow(trig_res_both$plot_df), nrow(trig_res_both$concordant_df)))
cat(sprintf("TG/HDL (STEP 1 only):    %d Bonf+nom, %d concordant\n",
            nrow(trig_res_s1$plot_df), nrow(trig_res_s1$concordant_df)))
cat(sprintf("TG/HDL (STEP 2 only):    %d Bonf+nom, %d concordant\n",
            nrow(trig_res_s2$plot_df), nrow(trig_res_s2$concordant_df)))
