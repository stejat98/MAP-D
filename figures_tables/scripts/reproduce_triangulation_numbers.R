.libPaths(c("/home/st320/R/x86_64-pc-linux-gnu-library/4.4",
            "/n/groups/patel/sivateja/R_libs", .libPaths()))
library(dplyr)

cat("================================================================\n")
cat("REPRODUCING TRIANGULATION NUMBERS - FULL AUDIT TRAIL\n")
cat("================================================================\n\n")

cat("INPUT FILES:\n")
cat("  1. STEP 1 trial:    /n/groups/patel/sivateja/STEP1_merged_results.csv\n")
cat("  2. STEP 2 trial:    /n/groups/patel/sivateja/STEP2_merged_results.csv\n")
cat("  3. Reverse obs:     /n/groups/patel/sivateja/UKB/PEWAS_results/Reverse_PEWAS_all_phenotypes_all_strata.csv\n")
cat("  4. Bidir MR (BMI):  .../Table_Bidirectional_MR_BMI_Full.csv\n")
cat("  5. Bidir MR (HbA1c):.../Table_Bidirectional_MR_HbA1c_Full.csv\n")
cat("  6. Bidir MR (TG/HDL):.../Table_Bidirectional_MR_TRIG_HDL_RATIO_Full.csv\n\n")

step1 <- read.csv("/n/groups/patel/sivateja/STEP1_merged_results.csv")
step2 <- read.csv("/n/groups/patel/sivateja/STEP2_merged_results.csv")
rev_obs <- read.csv("/n/groups/patel/sivateja/UKB/PEWAS_results/Reverse_PEWAS_all_phenotypes_all_strata.csv")
prot_map <- step1 %>% dplyr::select(Exposure, code) %>% distinct()

mr_dir <- "/n/groups/patel/sivateja/regenie_pipeline/results/twosampleMR/supplemental_tables"

cat("Raw table sizes:\n")
cat(sprintf("  STEP1: %d rows x %d cols\n", nrow(step1), ncol(step1)))
cat(sprintf("  STEP2: %d rows x %d cols\n", nrow(step2), ncol(step2)))
cat(sprintf("  Reverse obs: %d rows x %d cols\n", nrow(rev_obs), ncol(rev_obs)))
cat(sprintf("  Protein map (Exposure->code): %d unique mappings\n\n", nrow(prot_map)))

run_triangulation <- function(pheno, pheno_label, stratum, mr_file, step_subgroup) {

  cat("================================================================\n")
  cat(sprintf("  %s\n", pheno_label))
  cat("================================================================\n\n")

  # STEP A: Reverse observational
  cat("STEP A: Reverse observational (UKB)\n")
  cat(sprintf("  Filter: Phenotype == '%s', Stratum == '%s'\n", pheno, stratum))
  rev_sub <- rev_obs %>%
    filter(Phenotype == pheno, Stratum == stratum) %>%
    left_join(prot_map, by = c("Protein" = "Exposure")) %>%
    filter(!is.na(code))
  cat(sprintf("  Result: %d proteins with valid code mapping\n\n", nrow(rev_sub)))

  # STEP B: Bidirectional MR
  cat("STEP B: Bidirectional MR\n")
  cat(sprintf("  File: %s\n", basename(mr_file)))
  bidir <- read.csv(mr_file, check.names = FALSE)
  cat(sprintf("  Total rows: %d\n", nrow(bidir)))

  rev_beta_col <- grep("Rev.*IVW.*Beta|Rev.*Beta", names(bidir), value = TRUE)[1]
  rev_p_col <- grep("Rev.*IVW.*P|Rev.*P.value", names(bidir), value = TRUE)[1]
  cat(sprintf("  Rev beta column: '%s'\n", rev_beta_col))
  cat(sprintf("  Rev P column:    '%s'\n", rev_p_col))

  bidir$Rev_Beta <- as.numeric(bidir[[rev_beta_col]])
  bidir$Rev_P <- as.numeric(bidir[[rev_p_col]])

  rev_mr_sig <- bidir %>% filter(Rev_P < 0.05) %>%
    dplyr::select(Protein, Rev_Beta) %>%
    rename(code = Protein, Beta_revMR = Rev_Beta)
  cat(sprintf("  Filter: Rev_P < 0.05 => %d proteins with nominally sig reverse MR\n\n", nrow(rev_mr_sig)))

  # STEP C: STEP trial effects
  cat("STEP C: STEP 1/2 trial effects\n")
  cat(sprintf("  Filter: Phenotype == '%s', Subgroup == '%s'\n", pheno, step_subgroup))

  step1_sub <- step1 %>%
    filter(Phenotype == pheno, Subgroup == step_subgroup) %>%
    dplyr::select(code, effect_size, pvalue) %>%
    group_by(code) %>% arrange(pvalue) %>% slice_head(n = 1) %>% ungroup() %>%
    rename(Estimate_STEP1 = effect_size, Pval_STEP1 = pvalue)
  cat(sprintf("  STEP 1 proteins (after dedup): %d\n", nrow(step1_sub)))

  step2_sub <- step2 %>%
    filter(Phenotype == pheno, Subgroup == step_subgroup) %>%
    dplyr::select(code, effect_size, pvalue) %>%
    group_by(code) %>% arrange(pvalue) %>% slice_head(n = 1) %>% ungroup() %>%
    rename(Estimate_STEP2 = effect_size, Pval_STEP2 = pvalue)
  cat(sprintf("  STEP 2 proteins (after dedup): %d\n", nrow(step2_sub)))

  step_merged <- full_join(step1_sub, step2_sub, by = "code")
  step_merged$ESTIMATE_STEP_MEAN <- rowMeans(step_merged[, c("Estimate_STEP1", "Estimate_STEP2")], na.rm = TRUE)
  cat(sprintf("  Merged STEP 1+2: %d unique proteins\n\n", nrow(step_merged)))

  # STEP D: Merge and filter
  cat("STEP D: Merge + sequential filtering\n")
  rev_for_merge <- rev_sub %>%
    dplyr::select(code, Protein, Beta_obs_rev = estimate, Pval_obs_rev = p.value)

  n_rev <- nrow(rev_for_merge)
  bonf <- 0.05 / n_rev
  cat(sprintf("  Bonferroni threshold: 0.05 / %d = %.2e\n", n_rev, bonf))

  merged <- step_merged %>% inner_join(rev_for_merge, by = "code")
  cat(sprintf("  After inner_join (STEP x rev_obs):  %d proteins\n", nrow(merged)))

  after_bonf <- merged %>% filter(Pval_obs_rev < bonf)
  cat(sprintf("  After Pval_obs_rev < Bonf:          %d proteins\n", nrow(after_bonf)))

  after_step1 <- after_bonf %>% filter(Pval_STEP1 < 0.05)
  cat(sprintf("  After Pval_STEP1 < 0.05:            %d proteins  *** REPORTED N ***\n", nrow(after_step1)))

  plot_df <- after_step1 %>%
    left_join(rev_mr_sig, by = "code") %>%
    group_by(code) %>% slice_head(n = 1) %>% ungroup()

  plot_df$rev_mr_sig <- !is.na(plot_df$Beta_revMR)
  plot_df$concordant <- plot_df$rev_mr_sig &
    (sign(plot_df$Beta_revMR) == sign(plot_df$Beta_obs_rev))

  cat(sprintf("  With nominally sig reverse MR:      %d proteins\n", sum(plot_df$rev_mr_sig)))
  cat(sprintf("  CONCORDANT (obs+STEP+MR agree):     %d proteins  *** TRIPLE TRIANGULATION ***\n\n", sum(plot_df$concordant)))

  # STEP E: Top concordant hits
  concordant_df <- plot_df %>%
    filter(concordant) %>%
    mutate(combined_effect = abs(Beta_obs_rev) + abs(ESTIMATE_STEP_MEAN) + abs(Beta_revMR)) %>%
    arrange(desc(combined_effect))

  cat(sprintf("TOP %s CONCORDANT PROTEINS (by combined |obs|+|STEP|+|MR|):\n", pheno_label))
  cat("---------------------------------------------------------------\n")
  top <- concordant_df %>%
    select(code, Beta_obs_rev, ESTIMATE_STEP_MEAN, Beta_revMR, combined_effect) %>%
    head(10)
  print(as.data.frame(top), row.names = FALSE)

  cat(sprintf("\nQuadrant breakdown (all %d proteins in set):\n", nrow(plot_df)))
  plot_df$quadrant <- paste0(
    ifelse(plot_df$Beta_obs_rev > 0, "Inc_UKB", "Dec_UKB"), " / ",
    ifelse(plot_df$ESTIMATE_STEP_MEAN > 0, "Inc_STEP", "Dec_STEP"))
  quad_tbl <- plot_df %>%
    group_by(quadrant) %>%
    summarise(n_total = n(), n_concordant_MR = sum(concordant), .groups = "drop")
  print(as.data.frame(quad_tbl), row.names = FALSE)
  cat("\n\n")

  return(concordant_df)
}

bmi_conc <- run_triangulation(
  "BMI", "BMI", "non_T2D",
  file.path(mr_dir, "Table_Bidirectional_MR_BMI_Full.csv"),
  "non T2D")

hba1c_conc <- run_triangulation(
  "HbA1c", "HbA1c", "non_T2D",
  file.path(mr_dir, "Table_Bidirectional_MR_HbA1c_Full.csv"),
  "non T2D")

trig_conc <- run_triangulation(
  "TRIG_HDL_RATIO", "TG/HDL-C Ratio", "non_T2D",
  file.path(mr_dir, "Table_Bidirectional_MR_TRIG_HDL_RATIO_Full.csv"),
  "non T2D")

cat("================================================================\n")
cat("FINAL SUMMARY - MANUSCRIPT NUMBERS\n")
cat("================================================================\n")
cat(sprintf("BMI:          %d concordant triple-triangulated proteins\n", nrow(bmi_conc)))
cat(sprintf("HbA1c:        %d concordant triple-triangulated proteins\n", nrow(hba1c_conc)))
cat(sprintf("TG/HDL Ratio: %d concordant triple-triangulated proteins\n", nrow(trig_conc)))
