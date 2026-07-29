#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# 1) Table_Funnel_Manuscript_Reverse_Triangulation.csv — WIDE funnel (rows = stages,
#    columns = BMI or HbA1c × glycemic stratum). TG/HDL → supplemental CSV.
#    Also writes Table_Funnel_Manuscript_Reverse_Triangulation_long.csv (all hallmarks, long).
#    Matches Fig. 5 / validation_concordance_all_phenotypes.R:
#    - UKB strata: non_T2D, prediabetes, T2D
#    - Reverse observational: Bonferroni only (P < 0.05 / M per stratum), OR optional
#      global Bonferroni (P < 0.05 / N) with N = all reverse PEWAS tests across
#      BMI, HbA1c, TG/HDL × three glycemic strata (supplemental wide tables).
#    - STEP: BMI → STEP 1; HbA1c → STEP 2; TG/HDL → min(STEP1,STEP2) P<0.05
#    - Reverse MR: IVW P<0.05; concordant = sign(rev MR)=sign(UKB rev obs); STEP gates only
#
# 2) Table_Supplemental_HbA1c_Forward_Funnel_by_stratum.csv
#    Forward observational (HbA1c ~ protein) + STEP 2 + forward MR (Bonferroni or FDR P<0.05)
# -----------------------------------------------------------------------------

.libPaths(c("/n/groups/patel/sivateja/R_libs", .libPaths()))
suppressPackageStartupMessages(library(dplyr))

pewas_dir <- "/n/groups/patel/sivateja/UKB/PEWAS_results"
mr_dir <- "/n/groups/patel/sivateja/regenie_pipeline/results/twosampleMR/supplemental_tables"
out1 <- file.path(mr_dir, "Table_Funnel_Manuscript_Reverse_Triangulation.csv")
out1_long <- file.path(mr_dir, "Table_Funnel_Manuscript_Reverse_Triangulation_long.csv")
out_trig <- file.path(mr_dir, "Table_Funnel_Manuscript_Reverse_Triangulation_TG_HDL_supplemental.csv")
out1_global <- file.path(mr_dir, "Table_Funnel_Manuscript_Reverse_Triangulation_global_Bonferroni.csv")
out_trig_global <- file.path(
  mr_dir,
  "Table_Funnel_Manuscript_Reverse_Triangulation_TG_HDL_supplemental_global_Bonferroni.csv"
)
out_conc <- file.path(mr_dir, "Table_Concordant_Proteins_For_Drug_Overlap.csv")
out2 <- file.path(mr_dir, "Table_Supplemental_HbA1c_Forward_Funnel_by_stratum.csv")

step1 <- read.csv("/n/groups/patel/sivateja/STEP1_merged_results.csv")
step2 <- read.csv("/n/groups/patel/sivateja/STEP2_merged_results.csv")
rev_obs <- read.csv(file.path(pewas_dir, "Reverse_PEWAS_all_phenotypes_all_strata.csv"))
prot_map <- step1 %>% dplyr::select(Exposure, code) %>% distinct()

# ---- Shared: bidirectional MR file reader (nominal reverse MR IVW P < 0.05) ----
read_rev_mr_nominal <- function(pheno) {
  fn <- sprintf("Table_Bidirectional_MR_%s_Full.csv",
                ifelse(pheno == "TRIG_HDL_RATIO", "TRIG_HDL_RATIO", pheno))
  m <- read.csv(file.path(mr_dir, fn), check.names = FALSE)
  rev_p_col <- grep("Rev.*IVW.*P|Rev.*P.value", names(m), value = TRUE)[1]
  rev_beta_col <- grep("Rev.*IVW.*Beta|Rev.*Beta", names(m), value = TRUE)[1]
  tibble(
    code = m$Protein,
    Rev_P = as.numeric(m[[rev_p_col]]),
    Rev_Beta = as.numeric(m[[rev_beta_col]])
  ) %>%
    filter(!is.na(Rev_P), Rev_P < 0.05) %>%
    dplyr::select(code, Beta_revMR = Rev_Beta) %>%
    distinct(code, .keep_all = TRUE)
}

manuscript_reverse_row <- function(
    pheno,
    pheno_label,
    trial_mode,
    stratum,
    stratum_label,
    bonf_override = NULL,
    bonferroni_scope = "per_stratum") {
  step_subgroup <- "non T2D"

  rev_sub <- rev_obs %>%
    filter(Phenotype == pheno, Stratum == stratum) %>%
    left_join(prot_map, by = c("Protein" = "Exposure")) %>%
    filter(!is.na(code))

  rev_for_merge <- rev_sub %>%
    dplyr::select(code, Beta_obs_rev = estimate, Pval_obs_rev = p.value)

  n_rev <- nrow(rev_for_merge)
  bonf <- if (is.null(bonf_override)) {
    0.05 / max(n_rev, 1L)
  } else {
    bonf_override
  }

  n_bonf_only <- rev_for_merge %>%
    distinct(code, .keep_all = TRUE) %>%
    filter(Pval_obs_rev < bonf) %>%
    nrow()

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
  } else if (trial_mode == "step2") {
    step_merged$ESTIMATE_TRIAL <- step_merged$Estimate_STEP2
    step_merged$Pval_TRIAL <- step_merged$Pval_STEP2
  } else {
    step_merged$ESTIMATE_TRIAL <- rowMeans(step_merged[, c("Estimate_STEP1", "Estimate_STEP2")], na.rm = TRUE)
    step_merged$Pval_TRIAL <- pmin(step_merged$Pval_STEP1, step_merged$Pval_STEP2, na.rm = TRUE)
  }

  rev_mr_sig <- read_rev_mr_nominal(pheno)

  plot_df <- step_merged %>%
    inner_join(rev_for_merge, by = "code") %>%
    filter(Pval_obs_rev < bonf, Pval_TRIAL < 0.05) %>%
    left_join(rev_mr_sig, by = "code") %>%
    group_by(code) %>% slice_head(n = 1) %>% ungroup()

  # Concordance is ONLY reverse-MR vs UKB reverse-obs β (same sign). STEP is not in this
  # comparison; it only gates membership (Bonferroni rev obs + nominal STEP P).
  plot_df$rev_mr_sig <- !is.na(plot_df$Beta_revMR)
  plot_df$concordant <- plot_df$rev_mr_sig & (sign(plot_df$Beta_revMR) == sign(plot_df$Beta_obs_rev))

  summary_row <- tibble(
    Hallmark = pheno_label,
    Phenotype_code = pheno,
    UKB_glycemic_stratum = stratum_label,
    UKB_Stratum_code = stratum,
    M_proteins_reverse_PEWAS = n_rev,
    Bonferroni_P_threshold = bonf,
    N_UKB_reverse_obs_Bonferroni_only = n_bonf_only,
    N_Bonferroni_and_STEP_nominal = nrow(plot_df),
    N_with_nominal_reverse_MR = sum(plot_df$rev_mr_sig),
    N_concordant_reverse_obs_and_reverseMR = sum(plot_df$concordant),
    STEP_definition = dplyr::case_when(
      trial_mode == "step1" ~ "STEP 1 P < 0.05 (non-T2D trial subgroup)",
      trial_mode == "step2" ~ "STEP 2 P < 0.05 (non-T2D trial subgroup)",
      TRUE ~ "min(STEP1, STEP2) P < 0.05 (non-T2D trial subgroup)"
    ),
    Reverse_MR_definition = "Rev IVW P < 0.05; concordant = sign(Rev MR IVW β) = sign(UKB reverse obs β). STEP not used in concordance.",
    Bonferroni_scope = bonferroni_scope,
    Manuscript_notes = paste0(
      "STEP arms match Fig. 5 (BMI→STEP1, HbA1c→STEP2, TG/HDL→min P); MR not stratified by glycemia",
      if (bonferroni_scope != "per_stratum") {
        paste0(" | Reverse-obs Bonferroni: ", bonferroni_scope)
      } else {
        ""
      }
    )
  )

  concordant_detail <- plot_df %>%
    filter(concordant) %>%
    transmute(
      Hallmark = pheno_label,
      Phenotype_code = pheno,
      UKB_glycemic_stratum = stratum_label,
      UKB_Stratum_code = stratum,
      Gene = code,
      Beta_obs_rev,
      Beta_revMR,
      ESTIMATE_TRIAL
    )

  list(summary = summary_row, concordant_proteins = concordant_detail)
}

strata_codes <- c("non_T2D", "prediabetes", "T2D")
strata_labels <- c("Normoglycemic", "Prediabetes", "T2D")

phenos <- c("BMI", "HbA1c", "TRIG_HDL_RATIO")
pheno_labels <- c("BMI", "HbA1c", "TG/HDL")
trial_modes <- c("step1", "step2", "both")

manuscript_rows <- list()
concordant_stack <- list()
for (j in seq_along(strata_codes)) {
  for (i in seq_along(phenos)) {
    res <- manuscript_reverse_row(
      phenos[i], pheno_labels[i], trial_modes[i],
      strata_codes[j], strata_labels[j]
    )
    manuscript_rows[[length(manuscript_rows) + 1L]] <- res$summary
    if (nrow(res$concordant_proteins) > 0) {
      concordant_stack[[length(concordant_stack) + 1L]] <- res$concordant_proteins
    }
  }
}
tbl_manuscript <- bind_rows(manuscript_rows)
tbl_concordant_proteins <- if (length(concordant_stack)) bind_rows(concordant_stack) else tibble()
write.csv(tbl_concordant_proteins, out_conc, row.names = FALSE)
cat("Wrote:", out_conc, "(for DrugBank directional overlap)\n")

# --- Wide “funnel down” tables: rows = nested stages, cols = Hallmark | Stratum ---
build_funnel_wide <- function(tbl) {
  tbl <- tbl %>%
    arrange(match(Phenotype_code, c("BMI", "HbA1c", "TRIG_HDL_RATIO")),
            match(UKB_Stratum_code, strata_codes))
  stages <- c(
    "1. UKB reverse obs — Bonferroni only",
    "2. + STEP trial P < 0.05 (hallmark-specific arm)",
    "3. With nominal reverse MR (IVW P < 0.05)",
    "4. Concordant: sign(rev MR IVW) = sign(UKB reverse obs)"
  )
  out <- data.frame(Funnel_stage = stages, check.names = FALSE, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(tbl))) {
    cn <- sprintf("%s | %s", tbl$Hallmark[i], tbl$UKB_glycemic_stratum[i])
    out[[cn]] <- c(
      tbl$N_UKB_reverse_obs_Bonferroni_only[i],
      tbl$N_Bonferroni_and_STEP_nominal[i],
      tbl$N_with_nominal_reverse_MR[i],
      tbl$N_concordant_reverse_obs_and_reverseMR[i]
    )
  }
  out
}

tbl_main <- tbl_manuscript %>% filter(Phenotype_code %in% c("BMI", "HbA1c"))
tbl_trig <- tbl_manuscript %>% filter(Phenotype_code == "TRIG_HDL_RATIO")

funnel_main_wide <- build_funnel_wide(tbl_main)
funnel_trig_wide <- build_funnel_wide(tbl_trig)

write.csv(funnel_main_wide, out1, row.names = FALSE)
cat("Wrote (wide funnel, BMI + HbA1c):", out1, "\n")

write.csv(tbl_manuscript, out1_long, row.names = FALSE)
cat("Wrote (long, all hallmarks × strata):", out1_long, "\n")

write.csv(funnel_trig_wide, out_trig, row.names = FALSE)
cat("Wrote (wide funnel, TG/HDL supplemental):", out_trig, "\n")

# --- Global Bonferroni (optional supplemental): one α = 0.05/N across all reverse PEWAS tests ---
# N = distinct (Protein, Phenotype, Stratum) for BMI, HbA1c, TG/HDL × three glycemic strata.
N_global <- rev_obs %>%
  dplyr::filter(Phenotype %in% c("BMI", "HbA1c", "TRIG_HDL_RATIO")) %>%
  distinct(Protein, Phenotype, Stratum) %>%
  nrow()
bonf_global <- 0.05 / max(N_global, 1L)

manuscript_rows_global <- list()
for (j in seq_along(strata_codes)) {
  for (i in seq_along(phenos)) {
    res <- manuscript_reverse_row(
      phenos[i], pheno_labels[i], trial_modes[i],
      strata_codes[j], strata_labels[j],
      bonf_override = bonf_global,
      bonferroni_scope = sprintf("global N=%s (BMI+HbA1c+TG/HDL × 3 strata)", N_global)
    )
    manuscript_rows_global[[length(manuscript_rows_global) + 1L]] <- res$summary
  }
}
tbl_manuscript_global <- bind_rows(manuscript_rows_global)
tbl_main_global <- tbl_manuscript_global %>% filter(Phenotype_code %in% c("BMI", "HbA1c"))
tbl_trig_global <- tbl_manuscript_global %>% filter(Phenotype_code == "TRIG_HDL_RATIO")
funnel_main_wide_global <- build_funnel_wide(tbl_main_global)
funnel_trig_wide_global <- build_funnel_wide(tbl_trig_global)
write.csv(funnel_main_wide_global, out1_global, row.names = FALSE)
write.csv(funnel_trig_wide_global, out_trig_global, row.names = FALSE)
cat("Wrote (wide funnel, global Bonferroni, BMI+HbA1c):", out1_global, "\n")
cat("Wrote (wide funnel, global Bonferroni, TG/HDL):", out_trig_global, "\n")
cat("Global Bonferroni: N =", N_global, ", alpha =", bonf_global, "\n")

print(tbl_manuscript, width = 200)

# =============================================================================
# Supplemental: HbA1c FORWARD funnel × glycemic stratum
# =============================================================================

pheno <- "HbA1c"
mr_file <- file.path(mr_dir, "Table_Bidirectional_MR_HbA1c_Full.csv")
mr <- read.csv(mr_file, check.names = FALSE)
mr$Fwd_P_num <- as.numeric(mr[["Fwd: P-value"]])
mr$Fwd_FDR_BH <- p.adjust(mr$Fwd_P_num, method = "BH")
mr_fwd <- tibble(
  code = mr$Protein,
  fwd_bonf = mr[["Fwd: Bonferroni Sig"]] == "Yes",
  fwd_fdr_bh = mr$Fwd_FDR_BH < 0.05 & !is.na(mr$Fwd_FDR_BH)
) %>%
  mutate(fwd_mr_sig = fwd_bonf | fwd_fdr_bh) %>%
  distinct(code, .keep_all = TRUE)

strata <- c("non_T2D", "prediabetes", "T2D")
strata_lab <- c("Normoglycemic", "Prediabetes", "T2D")
suffix <- "adj_fasting_time.RDS"

step2_hba1c <- step2 %>%
  filter(Phenotype == pheno, Subgroup == "non T2D") %>%
  dplyr::select(code, effect_size, pvalue) %>%
  group_by(code) %>% arrange(pvalue) %>% slice_head(n = 1) %>% ungroup() %>%
  rename(Est_STEP2 = effect_size, P_STEP2 = pvalue)

forward_rows <- list()

for (i in seq_along(strata)) {
  st <- strata[i]
  fn <- file.path(pewas_dir, sprintf("HbA1c_Linear_regression_proteomic_lm_results_%s_%s", st, suffix))
  if (!file.exists(fn)) {
    warning("Missing: ", fn)
    next
  }
  fwd <- readRDS(fn) %>%
    filter(grepl("scale\\(protein", term)) %>%
    dplyr::select(Exposure, P = p.value) %>%
    left_join(prot_map, by = "Exposure") %>%
    filter(!is.na(code)) %>%
    distinct(code, .keep_all = TRUE)

  M <- nrow(fwd)
  bonf <- 0.05 / max(M, 1L)
  fwd$FDR <- p.adjust(fwd$P, method = "BH")

  fwd <- fwd %>%
    mutate(
      bonf_ok = !is.na(P) & P < bonf,
      fdr_ok = !is.na(FDR) & FDR < 0.05,
      obs_sig = bonf_ok | fdr_ok
    )

  n_obs <- sum(fwd$obs_sig, na.rm = TRUE)

  merged <- fwd %>%
    filter(obs_sig) %>%
    inner_join(step2_hba1c, by = "code") %>%
    filter(P_STEP2 < 0.05)

  n_step <- nrow(merged)

  merged_mr <- merged %>%
    left_join(mr_fwd, by = "code")

  n_mr <- sum(merged_mr$fwd_mr_sig %in% TRUE, na.rm = TRUE)

  forward_rows[[length(forward_rows) + 1L]] <- tibble(
    Glycemic_stratum = strata_lab[i],
    Glycemic_code = st,
    M_proteins_forward_PEWAS = M,
    Bonferroni_P_threshold = bonf,
    N_UKB_forward_obs_Bonferroni_or_FDR = n_obs,
    N_after_STEP2_P_lt_0.05 = n_step,
    N_after_forward_MR_Bonferroni_or_BH_FDR = n_mr,
    Forward_obs_note = "HbA1c ~ protein; FDR = BH across proteins in stratum",
    STEP_note = "STEP 2 P < 0.05, non-T2D trial subgroup (same as HbA1c concordance figure)",
    Forward_MR_note = "Fwd MR Bonferroni yes OR BH-FDR on Fwd IVW P across 349 tested proteins; MR not stratified by glycemia"
  )
}

tbl_fwd <- bind_rows(forward_rows)
write.csv(tbl_fwd, out2, row.names = FALSE)
cat("\nWrote:", out2, "\n")
print(tbl_fwd, width = 200)
