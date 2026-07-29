#!/usr/bin/env Rscript
# Glycemic stratum × metabolic hallmark: nested protein counts
# UKB reverse observational significance: GLOBAL Bonferroni only
#   p < 0.05 / N_total, where N_total = sum of M over all (phenotype × stratum) slices
#   for BMI, HbA1c, TRIG_HDL × (non_T2D, prediabetes, T2D).
# (b) STEP trial P < 0.05 (unchanged from funnel_table_obs_step_mr.R)
# (c) Reverse MR: nominal Rev IVW P < 0.05

.libPaths(c("/n/groups/patel/sivateja/R_libs", .libPaths()))
suppressPackageStartupMessages(library(dplyr))

pewas_dir <- "/n/groups/patel/sivateja/UKB/PEWAS_results"
mr_dir <- "/n/groups/patel/sivateja/regenie_pipeline/results/twosampleMR/supplemental_tables"
out_csv <- file.path(
  mr_dir,
  "Table_Funnel_ReverseObs_STEP_MR_by_stratum_hallmark_global_Bonferroni.csv"
)

traits <- c("BMI", "HbA1c", "TRIG_HDL_RATIO")
trait_labels <- c("BMI", "HbA1c", "TG/HDL")
strata <- c("non_T2D", "prediabetes", "T2D")
strata_labels <- c("Normoglycemic", "Prediabetes", "T2D")

step1 <- read.csv("/n/groups/patel/sivateja/STEP1_merged_results.csv")
step2 <- read.csv("/n/groups/patel/sivateja/STEP2_merged_results.csv")
rev_obs <- read.csv(file.path(pewas_dir, "Reverse_PEWAS_all_phenotypes_all_strata.csv"))
prot_map <- step1 %>% dplyr::select(Exposure, code) %>% distinct()

# ---- Global Bonferroni: one family = all reverse PEWAS tests (3 traits × 3 strata) ----
N_global <- 0L
for (pheno in traits) {
  for (st in strata) {
    rev_sub <- rev_obs %>%
      filter(Phenotype == pheno, Stratum == st) %>%
      left_join(prot_map, by = c("Protein" = "Exposure")) %>%
      filter(!is.na(code))
    N_global <- N_global + n_distinct(rev_sub$code)
  }
}
global_bonf <- 0.05 / max(N_global, 1L)

cat("Global reverse PEWAS tests (N_total):", N_global, "\n")
cat("Global Bonferroni threshold (0.05/N_total):", format(global_bonf, scientific = TRUE), "\n\n")

mr_file <- function(pheno) {
  fn <- sprintf(
    "Table_Bidirectional_MR_%s_Full.csv",
    ifelse(pheno == "TRIG_HDL_RATIO", "TRIG_HDL_RATIO", pheno)
  )
  file.path(mr_dir, fn)
}

load_mr_flags <- function(pheno) {
  m <- read.csv(mr_file(pheno), check.names = FALSE)
  rev_p <- suppressWarnings(as.numeric(m[["Rev: IVW P-value"]]))
  tibble(
    code = m$Protein,
    rev_mr_sig = !is.na(rev_p) & rev_p < 0.05,
    fwd_mr_bonf = m[["Fwd: Bonferroni Sig"]] == "Yes"
  ) %>%
    distinct(code, .keep_all = TRUE)
}

build_step_merged <- function(pheno, subgroup) {
  s1 <- step1 %>%
    filter(Phenotype == pheno, Subgroup == subgroup) %>%
    dplyr::select(code, effect_size, pvalue) %>%
    group_by(code) %>% arrange(pvalue) %>% slice_head(n = 1) %>% ungroup() %>%
    rename(Estimate_STEP1 = effect_size, Pval_STEP1 = pvalue)
  s2 <- step2 %>%
    filter(Phenotype == pheno, Subgroup == subgroup) %>%
    dplyr::select(code, effect_size, pvalue) %>%
    group_by(code) %>% arrange(pvalue) %>% slice_head(n = 1) %>% ungroup() %>%
    rename(Estimate_STEP2 = effect_size, Pval_STEP2 = pvalue)
  full_join(s1, s2, by = "code")
}

trial_pass_vec <- function(pheno, P1, P2) {
  if (pheno == "BMI") return(!is.na(P1) & P1 < 0.05)
  if (pheno == "HbA1c") return(!is.na(P2) & P2 < 0.05)
  pm <- suppressWarnings(pmin(P1, P2, na.rm = TRUE))
  ok <- pm < 0.05
  ok[is.infinite(pm)] <- FALSE
  ok
}

rows_out <- list()

for (ti in seq_along(traits)) {
  pheno <- traits[ti]
  mr <- load_mr_flags(pheno)
  step_merged <- build_step_merged(pheno, "non T2D")

  for (si in seq_along(strata)) {
    st <- strata[si]
    rev_sub <- rev_obs %>%
      filter(Phenotype == pheno, Stratum == st) %>%
      left_join(prot_map, by = c("Protein" = "Exposure")) %>%
      filter(!is.na(code))

    n_m <- n_distinct(rev_sub$code)

    rev_sub <- rev_sub %>%
      group_by(code) %>%
      slice_head(n = 1) %>%
      ungroup()

    rev_sub <- rev_sub %>%
      mutate(
        obs_sig_global_bonf = p.value < global_bonf
      )

    n_obs <- sum(rev_sub$obs_sig_global_bonf, na.rm = TRUE)

    merged <- rev_sub %>%
      filter(obs_sig_global_bonf) %>%
      inner_join(step_merged, by = "code") %>%
      mutate(step_ok = trial_pass_vec(pheno, Pval_STEP1, Pval_STEP2))

    n_step <- sum(merged$step_ok, na.rm = TRUE)

    merged_mr <- merged %>%
      filter(step_ok) %>%
      left_join(mr, by = "code")

    n_mr <- sum(merged_mr$rev_mr_sig %in% TRUE, na.rm = TRUE)

    rows_out[[length(rows_out) + 1L]] <- tibble(
      Glycemic_stratum = strata_labels[si],
      Glycemic_code = st,
      Hallmark = trait_labels[ti],
      Phenotype_code = pheno,
      M_proteins_tested_reverse_obs_slice = n_m,
      N_global_reverse_obs_tests = N_global,
      Global_Bonferroni_threshold = global_bonf,
      N_UKB_reverse_obs_global_Bonf = n_obs,
      N_after_STEP_p_lt_0.05 = n_step,
      N_after_reverse_MR_IVW_P_lt_0.05 = n_mr,
      STEP_arm_note = case_when(
        pheno == "BMI" ~ "STEP 1 P<0.05",
        pheno == "HbA1c" ~ "STEP 2 P<0.05",
        TRUE ~ "min(STEP1,STEP2) P<0.05"
      ),
      MR_note = paste0(
        "Rev MR: nominal IVW P < 0.05 (not Bonferroni/FDR); MR not stratified by glycemia. ",
        "UKB reverse obs: global Bonferroni p < 0.05/", N_global, " across ",
        "3 hallmarks × 3 strata (BH-FDR not used in this table)."
      )
    )
  }
}

out <- bind_rows(rows_out)
write.csv(out, out_csv, row.names = FALSE)
cat("Wrote:", out_csv, "\n")
print(out, n = 100)
