#!/usr/bin/env Rscript
# =============================================================================
# Master Integrated Table: Per-protein row with all evidence layers
#   - Forward observational (Phenotype ~ Protein)
#   - Reverse observational (Protein ~ Phenotype)
#   - Forward MR (Protein -> Phenotype via cis-pQTL)
#   - Reverse MR (Phenotype -> Protein via GWAS instruments)
#   - STEP 1/2 semaglutide treatment effects
#   - Colocalization PP.H4
# =============================================================================

.libPaths(c("/path/to/project/R_libs", .libPaths()))
library(dplyr)

setwd("/path/to/project/regenie_pipeline")
outdir <- "results/twosampleMR/supplemental_tables"

cat("=== Building Master Integrated Table ===\n\n")

# ---- Protein mapping (protein_x.N -> gene symbol) ----
step1 <- read.csv("/path/to/project/STEP1_merged_results.csv")
prot_map <- step1 %>%
  dplyr::select(Exposure, code) %>%
  distinct() %>%
  rename(Protein = Exposure, Gene = code)

cat("Protein mapping:", nrow(prot_map), "entries\n")

# ---- Phenotypes to process ----
phenotypes <- c("BMI", "HbA1c", "TRIG_HDL_RATIO")
phenotype_labels <- c("BMI", "HbA1c", "TG/HDL-C Ratio")

all_tables <- list()

for (pi in seq_along(phenotypes)) {
  pheno <- phenotypes[pi]
  pheno_label <- phenotype_labels[pi]
  cat(sprintf("\n--- Processing %s ---\n", pheno))

  # ---- 1. Forward observational PEWAS ----
  fwd_obs_file <- sprintf("/path/to/project/UKB/PEWAS_results/%s_Linear_regression_proteomic_lm_results_non_T2D_adj_fasting_time.RDS", pheno)
  if (file.exists(fwd_obs_file)) {
    fwd_obs <- readRDS(fwd_obs_file)
    fwd_obs <- fwd_obs %>%
      filter(grepl("scale\\(protein", term)) %>%
      dplyr::select(Protein = Exposure, Fwd_Obs_Beta = estimate, Fwd_Obs_SE = std.error, Fwd_Obs_P = p.value) %>%
      distinct(Protein, .keep_all = TRUE)
    cat(sprintf("  Forward obs: %d proteins\n", nrow(fwd_obs)))
  } else {
    fwd_obs <- data.frame(Protein = character(), Fwd_Obs_Beta = numeric(),
                          Fwd_Obs_SE = numeric(), Fwd_Obs_P = numeric())
    cat("  Forward obs: file not found\n")
  }

  # ---- 2. Reverse observational PEWAS ----
  rev_obs_file <- sprintf("/path/to/project/UKB/PEWAS_results/Reverse_%s_Linear_regression_proteomic_lm_results_non_T2D_adj_fasting_time.RDS", pheno)
  if (file.exists(rev_obs_file)) {
    rev_obs <- readRDS(rev_obs_file)
    rev_obs <- rev_obs %>%
      dplyr::select(Protein, Rev_Obs_Beta = estimate, Rev_Obs_SE = std.error, Rev_Obs_P = p.value) %>%
      distinct(Protein, .keep_all = TRUE)
    cat(sprintf("  Reverse obs: %d proteins\n", nrow(rev_obs)))
  } else {
    rev_obs <- data.frame(Protein = character(), Rev_Obs_Beta = numeric(),
                          Rev_Obs_SE = numeric(), Rev_Obs_P = numeric())
    cat("  Reverse obs: file not found\n")
  }

  # ---- 3. Bidirectional MR ----
  mr_file <- sprintf("results/twosampleMR/supplemental_tables/Table_Bidirectional_MR_%s_Full.csv",
                      ifelse(pheno == "TRIG_HDL_RATIO", "TRIG_HDL_RATIO", pheno))
  if (file.exists(mr_file)) {
    mr <- read.csv(mr_file, check.names = FALSE)
    mr_clean <- data.frame(
      Gene = mr$Protein,
      Fwd_MR_Beta = as.numeric(mr$`Fwd: Beta`),
      Fwd_MR_SE = as.numeric(mr$`Fwd: SE`),
      Fwd_MR_P = as.numeric(mr$`Fwd: P-value`),
      Fwd_MR_Bonf = mr$`Fwd: Bonferroni Sig`,
      Rev_MR_Beta = as.numeric(mr$`Rev: IVW Beta`),
      Rev_MR_SE = as.numeric(mr$`Rev: IVW SE`),
      Rev_MR_P = as.numeric(mr$`Rev: IVW P-value`),
      Rev_MR_Bonf = mr$`Rev: Bonferroni Sig`,
      MR_Concordant = mr$`Directional Concordance`,
      stringsAsFactors = FALSE
    )
    cat(sprintf("  Bidirectional MR: %d proteins\n", nrow(mr_clean)))
  } else {
    mr_clean <- data.frame(Gene = character(), stringsAsFactors = FALSE)
    cat("  Bidirectional MR: file not found\n")
  }

  # ---- 4. STEP 1/2 treatment effects ----
  step1_pheno <- step1 %>%
    filter(Phenotype == pheno, Subgroup == "non T2D") %>%
    dplyr::select(code, Estimate_STEP1 = effect_size, Pval_STEP1 = pvalue) %>%
    group_by(code) %>% arrange(Pval_STEP1) %>% slice_head(n = 1) %>% ungroup() %>%
    rename(Gene = code)

  step2 <- read.csv("/path/to/project/STEP2_merged_results.csv")
  step2_pheno <- step2 %>%
    filter(Phenotype == pheno, Subgroup == "non T2D") %>%
    dplyr::select(code, Estimate_STEP2 = effect_size, Pval_STEP2 = pvalue) %>%
    group_by(code) %>% arrange(Pval_STEP2) %>% slice_head(n = 1) %>% ungroup() %>%
    rename(Gene = code)

  step_merged <- full_join(step1_pheno, step2_pheno, by = "Gene")
  step_merged$STEP_Mean <- rowMeans(step_merged[, c("Estimate_STEP1", "Estimate_STEP2")], na.rm = TRUE)
  cat(sprintf("  STEP treatment effects: %d genes\n", nrow(step_merged)))

  # ---- 5. Colocalization ----
  if (pheno == "BMI") {
    coloc <- read.csv("merged_step1_2_ukb_obesity_coloc.csv")
    coloc_sub <- coloc %>%
      filter(Phenotype == "BMI", Subgroup == "non T2D") %>%
      dplyr::select(code_UKB, PP.H4.abf, coloc_sig) %>%
      group_by(code_UKB) %>% arrange(desc(PP.H4.abf)) %>% slice_head(n = 1) %>% ungroup() %>%
      rename(Gene = code_UKB, Coloc_PP_H4 = PP.H4.abf, Coloc_Sig = coloc_sig)
  } else if (pheno == "HbA1c") {
    coloc <- read.csv("merged_step1_2_ukb_T2D_coloc.csv")
    coloc_sub <- coloc %>%
      filter(Phenotype == "BMI", Subgroup == "non T2D") %>%
      dplyr::select(code_UKB, PP.H4.abf, T2D_coloc_sig) %>%
      group_by(code_UKB) %>% arrange(desc(PP.H4.abf)) %>% slice_head(n = 1) %>% ungroup() %>%
      rename(Gene = code_UKB, Coloc_PP_H4 = PP.H4.abf, Coloc_Sig = T2D_coloc_sig)
  } else {
    coloc_sub <- data.frame(Gene = character(), Coloc_PP_H4 = numeric(), Coloc_Sig = integer())
  }
  cat(sprintf("  Colocalization: %d genes\n", nrow(coloc_sub)))

  # ---- Merge everything ----
  # Start with protein map
  master <- prot_map

  # Add forward obs
  master <- master %>% left_join(fwd_obs, by = "Protein")

  # Add reverse obs
  master <- master %>% left_join(rev_obs, by = "Protein")

  # Add bidirectional MR
  master <- master %>% left_join(mr_clean, by = "Gene")

  # Add STEP
  master <- master %>% left_join(step_merged, by = "Gene")

  # Add coloc
  if (nrow(coloc_sub) > 0) {
    master <- master %>% left_join(coloc_sub, by = "Gene")
  } else {
    master$Coloc_PP_H4 <- NA
    master$Coloc_Sig <- NA
  }

  # Deduplicate: keep one row per gene (best forward obs p-value)
  master <- master %>%
    group_by(Gene) %>%
    arrange(Fwd_Obs_P) %>%
    slice_head(n = 1) %>%
    ungroup()

  # Classify causal direction
  master <- master %>% mutate(
    Fwd_MR_Sig = !is.na(Fwd_MR_Bonf) & Fwd_MR_Bonf == "Yes",
    Rev_MR_Sig = !is.na(Rev_MR_Bonf) & Rev_MR_Bonf == "Yes",
    Causal_Class = case_when(
      Fwd_MR_Sig & Rev_MR_Sig & MR_Concordant == TRUE ~ "Bidirectional_Concordant",
      Fwd_MR_Sig & Rev_MR_Sig & MR_Concordant == FALSE ~ "Bidirectional_Discordant",
      Fwd_MR_Sig & !Rev_MR_Sig ~ "Forward_Only",
      !Fwd_MR_Sig & Rev_MR_Sig ~ "Reverse_Only",
      TRUE ~ "Neither"
    )
  )

  master$Phenotype <- pheno_label

  cat(sprintf("  Master table: %d rows\n", nrow(master)))
  cat(sprintf("    Forward-only: %d | Reverse-only: %d | Bidirectional concordant: %d | Bidirectional discordant: %d\n",
              sum(master$Causal_Class == "Forward_Only", na.rm = TRUE),
              sum(master$Causal_Class == "Reverse_Only", na.rm = TRUE),
              sum(master$Causal_Class == "Bidirectional_Concordant", na.rm = TRUE),
              sum(master$Causal_Class == "Bidirectional_Discordant", na.rm = TRUE)))

  all_tables[[pheno]] <- master
}

# ---- Combine and save ----
final <- bind_rows(all_tables)

outfile <- file.path(outdir, "Master_Integrated_Table.csv")
write.csv(final, outfile, row.names = FALSE)
cat(sprintf("\n=== Saved: %s (%d rows) ===\n", outfile, nrow(final)))

# ---- Summary across phenotypes ----
cat("\n=== SUMMARY ===\n")
final %>%
  group_by(Phenotype, Causal_Class) %>%
  summarise(N = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Causal_Class, values_from = N, values_fill = 0) %>%
  print()

cat("\nDone.\n")
