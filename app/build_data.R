# =============================================================================
# build_data.R  ---  MAP-D Shiny app data builder
#
# Converts SUPP_TABLE_1_05_20_2026.xlsx (the manuscript supplementary workbook,
# the single source of truth) into a single fast-loading named list:  data/mapd.rds
#
# Run once at build time:   Rscript build_data.R
# The app (app.R) only ever reads data/mapd.rds -- it never touches the xlsx.
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
})

XLSX <- "SUPP_TABLE_1_05_20_2026.xlsx"
OUT  <- "data/mapd.rds"
if (!dir.exists("data")) dir.create("data")

message("Reading workbook: ", XLSX)

# ---- canonical hallmark naming ---------------------------------------------
# Master sheet uses "TG/HDL-C Ratio"; STEP / GTEx use "TRIG_HDL_RATIO".
# We canonicalise everything to: BMI, HbA1c, TG/HDL
canon_hallmark <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    x %in% c("TG/HDL-C Ratio", "TRIG_HDL_RATIO", "TG/HDL ratio", "TG_HDL_ratio") ~ "TG/HDL",
    TRUE ~ x
  )
}
HALLMARKS <- c("BMI", "HbA1c", "TG/HDL")
STAGES    <- c("Normoglycemic", "Prediabetes", "T2D")

# split "GENE;Full protein name" into gene + description
split_gene <- function(x) sub(";.*$", "", x)
split_desc <- function(x) ifelse(grepl(";", x), sub("^[^;]*;", "", x), x)

# =============================================================================
# 1. MASTER integrated evidence  (Fwd_Rev_Integrated_Evidence)
#    2923 proteins x 3 hallmarks = 8769 rows. Drives triangulation, MR, cards.
# =============================================================================
master <- read_excel(XLSX, sheet = "Fwd_Rev_Integrated_Evidence") %>%
  # some numeric cells are stored as text -> force numeric on all measure columns
  mutate(across(c(Fwd_Obs_Beta, Fwd_Obs_SE, Fwd_Obs_P, Rev_Obs_Beta, Rev_Obs_SE,
                  Rev_Obs_P, Fwd_MR_Beta, Fwd_MR_SE, Fwd_MR_P, Rev_MR_Beta, Rev_MR_SE,
                  Rev_MR_P, Estimate_STEP1, Pval_STEP1, Estimate_STEP2, Pval_STEP2,
                  STEP_Mean, Coloc_PP_H4),
                ~ suppressWarnings(as.numeric(.)))) %>%
  mutate(
    Hallmark   = canon_hallmark(Phenotype),
    Fwd_MR_Sig = as.logical(Fwd_MR_Sig),
    Rev_MR_Sig = as.logical(Rev_MR_Sig),
    # STEP effect used for the trial comparison at the protein level.
    # STEP1 = non-T2D trial (BMI-aligned), STEP2 = T2D trial (HbA1c-aligned).
    # UKB direction here = reverse observational beta (hallmark -> protein).
    UKB_dir  = sign(Rev_Obs_Beta),
    # Quadrant vs STEP mean (used on the MR / integrated views)
    STEP_dir = sign(STEP_Mean)
  )

# reversion / persistence quadrant, per Fig 4 legend, given a UKB beta & STEP beta
quadrant_label <- function(ukb, step) {
  dplyr::case_when(
    is.na(ukb) | is.na(step)      ~ NA_character_,
    ukb > 0 & step < 0 ~ "Reversion (UKB+ / GLP1RA−)",
    ukb < 0 & step > 0 ~ "Reversion (UKB− / GLP1RA+)",
    ukb > 0 & step > 0 ~ "Persistent (UKB+ / GLP1RA+)",
    ukb < 0 & step < 0 ~ "Persistent (UKB− / GLP1RA−)",
    TRUE ~ NA_character_
  )
}
quadrant_class <- function(ukb, step) {
  dplyr::case_when(
    is.na(ukb) | is.na(step)             ~ NA_character_,
    sign(ukb) != sign(step)              ~ "Reversion",
    TRUE                                 ~ "Persistent"
  )
}

# ---- protein full-name -> gene-symbol bridge --------------------------------
# STEP1_2_VALIDATION stores full protein names (no gene symbol). The GTEx and CAD
# sheets carry both, so we build a name->gene map to key STEP rows by gene symbol
# (needed to join MR / master evidence and to power the Protein Explorer).
gtex_map <- read_excel(XLSX, sheet = "MAP-D_GTEX") %>%
  transmute(name = UKB_Protein, Gene = Gene) %>% distinct()
cad_map <- read_excel(XLSX, sheet = "INCIDENT CAD REGRESSION OUTPUT") %>%
  transmute(name = split_desc(Protein), Gene = split_gene(Protein)) %>% distinct()
name2gene <- bind_rows(gtex_map, cad_map) %>%
  filter(!is.na(name), !is.na(Gene)) %>%
  distinct(name, .keep_all = TRUE)
map_gene <- function(nm) name2gene$Gene[match(nm, name2gene$name)]

# =============================================================================
# 2. STAGE-STRATIFIED validation  (STEP1_2_VALIDATION)
#    UKB reverse-obs beta per protein x hallmark x glycemic stage vs STEP betas.
# =============================================================================
step_stage <- read_excel(XLSX, sheet = "STEP1_2_VALIDATION") %>%
  mutate(across(c(Estimate_UKB, Significance_UKB, Estimate_STEP1, Significance_STEP1,
                  Estimate_STEP2, Significance_STEP2),
                ~ suppressWarnings(as.numeric(.)))) %>%
  mutate(Hallmark = canon_hallmark(Phenotype)) %>%
  filter(Hallmark %in% HALLMARKS) %>%
  mutate(
    Subgroup = factor(Subgroup, levels = STAGES),
    Gene     = map_gene(Protein_UKB),
    # pick the stage-appropriate trial: non-T2D stages -> STEP1, T2D -> STEP2
    STEP_used   = ifelse(Subgroup == "T2D", "STEP2", "STEP1"),
    STEP_beta   = ifelse(Subgroup == "T2D", Estimate_STEP2, Estimate_STEP1),
    STEP_pval   = ifelse(Subgroup == "T2D", Significance_STEP2, Significance_STEP1),
    Quadrant       = quadrant_label(Estimate_UKB, STEP_beta),
    Quadrant_class = quadrant_class(Estimate_UKB, STEP_beta)
  )

# =============================================================================
# 3. Evidence funnel  (Evidence_Funnel_Table)  -- long tidy form for plotting
# =============================================================================
funnel_raw <- read_excel(XLSX, sheet = "Evidence_Funnel_Table")
funnel <- funnel_raw %>%
  pivot_longer(-Funnel_stage, names_to = "col", values_to = "value") %>%
  separate(col, into = c("Hallmark", "Subgroup"), sep = " \\| ") %>%
  mutate(
    Hallmark = canon_hallmark(Hallmark),
    Subgroup = factor(Subgroup, levels = STAGES),
    stage_num = as.integer(str_extract(Funnel_stage, "^[0-9]+")),
    stage_lab = str_trim(str_remove(Funnel_stage, "^[0-9]+\\.\\s*")),
    # numeric count = leading integer of the (possibly annotated) value
    n = suppressWarnings(as.numeric(str_extract(as.character(value), "^[0-9]+")))
  ) %>%
  arrange(Hallmark, Subgroup, stage_num)

# short human labels for the 5 funnel stages
funnel_stage_short <- c(
  "1" = "UKB observational (Bonferroni)",
  "2" = "+ STEP trial (P < 0.05)",
  "3" = "+ reverse MR (IVW P < 0.05)",
  "4" = "+ MR concordant",
  "5" = "CAD DrugBank hits"
)
funnel$stage_short <- funnel_stage_short[as.character(funnel$stage_num)]

# =============================================================================
# 4. Incident CAD  (INCIDENT CAD REGRESSION OUTPUT)
# =============================================================================
cad <- read_excel(XLSX, sheet = "INCIDENT CAD REGRESSION OUTPUT") %>%
  mutate(across(c(estimate, std.error, statistic, p.value, FDR),
                ~ suppressWarnings(as.numeric(.)))) %>%
  mutate(
    Gene   = split_gene(Protein),
    Desc   = split_desc(Protein),
    OR     = exp(estimate),
    OR_lo  = exp(estimate - 1.96 * std.error),
    OR_hi  = exp(estimate + 1.96 * std.error),
    neglog10FDR = -log10(FDR),
    Sig    = FDR < 0.05,
    Direction = ifelse(estimate > 0, "Higher risk (OR>1)", "Lower risk (OR<1)")
  ) %>%
  arrange(FDR)

# =============================================================================
# 5. DrugBank hits  (DRUGBANK_HITS)
# =============================================================================
drugbank <- read_excel(XLSX, sheet = "DRUGBANK_HITS") %>%
  mutate(across(c(estimate, std.error, p.value, FDR, OR),
                ~ suppressWarnings(as.numeric(.)))) %>%
  mutate(
    OR = ifelse(!is.na(OR), OR, exp(estimate)),
    Direction = ifelse(estimate > 0, "Higher risk (OR>1)", "Lower risk (OR<1)")
  )

# =============================================================================
# 6. GTEx tissue enrichment  (MAP-D_GTEX)
# =============================================================================
gtex <- read_excel(XLSX, sheet = "MAP-D_GTEX") %>%
  mutate(Hallmark = canon_hallmark(Phenotype)) %>%
  filter(Hallmark %in% HALLMARKS) %>%
  mutate(Subgroup = factor(Subgroup, levels = STAGES))

# =============================================================================
# 7. LDSC genetic correlations  (LDSC_genetic_corrs)
# =============================================================================
ldsc <- read_excel(XLSX, sheet = "LDSC_genetic_corrs")

# =============================================================================
# 8. Analysis-layer model specs  (Analysis_Layer_Model_Specs)
# =============================================================================
specs <- read_excel(XLSX, sheet = "Analysis_Layer_Model_Specs")
names(specs) <- str_trim(names(specs))
# Reverse Observational PWAS did not adjust for smoking status (per manuscript methods)
specs$Covariates <- ifelse(
  specs[["Analysis Layer"]] == "Reverse Observational PWAS",
  str_squish(str_remove(specs$Covariates, regex(",?\\s*smoking status", ignore_case = TRUE))),
  specs$Covariates
)

# =============================================================================
# Derived summaries
# =============================================================================
# Causal-class counts per hallmark (BMI reverse-dominant, HbA1c forward, TG/HDL mixed)
causal_summary <- master %>%
  count(Hallmark, Causal_Class) %>%
  pivot_wider(names_from = Causal_Class, values_from = n, values_fill = 0)

mr_summary <- master %>%
  group_by(Hallmark) %>%
  summarise(
    Forward_MR_sig = sum(Fwd_MR_Sig, na.rm = TRUE),
    Reverse_MR_sig = sum(Rev_MR_Sig, na.rm = TRUE),
    Bidirectional  = sum(Fwd_MR_Sig & Rev_MR_Sig, na.rm = TRUE),
    .groups = "drop"
  )

# Phenotypic hallmark correlation of reverse-obs betas, per stage (Fig 3A analogue)
hallmark_corr <- step_stage %>%
  select(Protein_UKB, Hallmark, Subgroup, Estimate_UKB) %>%
  distinct(Protein_UKB, Hallmark, Subgroup, .keep_all = TRUE) %>%
  group_split(Subgroup) %>%
  map(function(d) {
    stage <- as.character(d$Subgroup[1])
    w <- d %>% select(Protein_UKB, Hallmark, Estimate_UKB) %>%
      pivot_wider(names_from = Hallmark, values_from = Estimate_UKB)
    hcols <- intersect(HALLMARKS, names(w))
    m <- suppressWarnings(cor(w[hcols], use = "pairwise.complete.obs"))
    as.data.frame(as.table(m)) %>%
      setNames(c("Hallmark_A", "Hallmark_B", "r")) %>%
      mutate(Subgroup = stage)
  }) %>%
  bind_rows()

# Proteome-wide cross-hallmark correlation of hallmark->protein (reverse-obs) betas,
# across all 2,923 proteins (master sheet). Headline for the Atlas architecture view.
hallmark_corr_all <- {
  w <- master %>%
    select(Gene, Hallmark, Rev_Obs_Beta) %>%
    pivot_wider(names_from = Hallmark, values_from = Rev_Obs_Beta)
  hcols <- intersect(HALLMARKS, names(w))
  m <- suppressWarnings(cor(w[hcols], use = "pairwise.complete.obs"))
  as.data.frame(as.table(m)) %>%
    setNames(c("Hallmark_A", "Hallmark_B", "r"))
}

# Protein <-> gene search index (union across every layer)
search_index <- sort(unique(c(master$Gene, step_stage$Gene, cad$Gene)))

meta <- list(
  n_proteins   = dplyr::n_distinct(master$Gene),
  n_hallmarks  = length(HALLMARKS),
  n_stages     = length(STAGES),
  hallmarks    = HALLMARKS,
  stages       = STAGES,
  n_cad_tested = nrow(cad),
  n_cad_sig    = sum(cad$Sig, na.rm = TRUE),
  n_drug_pairs = nrow(drugbank),
  n_drugs      = dplyr::n_distinct(drugbank$drug_name),
  n_drug_prot  = dplyr::n_distinct(drugbank$gene_symbol),
  built_on     = "2026-06-30",
  source_file  = XLSX
)

mapd <- list(
  master        = master,
  step_stage    = step_stage,
  funnel        = funnel,
  cad           = cad,
  drugbank      = drugbank,
  gtex          = gtex,
  ldsc          = ldsc,
  specs         = specs,
  causal_summary = causal_summary,
  mr_summary    = mr_summary,
  hallmark_corr = hallmark_corr,
  hallmark_corr_all = hallmark_corr_all,
  search_index  = search_index,
  meta          = meta
)

saveRDS(mapd, OUT)
message("Wrote ", OUT, " (", round(file.size(OUT)/1e6, 2), " MB)")

# =============================================================================
# Sanity checks -- assert the data reproduces the manuscript's headline numbers
# =============================================================================
message("\n---- SANITY CHECKS ----")
message("Proteins: ", meta$n_proteins, " (expect 2923)")
stopifnot(meta$n_proteins == 2923)

chk <- function(hallmark, fwd, rev) {
  row <- mr_summary[mr_summary$Hallmark == hallmark, ]
  message(sprintf("%-8s  forward-MR=%d (expect %d)  reverse-MR=%d (expect %d)",
                  hallmark, row$Forward_MR_sig, fwd, row$Reverse_MR_sig, rev))
  stopifnot(row$Forward_MR_sig == fwd, row$Reverse_MR_sig == rev)
}
chk("BMI",    10, 106)   # adiposity: reverse causation dominates
chk("HbA1c",  19, 0)     # glycemia: protein-first only
chk("TG/HDL", 17, 34)    # insulin resistance: bidirectional

message("CAD: ", meta$n_cad_sig, " FDR<0.05 of ", meta$n_cad_tested, " (expect 367 of 1008)")
stopifnot(meta$n_cad_tested == 1008, meta$n_cad_sig == 367)

message("DrugBank: ", meta$n_drug_pairs, " pairs, ", meta$n_drugs, " drugs, ",
        meta$n_drug_prot, " proteins")

message("\nAll sanity checks passed.")
