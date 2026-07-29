## Variance decomposition: BMI ~ circulating proteins (Olink Explore 3072)
## Stratified by Olink-SomaScan overlap ("triangulable") vs Olink-only
## Nature Medicine Perspective — Figure 2
##
## SomaScan reference: `inputs/somascan_v41_7k_human.csv` (columns uniprot, gene_symbol)
## if present; otherwise built from STEP1/STEP2 merged analyte lists (STEP-context proxy).
##
## Runtime / QC knobs (see variance_decomposition_methodology_notes.txt):
##   VARDECOMP_MAX_N          — max participants (default 30000)
##   VARDECOMP_MAX_PROT_MISS  — max missingness per protein, 0–1 (default 0.30); auto-relaxes if needed
##   VARDECOMP_NFOLDS         — glmnet CV folds (default 3)

.libPaths(c("/path/to/project/R_libs", .libPaths()))

library(glmnet)
library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)

RUN_SEED <- 2026L
set.seed(RUN_SEED)

out_dir <- "/path/to/project/regenie_pipeline"
setwd(out_dir)

cv_folds <- suppressWarnings(as.integer(Sys.getenv("VARDECOMP_NFOLDS", "3")))
if (is.na(cv_folds) || cv_folds < 3L) cv_folds <- 3L

max_prot_miss_target <- suppressWarnings(as.numeric(Sys.getenv("VARDECOMP_MAX_PROT_MISS", "0.30")))
if (is.na(max_prot_miss_target) || max_prot_miss_target <= 0 || max_prot_miss_target > 1) {
  max_prot_miss_target <- 0.30
}

## ── 1. Load data ──────────────────────────────────────────────────────────────

cat("Loading data...\n")
data <- readRDS("/path/to/project/UKB/PEWAS_results/data_plus_GLP_complications_glycemic_status_HbA1c_adjusted.RDS")
cat("  Raw dimensions:", dim(data), "\n")

coding <- read.delim("/path/to/project/UKB/PEWAS_results/coding143.tsv",
                     stringsAsFactors = FALSE)
cat("  Coding entries:", nrow(coding), "\n")

## ── 2. Build Olink–SomaScan overlap mapping ──────────────────────────────────

cat("Building protein overlap mapping...\n")

# Olink annotation from Eldjarn/Sun et al. 2023 supplementary ST3
st3 <- read_excel(file.path(out_dir, "41586_2023_6592_MOESM3_ESM.xlsx"),
                  sheet = "ST3", skip = 1)

# Map coding index to ST3 by gene symbol (Assay Target matches coding gene)
coding$gene_symbol <- sub(";.*", "", coding$meaning)
coding$protein_name <- sub("^[^;]+;", "", coding$meaning)

olink <- data.frame(
  protein_col = paste0("protein_x.", coding$coding),
  coding_idx  = coding$coding,
  gene_symbol = coding$gene_symbol,
  protein_name = coding$protein_name,
  stringsAsFactors = FALSE
)

# Merge with ST3 to get UniProt
st3_slim <- st3 %>%
  select(gene_symbol = `Gene symbol`, uniprot_id = UniProt, olink_id = `Olink ID`) %>%
  distinct(gene_symbol, .keep_all = TRUE)

olink <- olink %>% left_join(st3_slim, by = "gene_symbol")

# SomaScan analyte universe for Olink overlap (7k manifest if available; else STEP trial union)
load_or_build_somascan <- function(out_dir) {
  soma_path_7k <- file.path(out_dir, "inputs", "somascan_v41_7k_human.csv")
  step1 <- "/path/to/project/STEP1_merged_results.csv"
  step2 <- "/path/to/project/STEP2_merged_results.csv"
  if (file.exists(soma_path_7k)) {
    cat("  Using SomaScan manifest:", soma_path_7k, "\n")
    return(list(
      soma = read.csv(soma_path_7k, stringsAsFactors = FALSE),
      source = "SomaScan v4.1 7k manifest (inputs/somascan_v41_7k_human.csv)"
    ))
  }
  cat("  No inputs/somascan_v41_7k_human.csv — building analyte union from STEP merged tables.\n")
  if (!file.exists(step1) || !file.exists(step2)) {
    stop(
      "Need either inputs/somascan_v41_7k_human.csv or both STEP1_merged_results.csv and ",
      "STEP2_merged_results.csv under /path/to/project/"
    )
  }
  s1 <- read.csv(step1, stringsAsFactors = FALSE)
  s2 <- read.csv(step2, stringsAsFactors = FALSE)
  step <- rbind(s1, s2)
  up <- unique(unlist(strsplit(na.omit(as.character(step$UniProt)), "\\|")))
  up <- up[!is.na(up) & nzchar(up)]
  genes <- unique(c(
    na.omit(as.character(step$code)),
    na.omit(as.character(step$Target))
  ))
  genes <- genes[nzchar(genes) & !is.na(genes)]
  soma <- rbind(
    data.frame(uniprot = up, gene_symbol = "", stringsAsFactors = FALSE),
    data.frame(uniprot = "", gene_symbol = genes, stringsAsFactors = FALSE)
  )
  dir.create(file.path(out_dir, "inputs"), showWarnings = FALSE, recursive = TRUE)
  out_csv <- file.path(out_dir, "inputs", "somascan_step_union.csv")
  write.csv(soma, out_csv, row.names = FALSE)
  cat("  Wrote", nrow(soma), "reference rows to", out_csv, "\n")
  list(
    soma = soma,
    source = "STEP1/STEP2 merged analyte union (UniProt + gene; STEP-context proxy for SomaScan overlap)"
  )
}

soma_info <- load_or_build_somascan(out_dir)
soma <- soma_info$soma
use_gene_fallback <- !grepl("STEP1/STEP2", soma_info$source, fixed = TRUE)

# Expand SomaScan multi-UniProt entries
soma_uniprots <- unique(unlist(strsplit(na.omit(soma$uniprot), "\\|")))
soma_uniprots <- soma_uniprots[soma_uniprots != ""]
soma_genes    <- unique(unlist(strsplit(na.omit(soma$gene_symbol), "\\|")))
soma_genes    <- soma_genes[soma_genes != ""]

cat("  SomaScan unique UniProts:", length(soma_uniprots), "\n")
cat("  SomaScan unique genes:", length(soma_genes), "\n")
if (!use_gene_fallback) {
  cat("  (STEP-derived list: triangulation uses UniProt match only; gene fallback disabled.)\n")
}

# Match: UniProt first, then gene symbol fallback (7k manifest only; STEP union is UniProt-only)
olink$match_uniprot <- olink$uniprot_id %in% soma_uniprots
olink$match_gene    <- if (use_gene_fallback) {
  olink$gene_symbol %in% soma_genes
} else {
  rep(FALSE, nrow(olink))
}
olink$on_somascan_7k <- olink$match_uniprot | olink$match_gene
olink$somascan_match_method <- ifelse(
  olink$match_uniprot, "uniprot",
  ifelse(olink$match_gene, "gene_symbol", "none")
)

cat("  Olink proteins on SomaScan (triangulable):", sum(olink$on_somascan_7k), "\n")
cat("  Olink-only proteins:", sum(!olink$on_somascan_7k), "\n")

# Save mapping
mapping_out <- olink %>%
  select(olink_id, uniprot_id, gene_symbol, protein_name, on_somascan_7k, somascan_match_method)
write.csv(mapping_out, file.path(out_dir, "protein_overlap_mapping.csv"), row.names = FALSE)
cat("  Saved protein_overlap_mapping.csv\n")

## ── 3. Prepare analytic sample ───────────────────────────────────────────────

cat("Preparing analytic sample...\n")

max_n <- suppressWarnings(as.integer(Sys.getenv("VARDECOMP_MAX_N", unset = "30000")))
if (is.na(max_n) || max_n <= 0) max_n <- 30000L

# Sample participants *before* materialising the full European subset (avoids OOM)
eu_codes <- c("British", "Irish", "Any other white background")
eu_mask <- data$f.21000.0.0 %in% eu_codes
eid_eu <- unique(data$eid[eu_mask & !is.na(data$eid)])
cat("  European ancestry unique eids (full cohort):", length(eid_eu), "\n")

if (length(eid_eu) > max_n) {
  set.seed(RUN_SEED)
  eid_eu <- sample(eid_eu, max_n)
  cat("  Random sample of participants:", max_n, "(env VARDECOMP_MAX_N)\n")
}

keep <- eu_mask & data$eid %in% eid_eu
data_eur <- data[keep, , drop = FALSE]
rm(data)
gc(FALSE)
cat("  Rows after ancestry + participant filter:", nrow(data_eur), "\n")

data_eur <- data_eur %>% distinct(eid, .keep_all = TRUE)
cat("  After dedup:", nrow(data_eur), "\n")

# Protein columns (excluding ins_index)
protein_vars <- paste0("protein_x.", 1:2923)
protein_vars <- protein_vars[protein_vars %in% colnames(data_eur)]

# Covariate columns
covar_cols <- c("x.age", "x.sex", "f.54.0.0",
                paste0("f.22009.0.", 1:10))  # first 10 genetic PCs

# Check BMI availability
data_eur <- data_eur %>% filter(!is.na(BMI))
cat("  With non-missing BMI:", nrow(data_eur), "\n")

# Check covariate availability
for (cc in covar_cols) {
  if (!cc %in% colnames(data_eur)) {
    cat("  WARNING: covariate", cc, "not found\n")
  }
}
data_eur <- data_eur %>% drop_na(all_of(covar_cols))
cat("  After dropping missing covariates:", nrow(data_eur), "\n")

# Convert assessment centre to factor
data_eur$f.54.0.0 <- as.factor(data_eur$f.54.0.0)

## ── 4. Protein QC: missingness filter + imputation + rank normalisation ──────

cat("Protein QC...\n")

# Coerce protein columns to numeric (avoids silent all-NA behaviour if types differ)
for (pv in protein_vars) {
  if (!is.numeric(data_eur[[pv]])) {
    data_eur[[pv]] <- suppressWarnings(as.numeric(as.character(data_eur[[pv]])))
  }
}

# Missingness per protein (NA only)
prot_miss <- sapply(protein_vars, function(p) mean(is.na(data_eur[[p]])))

# Threshold: default 30% max missingness per analyte (common in proteomics QC); escalate if empty
relax_seq <- unique(
  c(max_prot_miss_target, seq(0.40, 0.95, by = 0.05))
)
relax_seq <- sort(relax_seq[relax_seq > 0 & relax_seq <= 1])

keep_prots <- character(0)
missingness_threshold_used <- NA_real_
for (thr in relax_seq) {
  kp <- names(prot_miss)[prot_miss <= thr]
  if (length(kp) > 0) {
    keep_prots <- kp
    missingness_threshold_used <- thr
    if (abs(thr - max_prot_miss_target) > 1e-6) {
      cat("  NOTE: auto-relaxed per-protein missingness ceiling to ", round(100 * thr, 1),
          "% (target was ", round(100 * max_prot_miss_target, 1), "%)\n", sep = "")
    }
    break
  }
}

if (length(keep_prots) == 0) {
  stop(
    "No proteins passed missingness filter even at 95%. Check protein columns in the RDS ",
    "(types, coding of missing values)."
  )
}

cat(
  "  Per-protein missingness ceiling: ", round(100 * missingness_threshold_used, 2),
  "% — retained ", length(keep_prots), " / ", length(protein_vars), " proteins\n",
  sep = ""
)

# Median imputation for remaining missingness
for (p in keep_prots) {
  nas <- is.na(data_eur[[p]])
  if (any(nas)) {
    data_eur[[p]][nas] <- median(data_eur[[p]], na.rm = TRUE)
  }
}

# Inverse-rank-normal transformation
rank_normal <- function(x) {
  r <- rank(x, na.last = "keep", ties.method = "average")
  qnorm((r - 0.5) / sum(!is.na(x)))
}

cat("  Inverse-rank normalising proteins...\n")
for (p in keep_prots) {
  data_eur[[p]] <- rank_normal(data_eur[[p]])
}

# Define protein sets
prot_triangulable <- olink %>%
  filter(on_somascan_7k, protein_col %in% keep_prots) %>%
  pull(protein_col)
prot_olink_only <- olink %>%
  filter(!on_somascan_7k, protein_col %in% keep_prots) %>%
  pull(protein_col)

cat("  Triangulable proteins (after QC):", length(prot_triangulable), "\n")
cat("  Olink-only proteins (after QC):", length(prot_olink_only), "\n")
cat("  Total proteins (after QC):", length(keep_prots), "\n")

## ── 5. Train/test split (80/20 stratified by sex) ───────────────────────────

cat("Train/test split...\n")

data_eur$split_group <- data_eur$x.sex
train_idx <- c()
for (g in unique(data_eur$split_group)) {
  idx <- which(data_eur$split_group == g)
  n_train <- round(0.8 * length(idx))
  train_idx <- c(train_idx, sample(idx, n_train))
}
test_idx <- setdiff(seq_len(nrow(data_eur)), train_idx)

cat("  Train N:", length(train_idx), "  Test N:", length(test_idx), "\n")

## ── 6. Build covariate matrix (assessment centre → dummies) ──────────────────

build_X <- function(df, covar_cols, prot_cols = NULL) {
  # Covariates: numeric + dummy-coded factor
  covar_parts <- list()
  for (cc in covar_cols) {
    if (is.factor(df[[cc]]) || is.character(df[[cc]])) {
      dummies <- model.matrix(~ df[[cc]] - 1)
      # drop first level to avoid collinearity
      if (ncol(dummies) > 1) dummies <- dummies[, -1, drop = FALSE]
      colnames(dummies) <- paste0(cc, "_", seq_len(ncol(dummies)))
      covar_parts[[cc]] <- dummies
    } else {
      covar_parts[[cc]] <- matrix(df[[cc]], ncol = 1, dimnames = list(NULL, cc))
    }
  }
  X_covar <- do.call(cbind, covar_parts)

  if (is.null(prot_cols) || length(prot_cols) == 0) {
    return(X_covar)
  }

  X_prot <- as.matrix(df[, prot_cols])
  cbind(X_covar, X_prot)
}

n_covar <- ncol(build_X(data_eur[1:5, ], covar_cols))
cat("  Number of covariate columns (with dummies):", n_covar, "\n")

## ── 7. Elastic net fitting function ──────────────────────────────────────────

fit_enet <- function(X_train, y_train, X_test, y_test,
                     penalty_factor = NULL, nfolds = cv_folds) {
  if (is.null(penalty_factor)) {
    penalty_factor <- rep(1, ncol(X_train))
  }

  cv_fit <- cv.glmnet(
    X_train, y_train, alpha = 0.5,
    penalty.factor = penalty_factor,
    nfolds = nfolds, type.measure = "mse", standardize = TRUE
  )

  pred_test <- predict(cv_fit, newx = X_test, s = "lambda.min")
  ss_res <- sum((y_test - pred_test)^2)
  ss_tot <- sum((y_test - mean(y_test))^2)
  r2 <- 1 - ss_res / ss_tot

  n_nonzero <- sum(coef(cv_fit, s = "lambda.min")[-1] != 0)

  list(r2 = r2, lambda = cv_fit$lambda.min, n_nonzero = n_nonzero, cv_fit = cv_fit)
}

## ── 8. Run three models for BMI ──────────────────────────────────────────────

y_train <- data_eur$BMI[train_idx]
y_test  <- data_eur$BMI[test_idx]

cat("\n=== MODEL 0: Covariates only ===\n")
X0_train <- build_X(data_eur[train_idx, ], covar_cols)
X0_test  <- build_X(data_eur[test_idx, ], covar_cols)
m0 <- fit_enet(X0_train, y_train, X0_test, y_test, nfolds = cv_folds)
cat("  R² =", round(m0$r2, 4), "\n")

cat("\n=== MODEL 1: Covariates + Triangulable proteins ===\n")
if (length(prot_triangulable) == 0) {
  stop("Triangulable protein set is empty after QC; check Olink–SomaScan mapping.")
}
X1_train <- build_X(data_eur[train_idx, ], covar_cols, prot_triangulable)
X1_test  <- build_X(data_eur[test_idx, ], covar_cols, prot_triangulable)
pf1 <- c(rep(0, n_covar), rep(1, length(prot_triangulable)))
m1 <- fit_enet(X1_train, y_train, X1_test, y_test, penalty_factor = pf1, nfolds = cv_folds)
cat("  R² =", round(m1$r2, 4), "\n")

cat("\n=== MODEL 2: Covariates + All Olink proteins ===\n")
X2_train <- build_X(data_eur[train_idx, ], covar_cols, keep_prots)
X2_test  <- build_X(data_eur[test_idx, ], covar_cols, keep_prots)
pf2 <- c(rep(0, n_covar), rep(1, length(keep_prots)))
m2 <- fit_enet(X2_train, y_train, X2_test, y_test, penalty_factor = pf2, nfolds = cv_folds)
cat("  R² =", round(m2$r2, 4), "\n")

## ── 9. Variance decomposition ────────────────────────────────────────────────

r2_baseline      <- max(m0$r2, 0)
r2_triangulable  <- max(m1$r2 - m0$r2, 0)
r2_olink_only    <- max(m2$r2 - m1$r2, 0)
r2_residual      <- max(1 - m2$r2, 0)

results <- data.frame(
  trait = "BMI",
  model = c("Model 0", "Model 1", "Model 2"),
  protein_set = c("Covariates only",
                  "Olink ∩ SomaScan",
                  "Full Olink"),
  n_proteins = c(0, length(prot_triangulable), length(keep_prots)),
  r2 = c(m0$r2, m1$r2, m2$r2),
  incremental_r2_over_baseline = c(NA, m1$r2 - m0$r2, m2$r2 - m0$r2),
  incremental_r2_over_triangulable = c(NA, NA, m2$r2 - m1$r2),
  stringsAsFactors = FALSE
)

cat("\n=== Results ===\n")
print(results)

write.csv(results, file.path(out_dir, "variance_decomposition_results.csv"), row.names = FALSE)

## ── 9b. Review summary — key R² for manuscript / triangulation perspective ───

key_summary <- data.frame(
  model = c(
    "covariates_only",
    "olink_intersect_somascan_triangulable",
    "full_olink_panel"
  ),
  description = c(
    "Covariates only (age, sex, assessment centre, 10 PCs)",
    "Covariates + proteins on both Olink and SomaScan (triangulable set)",
    "Covariates + full post-QC Olink panel"
  ),
  n_proteins = c(0L, length(prot_triangulable), length(keep_prots)),
  r2_test_set = c(m0$r2, m1$r2, m2$r2),
  incremental_r2_over_covariates = c(0, m1$r2 - m0$r2, m2$r2 - m0$r2),
  stringsAsFactors = FALSE
)
key_summary$somascan_reference <- soma_info$source

write.csv(key_summary, file.path(out_dir, "bmi_proteome_r2_summary.csv"), row.names = FALSE)

run_metadata <- data.frame(
  setting = c(
    "random_seed",
    "participant_cap_requested",
    "participants_after_filters",
    "train_n",
    "test_n",
    "per_protein_missingness_target",
    "per_protein_missingness_used",
    "glmnet_cv_folds",
    "elastic_net_alpha",
    "covariates",
    "soma_reference",
    "triangulation_match_rule"
  ),
  value = c(
    RUN_SEED,
    max_n,
    nrow(data_eur),
    length(train_idx),
    length(test_idx),
    max_prot_miss_target,
    missingness_threshold_used,
    cv_folds,
    0.5,
    "age, sex, assessment centre (factor), first 10 genetic PCs",
    soma_info$source,
    paste0(
      "Official 7k manifest: UniProt + gene-symbol fallback; ",
      "STEP-derived list: UniProt only (gene fallback disabled—see methodology notes)"
    )
  ),
  stringsAsFactors = FALSE
)
write.csv(run_metadata, file.path(out_dir, "variance_decomposition_run_metadata.csv"), row.names = FALSE)
cat("  Wrote variance_decomposition_run_metadata.csv\n")

summary_txt <- c(
  "BMI — proteome variance explained (UK Biobank European ancestry; out-of-sample R², elastic net)",
  paste0("Triangulation / SomaScan reference: ", soma_info$source),
  paste0(
    "Analytic settings: N≤", max_n, " EUR participants; per-protein missingness ≤",
    round(100 * missingness_threshold_used, 2), "%; glmnet ", cv_folds, "-fold CV; seed ", RUN_SEED
  ),
  "",
  paste0("1) Olink ∩ SomaScan (triangulable):  R² = ", sprintf("%.4f", m1$r2)),
  paste0("     N proteins (QC’d): ", length(prot_triangulable)),
  paste0("     ΔR² vs covariates only: ", sprintf("%.4f", m1$r2 - m0$r2)),
  "",
  paste0("2) Full Olink panel:                 R² = ", sprintf("%.4f", m2$r2)),
  paste0("     N proteins (QC’d): ", length(keep_prots)),
  paste0("     ΔR² vs covariates only: ", sprintf("%.4f", m2$r2 - m0$r2)),
  "",
  paste0("   Covariates-only baseline:         R² = ", sprintf("%.4f", m0$r2)),
  paste0("   Train N = ", length(train_idx), "; test N = ", length(test_idx)),
  ""
)
writeLines(summary_txt, file.path(out_dir, "bmi_proteome_r2_summary.txt"))
cat("\n", paste(summary_txt, collapse = "\n"), "\n", sep = "")
cat("  Wrote bmi_proteome_r2_summary.txt and bmi_proteome_r2_summary.csv\n")

## ── 10. Source data for figure ───────────────────────────────────────────────

source_data <- data.frame(
  trait = "BMI",
  component = c("Baseline covariates",
                "Olink ∩ SomaScan (triangulable)",
                "Olink-only increment",
                "Residual unexplained"),
  r2 = c(r2_baseline, r2_triangulable, r2_olink_only, r2_residual),
  sample_size_train = length(train_idx),
  sample_size_test = length(test_idx),
  n_proteins_triangulable = length(prot_triangulable),
  n_proteins_olink_only = length(prot_olink_only),
  n_proteins_total = length(keep_prots),
  stringsAsFactors = FALSE
)
write.csv(source_data, file.path(out_dir, "figure2_source_data.csv"), row.names = FALSE)

## ── 11. Figure 2 — Stacked bar chart ────────────────────────────────────────

cat("Generating figure...\n")

plot_df <- data.frame(
  trait = factor("BMI"),
  component = factor(
    c("Baseline covariates",
      "Triangulable proteins\n(Olink ∩ SomaScan)",
      "Olink-only proteins",
      "Residual unexplained"),
    levels = c("Residual unexplained",
               "Olink-only proteins",
               "Triangulable proteins\n(Olink ∩ SomaScan)",
               "Baseline covariates")
  ),
  value = c(r2_baseline, r2_triangulable, r2_olink_only, r2_residual)
)

# Cumulative positions for text labels
plot_df <- plot_df %>%
  arrange(trait, component) %>%
  group_by(trait) %>%
  mutate(
    ymax = cumsum(value),
    ymin = ymax - value,
    ymid = (ymin + ymax) / 2,
    label = ifelse(value >= 0.02, sprintf("%.3f", value), "")
  ) %>%
  ungroup()

# Text colour: white on dark fills, black on light/white
plot_df$text_col <- ifelse(
  plot_df$component %in% c("Baseline covariates",
                           "Triangulable proteins\n(Olink ∩ SomaScan)"),
  "white", "black"
)
# Override: Baseline grey might be light enough for black text
plot_df$text_col[plot_df$component == "Baseline covariates"] <- "black"

fill_colours <- c(
  "Baseline covariates"                           = "#BFBFBF",
  "Triangulable proteins\n(Olink ∩ SomaScan)"     = "#4E79A7",
  "Olink-only proteins"                           = "#F28E2B",
  "Residual unexplained"                          = "#FFFFFF"
)

p <- ggplot(plot_df, aes(x = trait, y = value, fill = component)) +
  geom_col(width = 0.55, colour = "black", linewidth = 0.3) +
  # Hatch pattern for residual: overlay with diagonal lines
  geom_rect(
    data = plot_df %>% filter(component == "Residual unexplained"),
    aes(xmin = as.numeric(trait) - 0.275,
        xmax = as.numeric(trait) + 0.275,
        ymin = ymin, ymax = ymax),
    fill = NA, colour = "#999999", linewidth = 0.3
  ) +
  geom_text(aes(y = ymid, label = label, colour = text_col),
            size = 6 / .pt, family = "sans", show.legend = FALSE) +
  scale_fill_manual(values = fill_colours,
                    breaks = rev(c("Baseline covariates",
                                   "Triangulable proteins\n(Olink ∩ SomaScan)",
                                   "Olink-only proteins",
                                   "Residual unexplained"))) +
  scale_colour_identity() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02)),
                     limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(y = expression(bold("Variance explained (R"^2*")")),
       x = NULL, fill = NULL) +
  theme_classic(base_family = "sans", base_size = 7) +
  theme(
    axis.title.y = element_text(size = 7, face = "bold"),
    axis.text = element_text(size = 7),
    axis.ticks.length = unit(1.5, "pt"),
    axis.ticks = element_line(linewidth = 0.3),
    axis.line = element_line(linewidth = 0.3),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 6),
    legend.key.size = unit(8, "pt"),
    legend.spacing.x = unit(3, "pt"),
    legend.margin = margin(t = 2),
    plot.margin = margin(t = 5, r = 5, b = 2, l = 5)
  ) +
  guides(fill = guide_legend(nrow = 2, reverse = FALSE))

ggsave(file.path(out_dir, "figure2_variance_decomposition.pdf"),
       p, width = 88, height = 100, units = "mm", device = cairo_pdf)
ggsave(file.path(out_dir, "figure2_variance_decomposition.png"),
       p, width = 88, height = 100, units = "mm", dpi = 300)

cat("  Saved figure2_variance_decomposition.pdf/png\n")

## ── 12. Methods paragraph ────────────────────────────────────────────────────

n_total <- length(train_idx) + length(test_idx)
methods_txt <- sprintf(
paste0(
"Variance decomposition analysis. To quantify the proportion of phenotypic variance ",
"in BMI explained by circulating proteins and to assess cross-platform measurability, ",
"we partitioned post-QC Olink proteins into a triangulable set (also measurable on ",
"SomaScan under our reference mapping; n = %d after QC) and Olink-only (n = %d). ",
"SomaScan reference: %s. ",
"We analysed %s UK Biobank participants of European ancestry (random subsample up to ",
"%s individuals for computational feasibility; seed %d). ",
"Per-protein missingness threshold: retain proteins with ≤%.1f%% missing values in this ",
"sample (median imputation for remaining missingness within retained proteins); ",
"inverse-rank-normalisation. Train/test split 80/20 stratified by sex. Elastic net ",
"(alpha = 0.5) with %d-fold cross-validation for lambda; covariates (age, sex, ",
"assessment centre, first 10 PCs) unpenalised. Out-of-sample R-squared on the test set. ",
"Incremental R-squared from triangulable and Olink-only protein blocks via sequential ",
"model comparison."
),
  length(prot_triangulable),
  length(prot_olink_only),
  soma_info$source,
  format(n_total, big.mark = ","),
  format(max_n, big.mark = ","),
  RUN_SEED,
  100 * missingness_threshold_used,
  cv_folds
)

writeLines(methods_txt, file.path(out_dir, "methods_paragraph.txt"))
cat("  Saved methods_paragraph.txt\n")

cat("\n=== DONE ===\n")
cat("Output files in:", out_dir, "\n")
