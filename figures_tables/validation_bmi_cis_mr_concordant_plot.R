.libPaths(c("/path/to/project/R_libs", .libPaths()))
library(dplyr)
library(ggplot2)
library(ggrepel)

# ---- Load STEP 1/2 trial data (Y-axis: treatment effect on protein) ----
step1 = read.delim('/path/to/project/STEP1_merged_results.csv', sep=',')
step2 = read.delim('/path/to/project/STEP2_merged_results.csv', sep=',')

# ---- Load reverse observational PEWAS: BMI -> protein (normoglycemic) ----
rev_obs <- read.csv('/path/to/project/UKB/PEWAS_results/Reverse_PEWAS_all_phenotypes_all_strata.csv')
rev_obs_bmi <- rev_obs %>% filter(Phenotype == "BMI", Stratum == "non_T2D")

# Map protein_x.N to gene symbols using STEP1 mapping
prot_map <- step1 %>% dplyr::select(Exposure, code) %>% distinct()
rev_obs_bmi <- rev_obs_bmi %>%
  left_join(prot_map, by = c("Protein" = "Exposure")) %>%
  filter(!is.na(code))

cat("Reverse observational BMI->protein (normoglycemic):", nrow(rev_obs_bmi), "proteins\n")

# ---- Load bidirectional MR: reverse direction (BMI -> protein) ----
bidir_mr = read.csv('/path/to/project/regenie_pipeline/results/twosampleMR/supplemental_tables/Table_Bidirectional_MR_BMI_Full.csv')
bidir_mr$Rev_Beta <- as.numeric(bidir_mr$Rev..IVW.Beta)
bidir_mr$Rev_P    <- as.numeric(bidir_mr$Rev..IVW.P.value)

rev_mr_sig <- bidir_mr %>% filter(Rev_P < 0.05)
cat("Nominally significant reverse-MR proteins (BMI->protein):", nrow(rev_mr_sig), "\n")

# ---- Build STEP 1/2 treatment effect estimates ----
# Deduplicate to one row per gene (code) BEFORE joining STEP1 and STEP2.
# For genes with multiple aptamers, keep the one with the best p-value.
step1val <- step1 %>%
  filter(!is.na(EntrezGeneID), Phenotype == "BMI", Subgroup == "non T2D") %>%
  dplyr::select(code, effect_size, pvalue) %>%
  group_by(code) %>%
  arrange(pvalue) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  rename(Estimate_STEP1 = effect_size, Pval_STEP1 = pvalue)

step2val <- step2 %>%
  filter(!is.na(EntrezGeneID), Phenotype == "BMI", Subgroup == "non T2D") %>%
  dplyr::select(code, effect_size, pvalue) %>%
  group_by(code) %>%
  arrange(pvalue) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  rename(Estimate_STEP2 = effect_size, Pval_STEP2 = pvalue)

cat("STEP1 unique genes:", nrow(step1val), "\n")
cat("STEP2 unique genes:", nrow(step2val), "\n")

step_merged <- full_join(step1val, step2val, by = "code")
step_merged$ESTIMATE_STEP_MEAN <- rowMeans(step_merged[, c("Estimate_STEP1", "Estimate_STEP2")], na.rm = TRUE)

cat("STEP treatment effects (unique genes):", nrow(step_merged), "\n")

# ---- Merge: reverse obs UKB + STEP treatment effect + reverse MR ----
# Use actual reverse PEWAS p-values (not forward model — they differ in practice)
rev_obs_for_merge <- rev_obs_bmi %>%
  dplyr::select(code, estimate, p.value) %>%
  rename(Beta_obs_rev = estimate, Pval_obs_rev = p.value)

n_rev_proteins <- nrow(rev_obs_for_merge)
bonf_thresh <- 0.05 / n_rev_proteins
cat("Reverse obs Bonferroni threshold (0.05 /", n_rev_proteins, "):", bonf_thresh, "\n")

rev_mr_for_merge <- rev_mr_sig %>%
  dplyr::select(Protein, Rev_Beta) %>%
  rename(code = Protein, Beta_revMR = Rev_Beta)

# Verify no duplicates before joining
stopifnot("Duplicate codes in step_merged" = !any(duplicated(step_merged$code)))
stopifnot("Duplicate codes in rev_obs" = !any(duplicated(rev_obs_for_merge$code)))
stopifnot("Duplicate codes in rev_mr" = !any(duplicated(rev_mr_for_merge$code)))

plot_df <- step_merged %>%
  inner_join(rev_obs_for_merge, by = "code")

cat("After merging STEP + reverse obs:", nrow(plot_df), "proteins\n")

# Filter: reverse obs UKB Bonferroni sig AND STEP1 p < 0.05
plot_df <- plot_df %>%
  filter(Pval_obs_rev < bonf_thresh, Pval_STEP1 < 0.05)

cat("After significance filter (rev obs Bonf < 0.05 & STEP1 p < 0.05):", nrow(plot_df), "proteins\n")

# Add reverse MR data (left join — not all will have MR)
plot_df <- plot_df %>% left_join(rev_mr_for_merge, by = "code")

stopifnot("Duplicates in final plot_df" = !any(duplicated(plot_df$code)))

# ---- Concordance: reverse MR significant AND same direction as reverse obs UKB ----
plot_df$rev_mr_sig <- !is.na(plot_df$Beta_revMR)
plot_df$concordant <- plot_df$rev_mr_sig & (sign(plot_df$Beta_revMR) == sign(plot_df$Beta_obs_rev))

cat("\nReverse-MR significant:", sum(plot_df$rev_mr_sig), "\n")
cat("Concordant (sig rev-MR + same sign as rev obs UKB):", sum(plot_df$concordant), "\n")
cat("Discordant (sig rev-MR + opposite sign):", sum(plot_df$rev_mr_sig & !plot_df$concordant), "\n")

# ---- Point type for plot ----
plot_df$highlight <- ifelse(plot_df$concordant, "Concordant Rev-MR Hit", "Other")

# ---- Quadrant labels ----
plot_df$quadrant <- factor(
  (plot_df$Beta_obs_rev > 0) + 2 * (plot_df$ESTIMATE_STEP_MEAN > 0),
  levels = c(0, 1, 2, 3),
  labels = c("Decreased UKB / Decreased Step 1", "Increased UKB / Decreased Step 1",
             "Decreased UKB / Increased Step 1", "Increased UKB / Increased Step 1")
)

# ---- Label top 5 concordant reverse-MR hits per quadrant ----
label_proteins <- plot_df %>%
  filter(concordant) %>%
  mutate(combined_effect = abs(Beta_obs_rev) + abs(ESTIMATE_STEP_MEAN) + abs(Beta_revMR)) %>%
  group_by(quadrant) %>%
  arrange(desc(combined_effect)) %>%
  slice_head(n = 5) %>%
  ungroup()

cat("\n=== TOP 5 CONCORDANT REVERSE-MR PROTEINS PER QUADRANT ===\n")
for (q in levels(plot_df$quadrant)) {
  qdat <- label_proteins %>% filter(quadrant == q)
  cat(sprintf("\n  %s (%d):\n", q, nrow(qdat)))
  for (i in 1:nrow(qdat)) {
    r <- qdat[i,]
    cat(sprintf("    %d. %s: RevObs=%.4f, STEP_mean=%.3f, RevMR=%.4f\n",
                i, r$code, r$Beta_obs_rev, r$ESTIMATE_STEP_MEAN, r$Beta_revMR))
  }
}

# ---- Quadrant annotations ----
quadrant_stats <- plot_df %>%
  group_by(quadrant) %>%
  summarise(
    n_total = n(),
    n_conc = sum(concordant),
    .groups = "drop"
  )

quadrant_stats$label <- sprintf(
  "MR: %d/%d (%.1f%%)",
  quadrant_stats$n_conc, quadrant_stats$n_total,
  100 * quadrant_stats$n_conc / quadrant_stats$n_total
)

x_range <- range(plot_df$Beta_obs_rev, na.rm = TRUE)
y_range <- range(plot_df$ESTIMATE_STEP_MEAN, na.rm = TRUE)
x_pad <- diff(x_range) * 0.02
y_pad <- diff(y_range) * 0.02

quadrant_positions <- data.frame(
  quadrant = c("Decreased UKB / Decreased Step 1",
               "Increased UKB / Decreased Step 1",
               "Decreased UKB / Increased Step 1",
               "Increased UKB / Increased Step 1"),
  x = c(x_range[1] + x_pad, x_range[2] - x_pad,
        x_range[1] + x_pad, x_range[2] - x_pad),
  y = c(y_range[1] + y_pad, y_range[1] + y_pad,
        y_range[2] - y_pad, y_range[2] - y_pad),
  hjust = c(0, 1, 0, 1),
  vjust = c(0, 0, 1, 1)
)

quadrant_labels <- merge(quadrant_stats, quadrant_positions, by = "quadrant")

cat("\n=== QUADRANT SUMMARY ===\n")
for (i in 1:nrow(quadrant_labels)) {
  cat(sprintf("  %s: %s\n", quadrant_labels$quadrant[i], quadrant_labels$label[i]))
}

# ---- Plot ----
p <- ggplot(plot_df, aes(x = Beta_obs_rev, y = ESTIMATE_STEP_MEAN, color = quadrant)) +
  geom_point(data = plot_df %>% filter(highlight == "Other"),
             shape = 16, size = 1.2, alpha = 0.15) +
  geom_point(data = plot_df %>% filter(highlight == "Concordant Rev-MR Hit"),
             aes(shape = "BMI-associated & reverted post-intervention"),
             size = 3.5, alpha = 0.9) +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'red') +
  scale_color_manual(values = c("blue", "darkgreen", "red", "darkorange"), guide = "none") +
  scale_shape_manual(name = NULL, values = c("BMI-associated & reverted post-intervention" = 15)) +
  geom_text_repel(data = label_proteins, aes(label = code),
                  size = 3.0, fontface = "bold", alpha = 0.85,
                  max.overlaps = 50, show.legend = FALSE) +
  geom_text(data = quadrant_labels,
            aes(x = x, y = y, label = label),
            hjust = quadrant_labels$hjust,
            vjust = quadrant_labels$vjust,
            size = 3.0, fontface = "bold", color = "grey30",
            inherit.aes = FALSE) +
  xlab('Estimate UKB') +
  ylab('Estimate Step 1') +
  ggtitle('Normoglycemic BMI Associations vs. Step 1 Results') +
  theme(legend.position = 'none')

ggsave('/path/to/project/validation_step1_points_BMI_concordant.pdf', plot = p, width = 10, height = 10)
cat("\nPlot saved to /path/to/project/validation_step1_points_BMI_concordant.pdf\n")
