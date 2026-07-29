#!/usr/bin/env Rscript
.libPaths(c("/n/groups/patel/sivateja/R_libs", .libPaths()))
library(dplyr)
library(ggplot2)
library(ggrepel)

cat("=== Forward Obs vs STEP Concordance Figures (Updated Labels) ===\n")

outdir <- "/n/groups/patel/sivateja/regenie_pipeline/results/twosampleMR/supplemental_tables"

step1 <- read.csv("/n/groups/patel/sivateja/STEP1_merged_results.csv")
step2 <- read.csv("/n/groups/patel/sivateja/STEP2_merged_results.csv")

step1val_ukb <- step1 %>%
  dplyr::select(Phenotype, Subgroup, EntrezGeneID, name, estimate, Bonferroni, Exposure, code) %>%
  setNames(c("Phenotype","Subgroup","EntrezGeneID","Protein_UKB","Estimate_UKB","Significance_UKB","Exposure_UKB","code_UKB"))
step2val_ukb <- step2 %>%
  dplyr::select(Phenotype, Subgroup, EntrezGeneID, name, estimate, Bonferroni, Exposure, code) %>%
  setNames(c("Phenotype","Subgroup","EntrezGeneID","Protein_UKB","Estimate_UKB","Significance_UKB","Exposure_UKB","code_UKB"))
ukb <- bind_rows(step1val_ukb, step2val_ukb) %>% distinct()

step1_eff <- step1 %>%
  filter(!is.na(EntrezGeneID)) %>%
  dplyr::select(ANALYTEID, Phenotype, Subgroup, EntrezGeneID, TargetFullName, effect_size, pvalue) %>%
  group_by(ANALYTEID, Phenotype, Subgroup, EntrezGeneID, TargetFullName) %>%
  summarise(effect_size = mean(effect_size), pvalue = mean(pvalue), .groups = "drop") %>%
  setNames(c("ANALYTEID","Phenotype","Subgroup","EntrezGeneID","Protein_STEP1","Estimate_STEP1","Significance_STEP1"))

step2_eff <- step2 %>%
  filter(!is.na(EntrezGeneID)) %>%
  dplyr::select(ANALYTEID, Phenotype, Subgroup, EntrezGeneID, TargetFullName, effect_size, pvalue) %>%
  group_by(ANALYTEID, Phenotype, Subgroup, EntrezGeneID, TargetFullName) %>%
  summarise(effect_size = mean(effect_size), pvalue = mean(pvalue), .groups = "drop") %>%
  setNames(c("ANALYTEID","Phenotype","Subgroup","EntrezGeneID","Protein_STEP2","Estimate_STEP2","Significance_STEP2"))

valdata <- full_join(step1_eff, step2_eff) %>% distinct()
merged_val <- full_join(ukb, valdata)
merged_val$ESTIMATE_STEP_MEAN <- rowMeans(merged_val[, c("Estimate_STEP1","Estimate_STEP2")], na.rm = TRUE)

merged_val_sig1 <- merged_val %>% filter(Significance_UKB < 0.05, Significance_STEP1 < 0.05)
merged_val_sig2 <- merged_val %>% filter(Significance_UKB < 0.05, Significance_STEP2 < 0.05)

quad_colors <- c("Q3" = "blue", "Q4" = "darkgreen", "Q2" = "red", "Q1" = "darkorange")

make_concordance_plot <- function(df, pheno, subgroup, trial_label, pheno_label) {
  sub <- df %>% filter(Phenotype == pheno, Subgroup == subgroup) %>%
    filter(!is.na(Estimate_UKB), !is.na(ESTIMATE_STEP_MEAN))

  if (nrow(sub) == 0) {
    cat(sprintf("  SKIP: %s / %s — no data\n", pheno_label, subgroup))
    return(NULL)
  }

  sub$quadrant <- factor(
    (sub$Estimate_UKB > 0) + 2 * (sub$ESTIMATE_STEP_MEAN > 0),
    levels = 0:3, labels = c("Q3","Q4","Q2","Q1"))

  top_proteins <- sub %>%
    group_by(quadrant) %>%
    arrange(desc(abs(Estimate_UKB) + abs(ESTIMATE_STEP_MEAN))) %>%
    distinct(Protein_UKB, .keep_all = TRUE) %>%
    slice_head(n = 5) %>%
    ungroup()

  quadrant_stats <- sub %>%
    group_by(quadrant) %>%
    summarise(n_total = n(), .groups = "drop")

  x_range <- range(sub$Estimate_UKB, na.rm = TRUE)
  y_range <- range(sub$ESTIMATE_STEP_MEAN, na.rm = TRUE)
  x_pad <- diff(x_range) * 0.03
  y_pad <- diff(y_range) * 0.03

  qlabels <- data.frame(
    quadrant = c("Q3","Q4","Q2","Q1"),
    dir_label = c("Persistent change\nafter intervention",
                  "Reversion\nafter intervention",
                  "Reversion\nafter intervention",
                  "Persistent change\nafter intervention"),
    x = c(x_range[1]+x_pad, x_range[2]-x_pad, x_range[1]+x_pad, x_range[2]-x_pad),
    y = c(y_range[1]+y_pad, y_range[1]+y_pad, y_range[2]-y_pad, y_range[2]-y_pad),
    hjust = c(0,1,0,1), vjust = c(0,0,1,1),
    stringsAsFactors = FALSE)
  qlabels <- merge(qlabels, quadrant_stats, by = "quadrant")
  qlabels$label <- sprintf("%s\nn = %d", qlabels$dir_label, qlabels$n_total)
  qlabels$color <- quad_colors[qlabels$quadrant]

  p <- ggplot(sub, aes(x = Estimate_UKB, y = ESTIMATE_STEP_MEAN, color = quadrant)) +
    geom_point(size = 1.2, alpha = 0.5) +
    theme_minimal(base_size = 13) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    scale_color_manual(values = quad_colors, guide = "none") +
    geom_text_repel(data = top_proteins, aes(label = Protein_UKB),
                    size = 2.8, fontface = "bold", alpha = 0.85,
                    max.overlaps = 50, show.legend = FALSE) +
    geom_text(data = qlabels, aes(x = x, y = y, label = label),
              hjust = qlabels$hjust, vjust = qlabels$vjust,
              size = 2.8, fontface = "bold", color = qlabels$color, inherit.aes = FALSE) +
    xlab("Forward Observational Effect Size (UKB)") +
    ylab(sprintf("Semaglutide Treatment Effect (%s)", trial_label)) +
    ggtitle(sprintf("%s Associations vs %s Results", pheno_label, trial_label)) +
    theme(legend.position = "none")

  cat(sprintf("  %s: %d proteins plotted\n", pheno_label, nrow(sub)))
  return(p)
}

configs <- list(
  list(df = "merged_val_sig1", pheno = "BMI", subgroup = "non T2D",
       trial = "STEP 1", label = "Normoglycemic BMI",
       fname = "Fig_Forward_Obs_vs_STEP1_BMI"),
  list(df = "merged_val_sig1", pheno = "TRIG_HDL_RATIO", subgroup = "non T2D",
       trial = "STEP 1", label = "Normoglycemic TG/HDL-C Ratio",
       fname = "Fig_Forward_Obs_vs_STEP1_TRIG_HDL"),
  list(df = "merged_val_sig2", pheno = "HbA1c", subgroup = "T2D",
       trial = "STEP 2", label = "T2D HbA1c",
       fname = "Fig_Forward_Obs_vs_STEP2_HbA1c"),
  list(df = "merged_val_sig2", pheno = "TRIG_HDL_RATIO", subgroup = "T2D",
       trial = "STEP 2", label = "T2D TG/HDL-C Ratio",
       fname = "Fig_Forward_Obs_vs_STEP2_TRIG_HDL"),
  list(df = "merged_val_sig1", pheno = "HbA1c", subgroup = "non T2D",
       trial = "STEP 1", label = "Normoglycemic HbA1c",
       fname = "Fig_Forward_Obs_vs_STEP1_HbA1c")
)

for (cfg in configs) {
  dat <- if (cfg$df == "merged_val_sig1") merged_val_sig1 else merged_val_sig2
  p <- make_concordance_plot(dat, cfg$pheno, cfg$subgroup, cfg$trial, cfg$label)
  if (!is.null(p)) {
    ggsave(file.path(outdir, paste0(cfg$fname, ".pdf")), p, width = 8, height = 8, dpi = 300)
    ggsave(file.path(outdir, paste0(cfg$fname, ".png")), p, width = 8, height = 8, dpi = 300)
    cat(sprintf("  Saved: %s.pdf/png\n", cfg$fname))
  }
}

cat("\n=== All forward obs vs STEP figures complete ===\n")
