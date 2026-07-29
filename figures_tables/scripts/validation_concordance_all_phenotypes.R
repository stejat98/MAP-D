#!/usr/bin/env Rscript
# =============================================================================
# Concordance Scatter Plots: UKB reverse obs vs STEP treatment effect
# with reverse-MR concordant hits highlighted
# For: BMI, HbA1c, TRIG_HDL_RATIO
# =============================================================================

.libPaths(c("/path/to/project/R_libs", .libPaths()))
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)

cat("=== Concordance Plots for All Phenotypes ===\n\n")

# ---- Shared data ----
step1 <- read.csv("/path/to/project/STEP1_merged_results.csv")
step2 <- read.csv("/path/to/project/STEP2_merged_results.csv")
rev_obs <- read.csv("/path/to/project/UKB/PEWAS_results/Reverse_PEWAS_all_phenotypes_all_strata.csv")
prot_map <- step1 %>% dplyr::select(Exposure, code) %>% distinct()

outdir <- "/path/to/project/regenie_pipeline/results/twosampleMR/supplemental_tables"

# ---- Data layer (shared by single-panel and multi-panel exports) ----
compute_concordance_data <- function(pheno, pheno_label, stratum, mr_file, step_subgroup, trial_mode = "both") {

  if (trial_mode == "step1") trial_label <- "STEP 1"
  else if (trial_mode == "step2") trial_label <- "STEP 2"
  else trial_label <- "STEP 1/2"

  axis_pheno <- switch(pheno, "TRIG_HDL_RATIO" = "TG/HDL-C Ratio", pheno)

  cat(sprintf("\n--- %s (trial: %s) ---\n", pheno_label, trial_label))

  rev_sub <- rev_obs %>%
    filter(Phenotype == pheno, Stratum == stratum) %>%
    left_join(prot_map, by = c("Protein" = "Exposure")) %>%
    filter(!is.na(code))
  cat("  Reverse obs proteins:", nrow(rev_sub), "\n")

  bidir <- read.csv(mr_file, check.names = FALSE)
  rev_beta_col <- grep("Rev.*IVW.*Beta|Rev.*Beta", names(bidir), value = TRUE)[1]
  rev_p_col <- grep("Rev.*IVW.*P|Rev.*P.value", names(bidir), value = TRUE)[1]

  if (is.na(rev_beta_col) || is.na(rev_p_col)) {
    cat("  WARNING: Could not find Rev Beta/P columns. Skipping.\n")
    return(NULL)
  }

  bidir$Rev_Beta <- as.numeric(bidir[[rev_beta_col]])
  bidir$Rev_P <- as.numeric(bidir[[rev_p_col]])
  rev_mr_tested <- bidir %>% filter(!is.na(Rev_P)) %>% pull(Protein) %>% unique()
  rev_mr_sig <- bidir %>% filter(Rev_P < 0.05) %>%
    dplyr::select(Protein, Rev_Beta) %>%
    rename(code = Protein, Beta_revMR = Rev_Beta)

  cat("  Reverse-MR tested:", length(rev_mr_tested), "\n")
  cat("  Nominally sig reverse-MR:", nrow(rev_mr_sig), "\n")

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

  rev_for_merge <- rev_sub %>%
    dplyr::select(code, Beta_obs_rev = estimate, Pval_obs_rev = p.value)

  n_rev <- nrow(rev_for_merge)
  bonf <- 0.05 / n_rev

  plot_df <- step_merged %>%
    inner_join(rev_for_merge, by = "code") %>%
    filter(Pval_obs_rev < bonf, Pval_TRIAL < 0.05) %>%
    left_join(rev_mr_sig, by = "code")

  plot_df <- plot_df %>%
    group_by(code) %>% slice_head(n = 1) %>% ungroup()

  plot_df$rev_mr_tested <- plot_df$code %in% rev_mr_tested
  plot_df$rev_mr_sig <- !is.na(plot_df$Beta_revMR)
  plot_df$concordant <- plot_df$rev_mr_sig & (sign(plot_df$Beta_revMR) == sign(plot_df$Beta_obs_rev))

  cat("  Plot proteins:", nrow(plot_df), "\n")
  cat("  Reverse-MR tested (in plot):", sum(plot_df$rev_mr_tested), "\n")
  cat("  Reverse-MR sig:", sum(plot_df$rev_mr_sig), "\n")
  cat("  Concordant:", sum(plot_df$concordant), "\n")

  if (nrow(plot_df) < 3) {
    cat("  Too few proteins for plot. Skipping.\n")
    return(NULL)
  }

  plot_df$highlight <- ifelse(plot_df$concordant, "Concordant Rev-MR", "Other")

  quad_colors <- c("Q3" = "blue", "Q4" = "darkgreen", "Q2" = "red", "Q1" = "darkorange")
  plot_df$quadrant <- factor(
    (plot_df$Beta_obs_rev > 0) + 2 * (plot_df$ESTIMATE_TRIAL > 0),
    levels = 0:3,
    labels = c("Q3", "Q4", "Q2", "Q1"))

  label_proteins <- plot_df %>%
    filter(concordant) %>%
    mutate(combined_effect = abs(Beta_obs_rev) + abs(ESTIMATE_TRIAL) + abs(Beta_revMR)) %>%
    group_by(quadrant) %>%
    arrange(desc(combined_effect)) %>%
    slice_head(n = 5) %>%
    ungroup()

  quadrant_stats <- plot_df %>%
    group_by(quadrant) %>%
    summarise(n_total = n(), n_tested = sum(rev_mr_tested), n_sig = sum(rev_mr_sig), n_conc = sum(concordant), .groups = "drop")

  cat("\n  === PER-QUADRANT BREAKDOWN ===\n")
  cat(sprintf("  %-12s %8s %10s %8s %10s\n", "Quadrant", "Total", "MR-tested", "MR-sig", "Concordant"))
  for (i in 1:nrow(quadrant_stats)) {
    qs <- quadrant_stats[i,]
    cat(sprintf("  %-12s %8d %10d %8d %10d\n", as.character(qs$quadrant), qs$n_total, qs$n_tested, qs$n_sig, qs$n_conc))
  }
  cat(sprintf("  %-12s %8d %10d %8d %10d\n", "TOTAL",
    sum(quadrant_stats$n_total), sum(quadrant_stats$n_tested),
    sum(quadrant_stats$n_sig), sum(quadrant_stats$n_conc)))
  cat(sprintf("  MR coverage: %d/%d (%.1f%%) of plot proteins were tested in reverse MR\n",
    sum(quadrant_stats$n_tested), sum(quadrant_stats$n_total),
    100 * sum(quadrant_stats$n_tested) / sum(quadrant_stats$n_total)))

  x_range <- range(plot_df$Beta_obs_rev, na.rm = TRUE)
  y_range <- range(plot_df$ESTIMATE_TRIAL, na.rm = TRUE)
  x_pad <- diff(x_range) * 0.02; y_pad <- diff(y_range) * 0.02

  qlabels <- data.frame(
    quadrant = c("Q3", "Q4", "Q2", "Q1"),
    dir_label = c("Persistent change\nafter intervention",
                  "Reversion\nafter intervention",
                  "Reversion\nafter intervention",
                  "Persistent change\nafter intervention"),
    x = c(x_range[1]+x_pad, x_range[2]-x_pad, x_range[1]+x_pad, x_range[2]-x_pad),
    y = c(y_range[1]+y_pad, y_range[1]+y_pad, y_range[2]-y_pad, y_range[2]-y_pad),
    hjust = c(0,1,0,1), vjust = c(0,0,1,1),
    stringsAsFactors = FALSE)
  qlabels <- merge(qlabels, quadrant_stats, by = "quadrant")
  qlabels$label <- sprintf("%s\nRev-MR: %d/%d (n=%d)", qlabels$dir_label, qlabels$n_conc, qlabels$n_tested, qlabels$n_total)
  qlabels$color <- quad_colors[qlabels$quadrant]

  list(
    plot_df = plot_df,
    label_proteins = label_proteins,
    qlabels = qlabels,
    quad_colors = quad_colors,
    axis_pheno = axis_pheno,
    trial_label = trial_label,
    pheno_label = pheno_label
  )
}

# ---- ggplot layer ----
build_concordance_ggplot <- function(dat, panel_mode = c("single", "combined")) {
  panel_mode <- match.arg(panel_mode)
  plot_df <- dat$plot_df
  label_proteins <- dat$label_proteins
  qlabels <- dat$qlabels
  quad_colors <- dat$quad_colors
  axis_pheno <- dat$axis_pheno
  trial_label <- dat$trial_label
  pheno_label <- dat$pheno_label
  intervention_lab <- "Intervention"

  p <- ggplot(plot_df, aes(x = Beta_obs_rev, y = ESTIMATE_TRIAL, color = quadrant)) +
    geom_point(data = plot_df %>% filter(highlight == "Other"),
               shape = 16, size = 1.2, alpha = 0.15) +
    geom_point(data = plot_df %>% filter(highlight == "Concordant Rev-MR"),
               shape = 15, size = 3.5, alpha = 0.9) +
    theme_minimal(base_size = 12) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    scale_color_manual(values = quad_colors, guide = "none") +
    geom_text_repel(data = label_proteins, aes(label = code),
                    size = 2.85, fontface = "bold", alpha = 0.85,
                    max.overlaps = 50, show.legend = FALSE) +
    geom_text(data = qlabels, aes(x = x, y = y, label = label),
              hjust = qlabels$hjust, vjust = qlabels$vjust,
              size = 3.5, fontface = "bold", color = qlabels$color, inherit.aes = FALSE) +
    xlab(bquote("UKB reverse association (" * .(axis_pheno) %->% "protein)")) +
    ylab(bquote("Semaglutide trial effect (" * .(intervention_lab) %->% "protein, " * .(trial_label) * ")")) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "grey85", linewidth = 0.35),
      axis.title = element_text(size = 11, face = "plain"),
      plot.title = element_text(size = 11.5, face = "bold"),
      plot.subtitle = element_text(size = 10, colour = "grey25", margin = margin(b = 6))
    )

  if (panel_mode == "single") {
    p <- p + ggtitle(sprintf(
      "%s: UKB reverse trait-to-protein associations versus %s semaglutide trial effects",
      pheno_label, trial_label))
  } else {
    p <- p + labs(
      title = sprintf("%s (%s)", pheno_label, trial_label),
      subtitle = NULL
    )
  }
  p
}

# ---- Wrapper: compute, plot, save ----
make_concordance_plot <- function(pheno, pheno_label, stratum, mr_file, step_subgroup, trial_mode = "both") {
  dat <- compute_concordance_data(pheno, pheno_label, stratum, mr_file, step_subgroup, trial_mode)
  if (is.null(dat)) return(invisible(NULL))

  p <- build_concordance_ggplot(dat, "single")

  fname <- sprintf("%s/Fig_Concordance_%s.pdf", outdir, pheno)
  ggsave(fname, plot = p, width = 10, height = 10, device = cairo_pdf)
  fname_png <- sprintf("%s/Fig_Concordance_%s.png", outdir, pheno)
  ggsave(fname_png, plot = p, width = 10, height = 10, dpi = 300)
  cat(sprintf("  Saved: %s\n", fname))

  invisible(list(plot = p, data = dat$plot_df, concordance = dat))
}

# ---- Run all three with trial-specific assignments ----
# BMI -> STEP 1 (obesity/non-T2D trial)
bmi_res <- make_concordance_plot(
  "BMI", "Normoglycemic BMI", "non_T2D",
  "/path/to/project/regenie_pipeline/results/twosampleMR/supplemental_tables/Table_Bidirectional_MR_BMI_Full.csv",
  "non T2D", trial_mode = "step1")

# HbA1c -> STEP 2 (T2D trial)
hba1c_res <- make_concordance_plot(
  "HbA1c", "T2D HbA1c", "non_T2D",
  "/path/to/project/regenie_pipeline/results/twosampleMR/supplemental_tables/Table_Bidirectional_MR_HbA1c_Full.csv",
  "non T2D", trial_mode = "step2")

# TG/HDL -> both trials (no single trial specifically targets dyslipidemia)
trig_res <- make_concordance_plot(
  "TRIG_HDL_RATIO", "TG/HDL-C Ratio", "non_T2D",
  "/path/to/project/regenie_pipeline/results/twosampleMR/supplemental_tables/Table_Bidirectional_MR_TRIG_HDL_RATIO_Full.csv",
  "non T2D", trial_mode = "both")

# ---- Combined vertical figure: BMI (A) + HbA1c (B), separate PDF for manuscript ----
if (!is.null(bmi_res) && !is.null(hba1c_res)) {
  p_a <- build_concordance_ggplot(bmi_res$concordance, "combined")
  p_b <- build_concordance_ggplot(hba1c_res$concordance, "combined")

  fig_combined <- p_a / p_b +
    plot_layout(ncol = 1, heights = c(1, 1)) +
    plot_annotation(
      tag_levels = "A",
      theme = theme(
        plot.tag = element_text(size = 13, face = "bold"),
        plot.margin = margin(10, 10, 10, 10)
      )
    )

  combined_pdf <- file.path(outdir, "Fig_Concordance_BMI_HbA1c_UKB_panels_AB.pdf")
  combined_png <- file.path(outdir, "Fig_Concordance_BMI_HbA1c_UKB_panels_AB.png")

  ggsave(combined_pdf, plot = fig_combined, width = 7.167, height = 14.25, device = cairo_pdf)
  ggsave(combined_png, plot = fig_combined, width = 7.167, height = 14.25, dpi = 600)

  cat("\n=== Combined panel figure (A: BMI, B: HbA1c) ===\n")
  cat("  Saved:", combined_pdf, "\n")
  cat("  Saved:", combined_png, "\n\n")

  caption_text <- paste(
    "Concordance of UK Biobank reverse trait-to-protein association estimates with",
    "semaglutide trial effects on the same proteins. Each point is a protein passing",
    "reverse-association significance in UKB (Bonferroni) and nominal trial significance;",
    "axis ranges are panel-specific. Squares highlight proteins with nominally significant",
    "reverse Mendelian randomization (MR) concordant in sign with the observational",
    "reverse association. (a) Normoglycemic BMI (non-T2D stratum) versus STEP 1.",
    "(b) HbA1c (non-T2D stratum) versus STEP 2. Dashed lines, null; quadrant shading,",
    "direction of UKB versus trial effects."
  )
  cat("--- Suggested figure caption (paste into manuscript) ---\n")
  cat(caption_text, "\n")
  cat("--- end caption ---\n")
}

cat("\n=== All concordance plots done ===\n")
