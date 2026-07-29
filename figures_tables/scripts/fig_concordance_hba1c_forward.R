#!/usr/bin/env Rscript
.libPaths(c("/n/groups/patel/sivateja/R_libs", .libPaths()))
library(dplyr)
library(ggplot2)
library(ggrepel)

cat("=== HbA1c Forward Concordance: Forward Obs + Forward MR vs STEP ===\n")

outdir <- "/n/groups/patel/sivateja/regenie_pipeline/results/twosampleMR/supplemental_tables"

step1 <- read.csv("/n/groups/patel/sivateja/STEP1_merged_results.csv")
step2 <- read.csv("/n/groups/patel/sivateja/STEP2_merged_results.csv")
prot_map <- step1 %>% dplyr::select(Exposure, code) %>% distinct()

bidir <- read.csv(file.path(outdir, "Table_Bidirectional_MR_HbA1c_Full.csv"), check.names = FALSE)
bidir$Fwd_Beta <- as.numeric(bidir[["Fwd: Beta"]])
bidir$Fwd_P <- as.numeric(bidir[["Fwd: P-value"]])
fwd_mr_sig <- bidir %>%
  filter(Fwd_P < 0.05) %>%
  dplyr::select(Protein, Fwd_Beta, Fwd_P) %>%
  rename(code = Protein, Beta_fwdMR = Fwd_Beta)

cat("Forward MR nominally sig:", nrow(fwd_mr_sig), "\n")
cat("Forward MR Bonferroni sig:", sum(bidir$Fwd_P < 0.05/349, na.rm = TRUE), "\n")

step1_fwd <- step1 %>%
  filter(Phenotype == "HbA1c", Subgroup == "non T2D") %>%
  dplyr::select(code, name, estimate, Bonferroni, effect_size, pvalue) %>%
  rename(Beta_fwd_obs = estimate, Bonf_fwd = Bonferroni,
         Estimate_STEP1 = effect_size, Pval_STEP1 = pvalue) %>%
  group_by(code) %>% arrange(Pval_STEP1) %>% slice_head(n = 1) %>% ungroup()

step2_fwd <- step2 %>%
  filter(Phenotype == "HbA1c", Subgroup == "T2D") %>%
  dplyr::select(code, effect_size, pvalue) %>%
  rename(Estimate_STEP2 = effect_size, Pval_STEP2 = pvalue) %>%
  group_by(code) %>% arrange(Pval_STEP2) %>% slice_head(n = 1) %>% ungroup()

plot_df <- step1_fwd %>%
  left_join(step2_fwd, by = "code")

plot_df$ESTIMATE_TRIAL <- plot_df$Estimate_STEP2
plot_df$Pval_TRIAL <- plot_df$Pval_STEP2
trial_label <- "STEP 2"

plot_df <- plot_df %>%
  filter(Bonf_fwd < 0.05, Pval_TRIAL < 0.05) %>%
  filter(!is.na(Beta_fwd_obs), !is.na(ESTIMATE_TRIAL))

plot_df <- plot_df %>%
  left_join(fwd_mr_sig, by = "code")

plot_df$fwd_mr_sig <- !is.na(plot_df$Beta_fwdMR)
plot_df$concordant <- plot_df$fwd_mr_sig &
  (sign(plot_df$Beta_fwd_obs) == sign(plot_df$Beta_fwdMR))

cat("Plot proteins:", nrow(plot_df), "\n")
cat("Forward-MR sig:", sum(plot_df$fwd_mr_sig), "\n")
cat("Concordant (fwd obs + fwd MR same sign):", sum(plot_df$concordant), "\n")

plot_df$highlight <- ifelse(plot_df$concordant, "Concordant Fwd-MR", "Other")

quad_colors <- c("Q3" = "blue", "Q4" = "darkgreen", "Q2" = "red", "Q1" = "darkorange")
plot_df$quadrant <- factor(
  (plot_df$Beta_fwd_obs > 0) + 2 * (plot_df$ESTIMATE_TRIAL > 0),
  levels = 0:3,
  labels = c("Q3", "Q4", "Q2", "Q1"))

label_proteins <- plot_df %>%
  filter(concordant) %>%
  mutate(combined_effect = abs(Beta_fwd_obs) + abs(ESTIMATE_TRIAL) + abs(Beta_fwdMR)) %>%
  group_by(quadrant) %>%
  arrange(desc(combined_effect)) %>%
  slice_head(n = 7) %>%
  ungroup()

if (nrow(label_proteins) < 5) {
  extra <- plot_df %>%
    filter(!concordant) %>%
    mutate(combined_effect = abs(Beta_fwd_obs) + abs(ESTIMATE_TRIAL)) %>%
    group_by(quadrant) %>%
    arrange(desc(combined_effect)) %>%
    slice_head(n = 3) %>%
    ungroup()
  label_proteins <- bind_rows(label_proteins, extra)
}

quadrant_stats <- plot_df %>%
  group_by(quadrant) %>%
  summarise(n_total = n(), n_conc = sum(concordant), .groups = "drop")

x_range <- range(plot_df$Beta_fwd_obs, na.rm = TRUE)
y_range <- range(plot_df$ESTIMATE_TRIAL, na.rm = TRUE)
x_pad <- diff(x_range) * 0.03; y_pad <- diff(y_range) * 0.03

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
qlabels$label <- sprintf("%s\nFwd-MR: %d/%d", qlabels$dir_label, qlabels$n_conc, qlabels$n_total)
qlabels$color <- quad_colors[qlabels$quadrant]

p <- ggplot(plot_df, aes(x = Beta_fwd_obs, y = ESTIMATE_TRIAL, color = quadrant)) +
  geom_point(data = plot_df %>% filter(highlight == "Other"),
             shape = 16, size = 1.5, alpha = 0.3) +
  geom_point(data = plot_df %>% filter(highlight == "Concordant Fwd-MR"),
             shape = 15, size = 4, alpha = 0.9) +
  theme_minimal(base_size = 13) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = quad_colors, guide = "none") +
  geom_text_repel(data = label_proteins, aes(label = code),
                  size = 3.0, fontface = "bold", alpha = 0.85,
                  max.overlaps = 50, show.legend = FALSE) +
  geom_text(data = qlabels, aes(x = x, y = y, label = label),
            hjust = qlabels$hjust, vjust = qlabels$vjust,
            size = 3.8, fontface = "bold", color = qlabels$color, inherit.aes = FALSE) +
  xlab("Forward Observational Effect Size (UKB)") +
  ylab(sprintf("Semaglutide Treatment Effect (%s)", trial_label)) +
  ggtitle("HbA1c Forward Associations vs STEP 2 Results") +
  theme(legend.position = "bottom")

ggsave(file.path(outdir, "Fig_Concordance_HbA1c_Forward.pdf"), p, width = 10, height = 10, dpi = 300)
ggsave(file.path(outdir, "Fig_Concordance_HbA1c_Forward.png"), p, width = 10, height = 10, dpi = 300)
cat("Saved: Fig_Concordance_HbA1c_Forward.pdf/png\n")
