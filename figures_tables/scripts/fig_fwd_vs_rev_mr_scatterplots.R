#!/usr/bin/env Rscript
.libPaths(c("/n/groups/patel/sivateja/R_libs", .libPaths()))
library(ggplot2)
library(dplyr)
library(ggrepel)

setwd("/n/groups/patel/sivateja/regenie_pipeline")
outdir <- "results/twosampleMR/supplemental_tables"

bmi <- read.csv("results/twosampleMR/bidirectional_bmi_proteins/bidirectional_MR_BMI_all_proteins.csv")
hba1c <- read.csv("results/twosampleMR/bidirectional_hba1c_proteins/bidirectional_MR_HbA1c_all_proteins.csv")

hba1c <- hba1c %>% rename(Rev_Beta = IVW_Beta, Rev_SE = IVW_SE, Rev_P = IVW_P)

n_prot <- nrow(bmi)
bonf <- 0.05 / n_prot

categorize <- function(df) {
  df %>% mutate(
    fwd_sig = !is.na(Fwd_P) & Fwd_P < bonf,
    rev_sig = !is.na(Rev_P) & Rev_P < bonf,
    cat = case_when(
      fwd_sig & rev_sig & Concordant == TRUE  ~ "Concordant",
      fwd_sig & rev_sig & Concordant == FALSE ~ "Discordant",
      fwd_sig & !rev_sig ~ "Forward only",
      !fwd_sig & rev_sig ~ "Reverse only",
      TRUE ~ "NS"
    ),
    cat = factor(cat, levels = c("Concordant", "Discordant", "Forward only", "Reverse only", "NS"))
  )
}

bmi <- categorize(bmi)
hba1c <- categorize(hba1c)

cat("BMI:\n"); print(table(bmi$cat))
cat("\nHbA1c:\n"); print(table(hba1c$cat))

cols <- c("Concordant" = "#C0392B", "Discordant" = "#8E44AD",
          "Forward only" = "#2980B9", "Reverse only" = "#27AE60", "NS" = "grey75")
shps <- c("Concordant" = 18, "Discordant" = 4,
          "Forward only" = 16, "Reverse only" = 16, "NS" = 16)

pick_labels <- function(df, n = 3) {
  df %>%
    filter(cat != "NS") %>%
    group_by(cat) %>%
    arrange(pmin(Fwd_P, Rev_P, na.rm = TRUE)) %>%
    slice_head(n = n) %>%
    ungroup()
}

# ---- Beta scatterplot ----
plot_beta <- function(df, trait) {
  df_complete <- df %>% filter(!is.na(Fwd_Beta), !is.na(Rev_Beta))
  labs_df <- pick_labels(df_complete)

  ggplot(df_complete, aes(Fwd_Beta, Rev_Beta)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_point(data = df_complete %>% filter(cat == "NS"),
               colour = "grey75", size = 0.6, alpha = 0.35) +
    geom_point(data = df_complete %>% filter(cat != "NS"),
               aes(colour = cat, shape = cat), size = 1.8, alpha = 0.85) +
    geom_text_repel(data = labs_df, aes(label = Protein),
                    size = 2.0, segment.size = 0.2, segment.colour = "grey50",
                    max.overlaps = 20, box.padding = 0.35, seed = 1,
                    min.segment.length = 0) +
    scale_colour_manual(values = cols, drop = FALSE) +
    scale_shape_manual(values = shps, drop = FALSE) +
    coord_cartesian(clip = "off") +
    labs(x = bquote("Forward MR" ~ hat(beta) ~ "(Protein" %->% .(trait) * ")"),
         y = bquote("Reverse MR" ~ hat(beta) ~ "(" * .(trait) %->% "Protein)")) +
    theme_classic(base_size = 7) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.3, colour = "black"),
      axis.text = element_text(colour = "black"),
      legend.position = "none",
      plot.margin = margin(8, 10, 4, 6)
    )
}

# ---- P-value scatterplot ----
plot_pval <- function(df, trait) {
  df_complete <- df %>%
    filter(!is.na(Fwd_P), !is.na(Rev_P)) %>%
    mutate(nlp_fwd = -log10(pmax(Fwd_P, 1e-300)),
           nlp_rev = -log10(pmax(Rev_P, 1e-300)))
  labs_df <- pick_labels(df_complete)
  bonf_line <- -log10(bonf)

  ggplot(df_complete, aes(nlp_fwd, nlp_rev)) +
    geom_hline(yintercept = bonf_line, linewidth = 0.3, linetype = "dashed", colour = "firebrick3") +
    geom_vline(xintercept = bonf_line, linewidth = 0.3, linetype = "dashed", colour = "firebrick3") +
    geom_point(data = df_complete %>% filter(cat == "NS"),
               colour = "grey75", size = 0.6, alpha = 0.35) +
    geom_point(data = df_complete %>% filter(cat != "NS"),
               aes(colour = cat, shape = cat), size = 1.8, alpha = 0.85) +
    geom_text_repel(data = labs_df, aes(label = Protein),
                    size = 2.0, segment.size = 0.2, segment.colour = "grey50",
                    max.overlaps = 20, box.padding = 0.35, seed = 1,
                    min.segment.length = 0) +
    scale_colour_manual(values = cols, drop = FALSE) +
    scale_shape_manual(values = shps, drop = FALSE) +
    coord_cartesian(clip = "off") +
    labs(x = bquote(-log[10] * "(P), Forward MR (Protein" %->% .(trait) * ")"),
         y = bquote(-log[10] * "(P), Reverse MR (" * .(trait) %->% "Protein)")) +
    theme_classic(base_size = 7) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.3, colour = "black"),
      axis.text = element_text(colour = "black"),
      legend.position = "none",
      plot.margin = margin(8, 10, 4, 6)
    )
}

# ---- Generate and save ----
p1 <- plot_beta(bmi, "BMI")
p2 <- plot_beta(hba1c, "HbA1c")
p3 <- plot_pval(bmi, "BMI")
p4 <- plot_pval(hba1c, "HbA1c")

ggsave(file.path(outdir, "Fig_FwdVsRev_Beta_BMI.pdf"),   p1, width = 4.5, height = 4.2, device = cairo_pdf)
ggsave(file.path(outdir, "Fig_FwdVsRev_Beta_HbA1c.pdf"), p2, width = 4.5, height = 4.2, device = cairo_pdf)
ggsave(file.path(outdir, "Fig_FwdVsRev_Pval_BMI.pdf"),   p3, width = 4.5, height = 4.2, device = cairo_pdf)
ggsave(file.path(outdir, "Fig_FwdVsRev_Pval_HbA1c.pdf"), p4, width = 4.5, height = 4.2, device = cairo_pdf)

cat("\nDone.\n")
