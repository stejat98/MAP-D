#!/usr/bin/env Rscript
# =============================================================================
# Figure: Bidirectional MR Summary across BMI, HbA1c, TG/HDL-C Ratio
# Shows the asymmetry in causal direction across phenotypes
# =============================================================================

.libPaths(c("/n/groups/patel/sivateja/R_libs", .libPaths()))
library(ggplot2)
library(dplyr)
library(tidyr)

setwd("/n/groups/patel/sivateja/regenie_pipeline")

cat("=== Generating Bidirectional MR Summary Figure ===\n\n")

# ---- Load bidirectional MR results ----
bmi <- read.csv("results/twosampleMR/bidirectional_bmi_proteins/bidirectional_MR_BMI_all_proteins.csv")
hba1c <- read.csv("results/twosampleMR/bidirectional_hba1c_proteins/bidirectional_MR_HbA1c_all_proteins.csv")
trig <- read.csv("results/twosampleMR/bidirectional_trig_hdl_ratio_proteins/bidirectional_MR_TRIG_HDL_RATIO_all_proteins.csv")

classify <- function(df, pheno_label) {
  n <- nrow(df)
  bonf <- 0.05 / n

  fwd_sig <- !is.na(df$Fwd_P) & as.numeric(df$Fwd_P) < bonf
  
  rev_p_col <- if ("Rev_P" %in% names(df)) "Rev_P" else if ("IVW_P" %in% names(df)) "IVW_P" else NA
  if (is.na(rev_p_col)) {
    rev_sig <- rep(FALSE, n)
  } else {
    rev_sig <- !is.na(df[[rev_p_col]]) & as.numeric(df[[rev_p_col]]) < bonf
  }

  conc_col <- if ("Concordant" %in% names(df)) "Concordant" else NULL
  concordant <- if (!is.null(conc_col)) !is.na(df[[conc_col]]) & df[[conc_col]] == TRUE else rep(FALSE, n)

  data.frame(
    Phenotype = pheno_label,
    Category = c("Forward-Only\n(Protein → Trait)",
                 "Reverse-Only\n(Trait → Protein)",
                 "Bidirectional\nConcordant",
                 "Bidirectional\nDiscordant"),
    Count = c(
      sum(fwd_sig & !rev_sig),
      sum(!fwd_sig & rev_sig),
      sum(fwd_sig & rev_sig & concordant),
      sum(fwd_sig & rev_sig & !concordant)
    ),
    stringsAsFactors = FALSE
  )
}

summary_df <- bind_rows(
  classify(bmi, "BMI"),
  classify(hba1c, "HbA1c"),
  classify(trig, "TG/HDL-C\nRatio")
)

summary_df$Phenotype <- factor(summary_df$Phenotype, levels = c("BMI", "HbA1c", "TG/HDL-C\nRatio"))
summary_df$Category <- factor(summary_df$Category,
  levels = c("Forward-Only\n(Protein → Trait)",
             "Reverse-Only\n(Trait → Protein)",
             "Bidirectional\nConcordant",
             "Bidirectional\nDiscordant"))

cat("Summary:\n")
print(summary_df)

# ---- Plot ----
p <- ggplot(summary_df, aes(x = Category, y = Count, fill = Phenotype)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = Count), position = position_dodge(width = 0.8),
            vjust = -0.5, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c("BMI" = "#E64B35", "HbA1c" = "#4DBBD5", "TG/HDL-C\nRatio" = "#00A087")) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    axis.text.x = element_text(size = 11),
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40")
  ) +
  labs(
    title = "Bidirectional Mendelian Randomization: Causal Direction Asymmetry",
    subtitle = "Bonferroni-significant proteins across 349 tested (split-sample UKB cis-pQTLs)",
    x = NULL,
    y = "Number of Proteins"
  ) +
  ylim(0, max(summary_df$Count) * 1.15)

ggsave("results/twosampleMR/supplemental_tables/Fig_Bidirectional_MR_Summary.pdf",
       plot = p, width = 10, height = 6)
ggsave("results/twosampleMR/supplemental_tables/Fig_Bidirectional_MR_Summary.png",
       plot = p, width = 10, height = 6, dpi = 300)

cat("\nFigure saved.\n")

# ---- Also make a stacked proportional version ----
prop_df <- summary_df %>%
  group_by(Phenotype) %>%
  mutate(Total = sum(Count), Pct = ifelse(Total > 0, 100 * Count / Total, 0)) %>%
  ungroup()

p2 <- ggplot(prop_df %>% filter(Count > 0),
             aes(x = Phenotype, y = Pct, fill = Category)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%d", Count)),
            position = position_stack(vjust = 0.5), size = 4, fontface = "bold", color = "white") +
  scale_fill_manual(values = c(
    "Forward-Only\n(Protein → Trait)" = "#3C5488",
    "Reverse-Only\n(Trait → Protein)" = "#E64B35",
    "Bidirectional\nConcordant" = "#F39B7F",
    "Bidirectional\nDiscordant" = "#8491B4"
  )) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5)
  ) +
  labs(
    title = "Causal Architecture Differs Across Metabolic Hallmarks",
    x = NULL,
    y = "Proportion of Significant Proteins (%)"
  )

ggsave("results/twosampleMR/supplemental_tables/Fig_Bidirectional_MR_Proportions.pdf",
       plot = p2, width = 9, height = 6)
ggsave("results/twosampleMR/supplemental_tables/Fig_Bidirectional_MR_Proportions.png",
       plot = p2, width = 9, height = 6, dpi = 300)

cat("Proportional figure saved.\nDone.\n")
