#!/usr/bin/env Rscript
.libPaths(c("/n/groups/patel/sivateja/R_libs", .libPaths()))
library(dplyr)
library(ggplot2)
library(ggrepel)

cat("=== Reverse PEWAS Summary Figure (Fig 1b) ===\n")

outdir <- "/n/groups/patel/sivateja/regenie_pipeline/results/twosampleMR/supplemental_tables"

rev_obs <- read.csv("/n/groups/patel/sivateja/UKB/PEWAS_results/Reverse_PEWAS_all_phenotypes_all_strata.csv")

step1_map <- read.csv("/n/groups/patel/sivateja/STEP1_merged_results.csv") %>%
  dplyr::select(Exposure, code) %>% distinct()
rev_obs <- rev_obs %>% left_join(step1_map, by = c("Protein" = "Exposure"))

pheno_map <- c("BMI" = "BMI", "HbA1c" = "HbA1c", "TRIG_HDL_RATIO" = "TG/HDL-C Ratio")
df <- rev_obs %>%
  filter(Stratum == "non_T2D", Phenotype %in% names(pheno_map)) %>%
  mutate(Phenotype_label = pheno_map[Phenotype])

df <- df %>%
  group_by(Phenotype) %>%
  mutate(Bonferroni = p.adjust(p.value, method = "bonferroni")) %>%
  ungroup()

df$sig <- df$Bonferroni < 0.05
df$neglog10p <- -log10(df$p.value)
df$neglog10p[df$neglog10p > 300] <- 300

cat("Summary per phenotype:\n")
for (ph in names(pheno_map)) {
  sub <- df %>% filter(Phenotype == ph)
  cat(sprintf("  %s: %d total, %d Bonferroni sig\n", ph, nrow(sub), sum(sub$sig)))
}

top_labels <- df %>%
  filter(sig) %>%
  group_by(Phenotype_label) %>%
  arrange(desc(abs(estimate))) %>%
  slice_head(n = 8) %>%
  ungroup()

top_labels$gene <- ifelse(!is.na(top_labels$code), top_labels$code, sub(";.*", "", top_labels$Protein))

df$Phenotype_label <- factor(df$Phenotype_label, levels = c("BMI", "HbA1c", "TG/HDL-C Ratio"))
top_labels$Phenotype_label <- factor(top_labels$Phenotype_label, levels = c("BMI", "HbA1c", "TG/HDL-C Ratio"))

p <- ggplot(df, aes(x = estimate, y = neglog10p)) +
  geom_point(data = df %>% filter(!sig), color = "grey70", size = 0.4, alpha = 0.4) +
  geom_point(data = df %>% filter(sig & estimate > 0), color = "#E64B35", size = 0.8, alpha = 0.6) +
  geom_point(data = df %>% filter(sig & estimate < 0), color = "#4DBBD5", size = 0.8, alpha = 0.6) +
  facet_wrap(~Phenotype_label, scales = "free") +
  geom_text_repel(data = top_labels, aes(label = gene),
                  size = 2.5, max.overlaps = 30, fontface = "bold",
                  color = "black", segment.color = "grey50") +
  geom_hline(yintercept = -log10(0.05 / nrow(df %>% filter(Phenotype == "BMI"))),
             linetype = "dashed", color = "red", linewidth = 0.3) +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(face = "bold", size = 12),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  xlab("Reverse Observational Effect Size") +
  ylab(expression(-log[10](P))) +
  ggtitle("Reverse PEWAS: Phenotype-to-Protein Associations (Normoglycemic)")

ggsave(file.path(outdir, "Fig_Reverse_PEWAS_Summary.pdf"), p, width = 14, height = 5, dpi = 300)
ggsave(file.path(outdir, "Fig_Reverse_PEWAS_Summary.png"), p, width = 14, height = 5, dpi = 300)
cat("Saved: Fig_Reverse_PEWAS_Summary.pdf/png\n")
