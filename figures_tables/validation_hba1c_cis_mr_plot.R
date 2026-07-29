.libPaths(c("/n/groups/patel/sivateja/R_libs", .libPaths()))
library(dplyr)
library(ggplot2)
library(ggrepel)

step1 = read.delim('/n/groups/patel/sivateja/STEP1_merged_results.csv', sep=',')
step2 = read.delim('/n/groups/patel/sivateja/STEP2_merged_results.csv', sep=',')

# Load full cis-MR results for HbA1c (all proteins, full cohort)
cis_mr = read.csv('/n/groups/patel/sivateja/regenie_pipeline/results/twosampleMR/supplemental_tables/Table_MR_Cis_pQTL_HbA1c_Full.csv')
# Get nominally significant (p < 0.05) cis-MR hit protein names for HbA1c
cis_mr_hits <- cis_mr %>%
  filter(P.0.05 == "Yes") %>%
  pull(Protein) %>%
  unique()

cat("Cis-MR nominal hits for HbA1c (Full):\n")
print(cis_mr_hits)

##### Build merged data (same as original code)
step1val_ukb = step1 %>% dplyr::select(Phenotype, Subgroup, EntrezGeneID, name, estimate, Bonferroni, Exposure, code)
step2val_ukb = step2 %>% dplyr::select(Phenotype, Subgroup, EntrezGeneID, name, estimate, Bonferroni, Exposure, code)
colnames(step1val_ukb) = c('Phenotype','Subgroup','EntrezGeneID','Protein_UKB','Estimate_UKB','Significance_UKB','Exposure_UKB','code_UKB')
colnames(step2val_ukb) = c('Phenotype','Subgroup','EntrezGeneID','Protein_UKB','Estimate_UKB','Significance_UKB','Exposure_UKB','code_UKB')

ukb = bind_rows(step1val_ukb, step2val_ukb) %>% distinct

step1val_step1 = step1 %>% filter(!is.na(EntrezGeneID)) %>%
  dplyr::select(ANALYTEID, Phenotype, Subgroup, EntrezGeneID, TargetFullName, effect_size, pvalue) %>%
  group_by(ANALYTEID, Phenotype, Subgroup, EntrezGeneID, TargetFullName) %>%
  summarise(effect_size = mean(effect_size), pvalue = mean(pvalue)) %>%
  ungroup

step2val_step2 = step2 %>% filter(!is.na(EntrezGeneID)) %>%
  dplyr::select(ANALYTEID, Phenotype, Subgroup, EntrezGeneID, TargetFullName, effect_size, pvalue) %>%
  group_by(ANALYTEID, Phenotype, Subgroup, EntrezGeneID, TargetFullName) %>%
  summarise(effect_size = mean(effect_size), pvalue = mean(pvalue)) %>%
  ungroup

colnames(step1val_step1) = c('ANALYTEID','Phenotype','Subgroup','EntrezGeneID','Protein_STEP1','Estimate_STEP1','Significance_STEP1')
colnames(step2val_step2) = c('ANALYTEID','Phenotype','Subgroup','EntrezGeneID','Protein_STEP2','Estimate_STEP2','Significance_STEP2')

valdata = full_join(step1val_step1, step2val_step2) %>% distinct
merged_val = full_join(ukb, valdata)
merged_val$ESTIMATE_STEP_MEAN <- rowMeans(merged_val[, c("Estimate_STEP1", "Estimate_STEP2")], na.rm = TRUE)
merged_val_sig2 = merged_val %>% filter(Significance_UKB < 0.05, Significance_STEP2 < 0.05)

##### HbA1c T2D PLOT WITH CIS-MR SHAPE AESTHETIC

# Assign quadrant labels
merged_val_sig2$quadrant <- factor(
  (merged_val_sig2$Estimate_UKB > 0) + 2 * (merged_val_sig2$ESTIMATE_STEP_MEAN > 0),
  levels = c(0, 1, 2, 3),
  labels = c("Decreased UKB / Decreased Step 2", "Increased UKB / Decreased Step 2",
             "Decreased UKB / Increased Step 2", "Increased UKB / Increased Step 2")
)

# Add cis-MR hit flag (match by gene symbol via code_UKB)
merged_val_sig2$cis_MR_hit <- ifelse(merged_val_sig2$code_UKB %in% cis_mr_hits, "Cis-MR Hit", "Not Cis-MR")

# Filter to HbA1c, T2D for the plot
plot_data <- merged_val_sig2 %>% ungroup %>% filter(Subgroup == 'T2D', Phenotype == 'HbA1c')

cat("\nCis-MR hits in plot data:\n")
print(table(plot_data$cis_MR_hit))

# Select top cis-MR hit proteins for labeling
top_proteins <- plot_data %>% filter(cis_MR_hit == "Cis-MR Hit") %>%
  group_by(Phenotype, Subgroup, quadrant) %>%
  arrange(desc(abs(Estimate_UKB) + abs(ESTIMATE_STEP_MEAN))) %>%
  distinct(Protein_UKB, .keep_all = TRUE) %>%
  slice_head(n = 7) %>%
  ungroup()

# Scatter plot: faded circles for non-cis-MR, prominent triangles for cis-MR hits
p <- ggplot(plot_data, aes(x = Estimate_UKB, y = ESTIMATE_STEP_MEAN, color = quadrant)) +
  geom_point(data = plot_data %>% filter(cis_MR_hit == "Not Cis-MR"),
             shape = 16, size = 1.2, alpha = 0.15) +
  geom_point(data = plot_data %>% filter(cis_MR_hit == "Cis-MR Hit"),
             aes(shape = "HbA1c-associated & reverted post-intervention"),
             size = 3.5, alpha = 0.9) +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'red') +
  xlim(-1, 1) + ylim(-1, 1) +
  scale_color_manual(values = c("blue", "darkgreen", "red", "darkorange"), guide = "none") +
  scale_shape_manual(name = NULL, values = c("HbA1c-associated & reverted post-intervention" = 15)) +
  geom_text_repel(alpha = 0.7, data = top_proteins, aes(label = Protein_UKB),
                  max.overlaps = 100, show.legend = FALSE) +
  xlab('Estimate UKB') + ylab('Estimate Step 2') +
  ggtitle('T2D HbA1c Associations vs. Step 2 Results') +
  theme(legend.position = 'none')

ggsave('/n/groups/patel/sivateja/validation_step2_points_HBA1C.pdf', plot = p, width = 8, height = 8)
cat("\nPlot saved to /n/groups/patel/sivateja/validation_step2_points_HBA1C.pdf\n")
