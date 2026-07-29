.libPaths(c("/n/groups/patel/sivateja/R_libs", .libPaths()))
library(dplyr)
library(ggplot2)
library(ggrepel)

step1 = read.delim('/n/groups/patel/sivateja/STEP1_merged_results.csv', sep=',')
step2 = read.delim('/n/groups/patel/sivateja/STEP2_merged_results.csv', sep=',')

# Load cis-MR results for BMI
cis_mr = read.csv('/n/groups/patel/sivateja/regenie_pipeline/results/twosampleMR/supplemental_tables/Table_MR_Cis_pQTL_BMI_Full_merged_filtered.csv')
# Get nominally significant (p < 0.05) cis-MR hit protein names for BMI, non T2D
cis_mr_hits <- cis_mr %>%
  filter(Cis_pQTL_Plt0.05 == "Yes", Phenotype == "BMI", Subgroup == "non T2D") %>%
  pull(Protein_UKB) %>%
  unique()

cat("Cis-MR nominal hits for BMI (non T2D):\n")
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
merged_val_sig1 = merged_val %>% filter(Significance_UKB < 0.05, Significance_STEP1 < 0.05)

##### BMI PLOT WITH CIS-MR SHAPE AESTHETIC

# Assign quadrant labels
merged_val_sig1$quadrant <- factor(
  (merged_val_sig1$Estimate_UKB > 0) + 2 * (merged_val_sig1$ESTIMATE_STEP_MEAN > 0),
  levels = c(0, 1, 2, 3),
  labels = c("Decreased UKB / Decreased Step 1", "Increased UKB / Decreased Step 1",
             "Decreased UKB / Increased Step 1", "Increased UKB / Increased Step 1")
)

# Add cis-MR hit flag
merged_val_sig1$cis_MR_hit <- ifelse(merged_val_sig1$Protein_UKB %in% cis_mr_hits, "Cis-MR Hit", "Not Cis-MR")

# Filter to BMI, non T2D for the plot
plot_data <- merged_val_sig1 %>% ungroup %>% filter(Subgroup == 'non T2D', Phenotype == 'BMI')

cat("\nCis-MR hits in plot data:\n")
print(table(plot_data$cis_MR_hit))

# Select top cis-MR hit proteins for labeling
top_proteins <- merged_val_sig1 %>% filter(Subgroup == 'non T2D', Phenotype == 'BMI', cis_MR_hit == "Cis-MR Hit") %>%
  group_by(Phenotype, Subgroup, quadrant) %>%
  arrange(desc(abs(Estimate_UKB) + abs(ESTIMATE_STEP_MEAN))) %>%
  distinct(Protein_UKB, .keep_all = TRUE) %>%
  slice_head(n = 7) %>%
  ungroup()

# Scatter plot with quadrant colors, cis-MR shape, and protein labels
# Plot non-cis-MR points first (faded), then cis-MR on top (prominent)
p <- ggplot(plot_data, aes(x = Estimate_UKB, y = ESTIMATE_STEP_MEAN, color = quadrant)) +
  geom_point(data = plot_data %>% filter(cis_MR_hit == "Not Cis-MR"),
             shape = 16, size = 1.2, alpha = 0.15) +
  geom_point(data = plot_data %>% filter(cis_MR_hit == "Cis-MR Hit"),
             aes(shape = "BMI-associated & reverted post-intervention"),
             size = 3.5, alpha = 0.9) +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'red') +
  xlim(-1, 1) + ylim(-1, 1) +
  scale_color_manual(values = c("blue", "darkgreen", "red", "darkorange"), guide = "none") +
  scale_shape_manual(name = NULL, values = c("BMI-associated & reverted post-intervention" = 15)) +
  geom_text_repel(alpha = 0.7, data = top_proteins, aes(label = Protein_UKB), show.legend = FALSE) +
  xlab('Estimate UKB') + ylab('Estimate Step 1') +
  ggtitle('Normoglycemic BMI Associations vs. Step 1 Results') +
  theme(legend.position = 'none')

ggsave('/n/groups/patel/sivateja/validation_step1_points_BMI.pdf', plot = p, width = 8, height = 8)
cat("\nPlot saved to /n/groups/patel/sivateja/validation_step1_points_BMI.pdf\n")
