# TwoSampleMR Analysis: Validated Proteins -> BMI (Full Strata)
# This script performs two-sample Mendelian Randomization to test
# causal effects of Step 1 & 2 validated proteins on BMI

# Set library path
.libPaths("/path/to/project/.R/library")

# Load required packages
library(TwoSampleMR)
library(readxl)
library(dplyr)
library(data.table)

cat("=== TwoSampleMR: Validated Proteins -> BMI (Full Strata) ===\n\n")

# Set working directory
setwd("/path/to/project/regenie_pipeline")

# -----------------------------------------------------------------------------
# 1. Load validated proteins from merged_step1_2_ukb_T2D_coloc.csv
# -----------------------------------------------------------------------------
cat("Step 1: Loading validated proteins...\n")

merged_data <- fread("merged_step1_2_ukb_T2D_coloc.csv")
cat(sprintf("  Total rows in merged file: %d\n", nrow(merged_data)))

# Get unique proteins (code_UKB column)
unique_proteins <- unique(merged_data$code_UKB)
cat(sprintf("  Unique proteins (code_UKB): %d\n", length(unique_proteins)))

# -----------------------------------------------------------------------------
# 2. Load pQTL data from ST10 sheet (SNP-protein associations)
# -----------------------------------------------------------------------------
cat("\nStep 2: Loading pQTL data from ST10...\n")

st10 <- read_excel("41586_2023_6592_MOESM3_ESM.xlsx", sheet = "ST10", skip = 3)
cat(sprintf("  Total pQTLs in ST10: %d\n", nrow(st10)))

# Clean column names for easier access
colnames(st10) <- gsub(" ", "_", colnames(st10))
colnames(st10) <- gsub("[()]", "", colnames(st10))
colnames(st10) <- gsub(":", "_", colnames(st10))

# View column names
cat("  Column names:\n")
print(colnames(st10))

# -----------------------------------------------------------------------------
# 3. Map proteins: code_UKB -> Assay_Target in ST10
# -----------------------------------------------------------------------------
cat("\nStep 3: Mapping proteins to pQTLs...\n")

# Check mapping between code_UKB and Assay_Target
# First, let's see what Assay_Target values look like
cat("  Sample Assay_Target values:\n")
print(head(unique(st10$Assay_Target), 20))

# Filter ST10 for proteins in our validated list
# Note: code_UKB might be gene symbols, Assay_Target might be protein names
# We need to check the mapping

# Let's also check UKBPPP_ProteinID which might contain the gene symbol
cat("\n  Sample UKBPPP_ProteinID values:\n")
print(head(unique(st10$UKBPPP_ProteinID), 20))

# Extract gene symbol from UKBPPP_ProteinID (format: "GENE:seq_id")
st10$Gene_Symbol <- sapply(strsplit(st10$UKBPPP_ProteinID, ":"), function(x) x[1])

# Check overlap
overlap_genes <- intersect(unique_proteins, unique(st10$Gene_Symbol))
cat(sprintf("\n  Proteins found in ST10 (by gene symbol): %d out of %d\n", 
            length(overlap_genes), length(unique_proteins)))

# Filter ST10 for our validated proteins
st10_filtered <- st10 %>%
  filter(Gene_Symbol %in% unique_proteins)
cat(sprintf("  Total pQTLs for validated proteins: %d\n", nrow(st10_filtered)))

# Check how many unique proteins have pQTLs
proteins_with_pqtls <- unique(st10_filtered$Gene_Symbol)
cat(sprintf("  Proteins with at least one pQTL: %d\n", length(proteins_with_pqtls)))

# -----------------------------------------------------------------------------
# 4. Prepare exposure data for TwoSampleMR
# -----------------------------------------------------------------------------
cat("\nStep 4: Preparing exposure data for TwoSampleMR...\n")

# For TwoSampleMR, we need:
# SNP, effect_allele, other_allele, eaf, beta, se, pval, exposure

# Parse Variant ID to extract chr, pos, ref, alt
# Format: "CHROM:GENPOS (hg37):A0:A1:imp:v1"
variant_parts <- strsplit(st10_filtered$`Variant_ID_CHROM_GENPOS_hg37_A0_A1_imp_v1`, ":")

st10_filtered$chr <- sapply(variant_parts, function(x) x[1])
st10_filtered$pos_hg37 <- sapply(variant_parts, function(x) x[2])
st10_filtered$other_allele <- sapply(variant_parts, function(x) x[3])
st10_filtered$effect_allele <- sapply(variant_parts, function(x) x[4])

# Create exposure dataframe
exposure_data <- data.frame(
  SNP = st10_filtered$rsID,
  chr = as.numeric(st10_filtered$chr),
  pos = as.numeric(st10_filtered$pos_hg37),
  effect_allele = st10_filtered$effect_allele,
  other_allele = st10_filtered$other_allele,
  eaf = st10_filtered$A1FREQ,
  beta = st10_filtered$BETA,
  se = st10_filtered$SE,
  pval = 10^(-st10_filtered$`log10p`),
  exposure = st10_filtered$Gene_Symbol,
  id.exposure = st10_filtered$Gene_Symbol,
  stringsAsFactors = FALSE
)

# Remove rows with missing rsIDs
exposure_data <- exposure_data %>% filter(!is.na(SNP) & SNP != "")
cat(sprintf("  Exposure data rows (with valid rsIDs): %d\n", nrow(exposure_data)))

# Format for TwoSampleMR
exposure_dat <- format_data(
  exposure_data,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  phenotype_col = "exposure"
)

cat(sprintf("  Formatted exposure data: %d SNPs\n", nrow(exposure_dat)))

# -----------------------------------------------------------------------------
# 5. Load BMI GWAS outcome data (Full Strata)
# -----------------------------------------------------------------------------
cat("\nStep 5: Loading BMI GWAS outcome data...\n")

# Read BMI GWAS results using gzfile for gzipped data
bmi_gwas <- fread(cmd = "zcat results/BMI_full_all_chr.regenie.gz")
cat(sprintf("  Total SNPs in BMI GWAS: %d\n", nrow(bmi_gwas)))

# Get SNPs we need from exposure
snps_needed <- unique(exposure_dat$SNP)
cat(sprintf("  SNPs needed from outcome: %d\n", length(snps_needed)))

# Filter BMI GWAS for SNPs in exposure
bmi_filtered <- bmi_gwas %>% filter(ID %in% snps_needed)
cat(sprintf("  SNPs found in BMI GWAS: %d\n", nrow(bmi_filtered)))

# Prepare outcome data
outcome_data <- data.frame(
  SNP = bmi_filtered$ID,
  chr = bmi_filtered$CHROM,
  pos = bmi_filtered$GENPOS,
  effect_allele = bmi_filtered$ALLELE1,
  other_allele = bmi_filtered$ALLELE0,
  eaf = bmi_filtered$A1FREQ,
  beta = bmi_filtered$BETA,
  se = bmi_filtered$SE,
  pval = 10^(-bmi_filtered$LOG10P),
  outcome = "BMI_full_strata",
  id.outcome = "BMI_full_strata",
  samplesize = bmi_filtered$N,
  stringsAsFactors = FALSE
)

outcome_dat <- format_data(
  outcome_data,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  phenotype_col = "outcome",
  samplesize_col = "samplesize"
)

cat(sprintf("  Formatted outcome data: %d SNPs\n", nrow(outcome_dat)))

# -----------------------------------------------------------------------------
# 6. Harmonize exposure and outcome data
# -----------------------------------------------------------------------------
cat("\nStep 6: Harmonizing exposure and outcome data...\n")

dat_harmonized <- harmonise_data(
  exposure_dat = exposure_dat,
  outcome_dat = outcome_dat,
  action = 2  # Try to infer positive strand alleles
)

cat(sprintf("  Harmonized SNPs: %d\n", nrow(dat_harmonized)))
cat(sprintf("  Unique exposures with harmonized data: %d\n", 
            length(unique(dat_harmonized$id.exposure))))

# -----------------------------------------------------------------------------
# 7. Perform MR analysis
# -----------------------------------------------------------------------------
cat("\nStep 7: Performing MR analysis...\n")

# Run MR for each protein
mr_results <- mr(dat_harmonized, method_list = c(
  "mr_ivw",
  "mr_egger_regression",
  "mr_weighted_median",
  "mr_weighted_mode"
))

cat(sprintf("  MR results: %d rows\n", nrow(mr_results)))

# Generate OR and CI
mr_results$OR <- exp(mr_results$b)
mr_results$OR_lci95 <- exp(mr_results$b - 1.96 * mr_results$se)
mr_results$OR_uci95 <- exp(mr_results$b + 1.96 * mr_results$se)

# -----------------------------------------------------------------------------
# 8. Sensitivity analyses
# -----------------------------------------------------------------------------
cat("\nStep 8: Running sensitivity analyses...\n")

# Heterogeneity test
het_results <- mr_heterogeneity(dat_harmonized)
cat(sprintf("  Heterogeneity tests: %d\n", nrow(het_results)))

# Pleiotropy test (MR-Egger intercept)
pleio_results <- mr_pleiotropy_test(dat_harmonized)
cat(sprintf("  Pleiotropy tests: %d\n", nrow(pleio_results)))

# Leave-one-out analysis (for proteins with multiple SNPs)
loo_results <- tryCatch({
  mr_leaveoneout(dat_harmonized)
}, error = function(e) {
  cat("  Note: Leave-one-out analysis skipped (possible single SNP exposures)\n")
  NULL
})

# -----------------------------------------------------------------------------
# 9. Save results
# -----------------------------------------------------------------------------
cat("\nStep 9: Saving results...\n")

output_dir <- "results/twosampleMR"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Save main MR results
write.csv(mr_results, file.path(output_dir, "MR_proteins_BMI_full.csv"), row.names = FALSE)
cat(sprintf("  Saved: %s\n", file.path(output_dir, "MR_proteins_BMI_full.csv")))

# Save heterogeneity results
write.csv(het_results, file.path(output_dir, "heterogeneity_proteins_BMI_full.csv"), row.names = FALSE)
cat(sprintf("  Saved: %s\n", file.path(output_dir, "heterogeneity_proteins_BMI_full.csv")))

# Save pleiotropy results
write.csv(pleio_results, file.path(output_dir, "pleiotropy_proteins_BMI_full.csv"), row.names = FALSE)
cat(sprintf("  Saved: %s\n", file.path(output_dir, "pleiotropy_proteins_BMI_full.csv")))

# Save harmonized data for reference
write.csv(dat_harmonized, file.path(output_dir, "harmonized_data_proteins_BMI_full.csv"), row.names = FALSE)
cat(sprintf("  Saved: %s\n", file.path(output_dir, "harmonized_data_proteins_BMI_full.csv")))

# Save leave-one-out if available
if (!is.null(loo_results)) {
  write.csv(loo_results, file.path(output_dir, "leaveoneout_proteins_BMI_full.csv"), row.names = FALSE)
  cat(sprintf("  Saved: %s\n", file.path(output_dir, "leaveoneout_proteins_BMI_full.csv")))
}

# -----------------------------------------------------------------------------
# 10. Summary statistics
# -----------------------------------------------------------------------------
cat("\n=== Summary ===\n")

# Count significant results (p < 0.05) by method
sig_results <- mr_results %>%
  filter(pval < 0.05) %>%
  group_by(method) %>%
  summarise(
    n_significant = n(),
    proteins = paste(unique(exposure), collapse = ", ")
  )

cat("\nSignificant results (p < 0.05) by method:\n")
print(sig_results)

# Top results by IVW p-value
cat("\nTop 20 results (IVW method):\n")
ivw_results <- mr_results %>%
  filter(method == "Inverse variance weighted") %>%
  arrange(pval) %>%
  head(20) %>%
  select(exposure, b, se, pval, OR, OR_lci95, OR_uci95, nsnp)

print(ivw_results)

cat("\n=== Analysis Complete ===\n")
