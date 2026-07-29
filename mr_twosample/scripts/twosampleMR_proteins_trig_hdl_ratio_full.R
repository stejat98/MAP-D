# TwoSampleMR Analysis: Validated Proteins -> TRIG_HDL_RATIO (Full Strata)
# This script performs two-sample Mendelian Randomization to test
# causal effects of Step 1 & 2 validated proteins on TRIG_HDL_RATIO

# Set library path
.libPaths("/n/groups/patel/sivateja/.R/library")

# Load required packages
library(TwoSampleMR)
library(readxl)
library(dplyr)
library(data.table)

cat("=== TwoSampleMR: Validated Proteins -> TRIG_HDL_RATIO (Full Strata) ===\n\n")

# Set working directory
setwd("/n/groups/patel/sivateja/regenie_pipeline")

# -----------------------------------------------------------------------------
# 1. Load validated proteins from merged_step1_2_ukb_T2D_coloc.csv
# -----------------------------------------------------------------------------
cat("Step 1: Loading validated proteins...\n")

merged_data <- fread("merged_step1_2_ukb_T2D_coloc.csv")
cat(sprintf("  Total rows in merged file: %d\n", nrow(merged_data)))

unique_proteins <- unique(merged_data$code_UKB)
cat(sprintf("  Unique proteins (code_UKB): %d\n", length(unique_proteins)))

# -----------------------------------------------------------------------------
# 2. Load pQTL data from ST10 sheet (SNP-protein associations)
# -----------------------------------------------------------------------------
cat("\nStep 2: Loading pQTL data from ST10...\n")

st10 <- read_excel("41586_2023_6592_MOESM3_ESM.xlsx", sheet = "ST10", skip = 3)
cat(sprintf("  Total pQTLs in ST10: %d\n", nrow(st10)))

colnames(st10) <- gsub(" ", "_", colnames(st10))
colnames(st10) <- gsub("[()]", "", colnames(st10))
colnames(st10) <- gsub(":", "_", colnames(st10))

# -----------------------------------------------------------------------------
# 3. Map proteins: code_UKB -> Assay_Target in ST10
# -----------------------------------------------------------------------------
cat("\nStep 3: Mapping proteins to pQTLs...\n")

st10$Gene_Symbol <- sapply(strsplit(st10$UKBPPP_ProteinID, ":"), function(x) x[1])

overlap_genes <- intersect(unique_proteins, unique(st10$Gene_Symbol))
cat(sprintf("  Proteins found in ST10 (by gene symbol): %d out of %d\n", 
            length(overlap_genes), length(unique_proteins)))

st10_filtered <- st10 %>%
  filter(Gene_Symbol %in% unique_proteins)
cat(sprintf("  Total pQTLs for validated proteins: %d\n", nrow(st10_filtered)))

proteins_with_pqtls <- unique(st10_filtered$Gene_Symbol)
cat(sprintf("  Proteins with at least one pQTL: %d\n", length(proteins_with_pqtls)))

# -----------------------------------------------------------------------------
# 4. Prepare exposure data for TwoSampleMR
# -----------------------------------------------------------------------------
cat("\nStep 4: Preparing exposure data for TwoSampleMR...\n")

variant_parts <- strsplit(st10_filtered$`Variant_ID_CHROM_GENPOS_hg37_A0_A1_imp_v1`, ":")

st10_filtered$chr <- sapply(variant_parts, function(x) x[1])
st10_filtered$pos_hg37 <- sapply(variant_parts, function(x) x[2])
st10_filtered$other_allele <- sapply(variant_parts, function(x) x[3])
st10_filtered$effect_allele <- sapply(variant_parts, function(x) x[4])

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

exposure_data <- exposure_data %>% filter(!is.na(SNP) & SNP != "")
cat(sprintf("  Exposure data rows (with valid rsIDs): %d\n", nrow(exposure_data)))

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
# 5. Load TRIG_HDL_RATIO GWAS outcome data (Full Strata)
# -----------------------------------------------------------------------------
cat("\nStep 5: Loading TRIG_HDL_RATIO GWAS outcome data (Full Strata)...\n")

outcome_gwas <- fread(cmd = "zcat results/TRIG_HDL_RATIO_full_all_chr.regenie.gz")
cat(sprintf("  Total SNPs in TRIG_HDL_RATIO GWAS: %d\n", nrow(outcome_gwas)))

snps_needed <- unique(exposure_dat$SNP)
cat(sprintf("  SNPs needed from outcome: %d\n", length(snps_needed)))

outcome_filtered <- outcome_gwas %>% filter(ID %in% snps_needed)
cat(sprintf("  SNPs found in TRIG_HDL_RATIO GWAS: %d\n", nrow(outcome_filtered)))

outcome_data <- data.frame(
  SNP = outcome_filtered$ID,
  chr = outcome_filtered$CHROM,
  pos = outcome_filtered$GENPOS,
  effect_allele = outcome_filtered$ALLELE1,
  other_allele = outcome_filtered$ALLELE0,
  eaf = outcome_filtered$A1FREQ,
  beta = outcome_filtered$BETA,
  se = outcome_filtered$SE,
  pval = 10^(-outcome_filtered$LOG10P),
  outcome = "TRIG_HDL_RATIO_full_strata",
  id.outcome = "TRIG_HDL_RATIO_full_strata",
  samplesize = outcome_filtered$N,
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
  action = 2
)

cat(sprintf("  Harmonized SNPs: %d\n", nrow(dat_harmonized)))
cat(sprintf("  Unique exposures with harmonized data: %d\n", 
            length(unique(dat_harmonized$id.exposure))))

# -----------------------------------------------------------------------------
# 7. Perform MR analysis
# -----------------------------------------------------------------------------
cat("\nStep 7: Performing MR analysis...\n")

mr_results <- mr(dat_harmonized, method_list = c(
  "mr_ivw",
  "mr_egger_regression",
  "mr_weighted_median",
  "mr_weighted_mode"
))

cat(sprintf("  MR results: %d rows\n", nrow(mr_results)))

mr_results$OR <- exp(mr_results$b)
mr_results$OR_lci95 <- exp(mr_results$b - 1.96 * mr_results$se)
mr_results$OR_uci95 <- exp(mr_results$b + 1.96 * mr_results$se)

# -----------------------------------------------------------------------------
# 8. Sensitivity analyses
# -----------------------------------------------------------------------------
cat("\nStep 8: Running sensitivity analyses...\n")

het_results <- mr_heterogeneity(dat_harmonized)
cat(sprintf("  Heterogeneity tests: %d\n", nrow(het_results)))

pleio_results <- mr_pleiotropy_test(dat_harmonized)
cat(sprintf("  Pleiotropy tests: %d\n", nrow(pleio_results)))

loo_results <- tryCatch({
  mr_leaveoneout(dat_harmonized)
}, error = function(e) {
  cat("  Note: Leave-one-out analysis skipped\n")
  NULL
})

# -----------------------------------------------------------------------------
# 9. Save results
# -----------------------------------------------------------------------------
cat("\nStep 9: Saving results...\n")

output_dir <- "results/twosampleMR"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(mr_results, file.path(output_dir, "MR_proteins_TRIG_HDL_RATIO_full.csv"), row.names = FALSE)
cat(sprintf("  Saved: %s\n", file.path(output_dir, "MR_proteins_TRIG_HDL_RATIO_full.csv")))

write.csv(het_results, file.path(output_dir, "heterogeneity_proteins_TRIG_HDL_RATIO_full.csv"), row.names = FALSE)
cat(sprintf("  Saved: %s\n", file.path(output_dir, "heterogeneity_proteins_TRIG_HDL_RATIO_full.csv")))

write.csv(pleio_results, file.path(output_dir, "pleiotropy_proteins_TRIG_HDL_RATIO_full.csv"), row.names = FALSE)
cat(sprintf("  Saved: %s\n", file.path(output_dir, "pleiotropy_proteins_TRIG_HDL_RATIO_full.csv")))

write.csv(dat_harmonized, file.path(output_dir, "harmonized_data_proteins_TRIG_HDL_RATIO_full.csv"), row.names = FALSE)
cat(sprintf("  Saved: %s\n", file.path(output_dir, "harmonized_data_proteins_TRIG_HDL_RATIO_full.csv")))

if (!is.null(loo_results)) {
  write.csv(loo_results, file.path(output_dir, "leaveoneout_proteins_TRIG_HDL_RATIO_full.csv"), row.names = FALSE)
  cat(sprintf("  Saved: %s\n", file.path(output_dir, "leaveoneout_proteins_TRIG_HDL_RATIO_full.csv")))
}

# -----------------------------------------------------------------------------
# 10. Summary statistics
# -----------------------------------------------------------------------------
cat("\n=== Summary ===\n")

sig_results <- mr_results %>%
  filter(pval < 0.05) %>%
  group_by(method) %>%
  summarise(
    n_significant = n(),
    proteins = paste(unique(exposure), collapse = ", ")
  )

cat("\nSignificant results (p < 0.05) by method:\n")
print(sig_results)

cat("\nTop 20 results (IVW method):\n")
ivw_results <- mr_results %>%
  filter(method == "Inverse variance weighted") %>%
  arrange(pval) %>%
  head(20) %>%
  select(exposure, b, se, pval, OR, OR_lci95, OR_uci95, nsnp)

print(ivw_results)

cat("\n=== Analysis Complete ===\n")
