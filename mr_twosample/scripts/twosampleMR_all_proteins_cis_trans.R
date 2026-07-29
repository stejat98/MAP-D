# TwoSampleMR Analysis: ALL Proteins -> ALL Hallmarks
# Separate analysis for cis vs trans pQTLs

# Set library path
.libPaths("/path/to/project/.R/library")

# Load required packages
library(TwoSampleMR)
library(readxl)
library(dplyr)
library(data.table)
library(ggplot2)

cat("=== TwoSampleMR: All Proteins Cis vs Trans Analysis ===\n\n")

# Set working directory
setwd("/path/to/project/regenie_pipeline")

# Create output directory
output_dir <- "results/twosampleMR/cis_trans"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Load validated proteins
# -----------------------------------------------------------------------------
cat("Step 1: Loading validated proteins...\n")

merged_data <- fread("merged_step1_2_ukb_T2D_coloc.csv")
unique_proteins <- unique(merged_data$code_UKB)
cat(sprintf("  Unique validated proteins: %d\n", length(unique_proteins)))

# -----------------------------------------------------------------------------
# 2. Load pQTL data from ST10 sheet with cis/trans info
# -----------------------------------------------------------------------------
cat("\nStep 2: Loading pQTL data from ST10...\n")

st10 <- read_excel("41586_2023_6592_MOESM3_ESM.xlsx", sheet = "ST10", skip = 3)

# Clean column names
colnames(st10) <- gsub(" ", "_", colnames(st10))
colnames(st10) <- gsub("[()]", "", colnames(st10))
colnames(st10) <- gsub(":", "_", colnames(st10))
colnames(st10) <- gsub("/", "_", colnames(st10))

# Extract gene symbol
st10$Gene_Symbol <- sapply(strsplit(st10$UKBPPP_ProteinID, ":"), function(x) x[1])

# Filter for validated proteins
st10_filtered <- st10 %>% filter(Gene_Symbol %in% unique_proteins)
cat(sprintf("  Total pQTLs for validated proteins: %d\n", nrow(st10_filtered)))

# Parse Variant ID
variant_parts <- strsplit(st10_filtered$`Variant_ID_CHROM_GENPOS_hg37_A0_A1_imp_v1`, ":")
st10_filtered$chr <- sapply(variant_parts, function(x) x[1])
st10_filtered$pos_hg37 <- sapply(variant_parts, function(x) x[2])
st10_filtered$other_allele <- sapply(variant_parts, function(x) x[3])
st10_filtered$effect_allele <- sapply(variant_parts, function(x) x[4])

# Summary of cis/trans by protein
cis_trans_summary <- st10_filtered %>%
  group_by(Gene_Symbol) %>%
  summarise(
    n_cis = sum(cis_trans == "cis"),
    n_trans = sum(cis_trans == "trans"),
    n_total = n()
  )
cat("\nCis/Trans summary:\n")
cat(sprintf("  Proteins with cis pQTLs: %d\n", sum(cis_trans_summary$n_cis > 0)))
cat(sprintf("  Proteins with trans pQTLs: %d\n", sum(cis_trans_summary$n_trans > 0)))
cat(sprintf("  Proteins with both: %d\n", sum(cis_trans_summary$n_cis > 0 & cis_trans_summary$n_trans > 0)))

# Save cis/trans summary
write.csv(cis_trans_summary, file.path(output_dir, "protein_cis_trans_summary.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 3. Define outcomes
# -----------------------------------------------------------------------------
outcomes <- list(
  list(phenotype = "BMI", strata = c("full")),
  list(phenotype = "HbA1c", strata = c("full")),
  list(phenotype = "TRIG_HDL_RATIO", strata = c("full"))
)

# -----------------------------------------------------------------------------
# 4. Function to prepare exposure data
# -----------------------------------------------------------------------------
prepare_exposure <- function(pqtl_data, protein_name, qtl_type) {
  if (nrow(pqtl_data) == 0) return(NULL)
  
  exposure_data <- data.frame(
    SNP = pqtl_data$rsID,
    chr = as.numeric(pqtl_data$chr),
    pos = as.numeric(pqtl_data$pos_hg37),
    effect_allele = pqtl_data$effect_allele,
    other_allele = pqtl_data$other_allele,
    eaf = pqtl_data$A1FREQ,
    beta = pqtl_data$BETA,
    se = pqtl_data$SE,
    pval = 10^(-pqtl_data$log10p),
    exposure = protein_name,
    id.exposure = paste0(protein_name, "_", qtl_type),
    stringsAsFactors = FALSE
  )
  
  exposure_data <- exposure_data %>% filter(!is.na(SNP) & SNP != "")
  if (nrow(exposure_data) == 0) return(NULL)
  
  tryCatch({
    format_data(
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
  }, error = function(e) NULL)
}

# -----------------------------------------------------------------------------
# 5. Load and cache outcome data
# -----------------------------------------------------------------------------
cat("\nStep 3: Loading outcome GWAS data...\n")

outcome_cache <- list()

for (outcome_info in outcomes) {
  phenotype <- outcome_info$phenotype
  
  for (stratum in outcome_info$strata) {
    key <- sprintf("%s_%s", phenotype, stratum)
    file_path <- sprintf("results/%s_all_chr.regenie.gz", key)
    
    if (file.exists(file_path)) {
      cat(sprintf("  Loading %s...\n", key))
      gwas <- fread(cmd = sprintf("zcat %s", file_path))
      outcome_cache[[key]] <- gwas
    } else {
      cat(sprintf("  Warning: File not found: %s\n", file_path))
    }
  }
}

# -----------------------------------------------------------------------------
# 6. Run MR for all proteins
# -----------------------------------------------------------------------------
cat("\nStep 4: Running MR analyses...\n")

all_results <- list()
proteins_with_pqtls <- unique(st10_filtered$Gene_Symbol)

total_proteins <- length(proteins_with_pqtls)
counter <- 0

for (protein in proteins_with_pqtls) {
  counter <- counter + 1
  if (counter %% 10 == 0) {
    cat(sprintf("  Processing protein %d/%d: %s\n", counter, total_proteins, protein))
  }
  
  protein_pqtls <- st10_filtered %>% filter(Gene_Symbol == protein)
  
  for (qtl_type in c("cis", "trans")) {
    # Filter by qtl type
    if (qtl_type == "cis") {
      pqtl_subset <- protein_pqtls %>% filter(cis_trans == "cis")
    } else {
      pqtl_subset <- protein_pqtls %>% filter(cis_trans == "trans")
    }
    
    if (nrow(pqtl_subset) == 0) next
    
    # Prepare exposure
    exposure_dat <- prepare_exposure(pqtl_subset, protein, qtl_type)
    if (is.null(exposure_dat) || nrow(exposure_dat) == 0) next
    
    snps_needed <- unique(exposure_dat$SNP)
    
    # Run against each outcome
    for (outcome_info in outcomes) {
      phenotype <- outcome_info$phenotype
      
      for (stratum in outcome_info$strata) {
        outcome_key <- sprintf("%s_%s", phenotype, stratum)
        
        if (!(outcome_key %in% names(outcome_cache))) next
        
        gwas <- outcome_cache[[outcome_key]]
        
        # Filter for needed SNPs
        gwas_filtered <- gwas[ID %in% snps_needed]
        if (nrow(gwas_filtered) == 0) next
        
        # Format outcome
        outcome_data <- data.frame(
          SNP = gwas_filtered$ID,
          effect_allele = gwas_filtered$ALLELE1,
          other_allele = gwas_filtered$ALLELE0,
          eaf = gwas_filtered$A1FREQ,
          beta = gwas_filtered$BETA,
          se = gwas_filtered$SE,
          pval = 10^(-gwas_filtered$LOG10P),
          outcome = outcome_key,
          id.outcome = outcome_key,
          samplesize = gwas_filtered$N,
          stringsAsFactors = FALSE
        )
        
        outcome_dat <- tryCatch({
          format_data(
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
        }, error = function(e) NULL)
        
        if (is.null(outcome_dat) || nrow(outcome_dat) == 0) next
        
        # Harmonize
        dat_harmonized <- tryCatch({
          harmonise_data(
            exposure_dat = exposure_dat,
            outcome_dat = outcome_dat,
            action = 2
          )
        }, error = function(e) NULL)
        
        if (is.null(dat_harmonized) || nrow(dat_harmonized) == 0) next
        
        # Run MR
        if (nrow(dat_harmonized) == 1) {
          mr_res <- tryCatch({
            mr(dat_harmonized, method_list = c("mr_wald_ratio"))
          }, error = function(e) NULL)
        } else {
          mr_res <- tryCatch({
            mr(dat_harmonized, method_list = c(
              "mr_ivw",
              "mr_egger_regression",
              "mr_weighted_median",
              "mr_weighted_mode"
            ))
          }, error = function(e) NULL)
        }
        
        if (!is.null(mr_res) && nrow(mr_res) > 0) {
          mr_res$qtl_type <- qtl_type
          mr_res$phenotype <- phenotype
          mr_res$stratum <- stratum
          mr_res$protein <- protein
          
          result_key <- sprintf("%s_%s_%s_%s", protein, qtl_type, phenotype, stratum)
          all_results[[result_key]] <- mr_res
        }
      }
    }
  }
}

# -----------------------------------------------------------------------------
# 7. Combine and save results
# -----------------------------------------------------------------------------
cat("\nStep 5: Saving results...\n")

if (length(all_results) > 0) {
  combined_results <- do.call(rbind, all_results)
  
  # Add OR columns
  combined_results$OR <- exp(combined_results$b)
  combined_results$OR_lci95 <- exp(combined_results$b - 1.96 * combined_results$se)
  combined_results$OR_uci95 <- exp(combined_results$b + 1.96 * combined_results$se)
  
  # Save full results
  write.csv(combined_results, file.path(output_dir, "MR_all_proteins_cis_trans.csv"), row.names = FALSE)
  cat(sprintf("  Saved: %s\n", file.path(output_dir, "MR_all_proteins_cis_trans.csv")))
  
  # Create summary (IVW/Wald ratio only)
  summary_results <- combined_results %>%
    filter(method %in% c("Inverse variance weighted", "Wald ratio")) %>%
    select(protein, qtl_type, phenotype, stratum, nsnp, b, se, pval, OR, OR_lci95, OR_uci95) %>%
    arrange(phenotype, stratum, protein, qtl_type)
  
  write.csv(summary_results, file.path(output_dir, "MR_all_proteins_cis_trans_IVW_summary.csv"), row.names = FALSE)
  cat(sprintf("  Saved: %s\n", file.path(output_dir, "MR_all_proteins_cis_trans_IVW_summary.csv")))
  
  # Summary statistics
  cat("\n=== Summary ===\n")
  cat(sprintf("  Total MR results: %d\n", nrow(combined_results)))
  cat(sprintf("  Unique proteins: %d\n", length(unique(combined_results$protein))))
  
  sig_cis <- summary_results %>% filter(qtl_type == "cis" & pval < 0.05)
  sig_trans <- summary_results %>% filter(qtl_type == "trans" & pval < 0.05)
  
  cat(sprintf("  Significant cis results (p<0.05): %d\n", nrow(sig_cis)))
  cat(sprintf("  Significant trans results (p<0.05): %d\n", nrow(sig_trans)))
  
  # Save significant results
  write.csv(sig_cis, file.path(output_dir, "MR_significant_cis_p05.csv"), row.names = FALSE)
  write.csv(sig_trans, file.path(output_dir, "MR_significant_trans_p05.csv"), row.names = FALSE)
  
} else {
  cat("  No results to save.\n")
}

cat("\n=== Analysis Complete ===\n")
cat(sprintf("Results saved to: %s\n", output_dir))
