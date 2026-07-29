# TwoSampleMR Analysis: LEP (Leptin) -> BMI & HbA1c
# Separate analysis for cis vs trans pQTLs with forest plots

# Set library path
.libPaths("/n/groups/patel/sivateja/.R/library")

# Load required packages
library(TwoSampleMR)
library(readxl)
library(dplyr)
library(data.table)
library(ggplot2)
library(gridExtra)

cat("=== TwoSampleMR: LEP (Leptin) Cis vs Trans Analysis ===\n\n")

# Set working directory
setwd("/n/groups/patel/sivateja/regenie_pipeline")

# Create output directory
output_dir <- "results/twosampleMR/leptin_cis_trans"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Load pQTL data from ST10 sheet
# -----------------------------------------------------------------------------
cat("Step 1: Loading pQTL data from ST10...\n")

st10 <- read_excel("41586_2023_6592_MOESM3_ESM.xlsx", sheet = "ST10", skip = 3)

# Clean column names
colnames(st10) <- gsub(" ", "_", colnames(st10))
colnames(st10) <- gsub("[()]", "", colnames(st10))
colnames(st10) <- gsub(":", "_", colnames(st10))
colnames(st10) <- gsub("/", "_", colnames(st10))

# Extract gene symbol
st10$Gene_Symbol <- sapply(strsplit(st10$UKBPPP_ProteinID, ":"), function(x) x[1])

# Filter for LEP only
lep_pqtls <- st10 %>% filter(Gene_Symbol == "LEP")
cat(sprintf("  Total LEP pQTLs: %d\n", nrow(lep_pqtls)))
cat(sprintf("  Cis pQTLs: %d\n", sum(lep_pqtls$cis_trans == "cis")))
cat(sprintf("  Trans pQTLs: %d\n", sum(lep_pqtls$cis_trans == "trans")))

# Parse Variant ID
variant_parts <- strsplit(lep_pqtls$`Variant_ID_CHROM_GENPOS_hg37_A0_A1_imp_v1`, ":")
lep_pqtls$chr <- sapply(variant_parts, function(x) x[1])
lep_pqtls$pos_hg37 <- sapply(variant_parts, function(x) x[2])
lep_pqtls$other_allele <- sapply(variant_parts, function(x) x[3])
lep_pqtls$effect_allele <- sapply(variant_parts, function(x) x[4])

# -----------------------------------------------------------------------------
# 2. Define outcomes to test
# -----------------------------------------------------------------------------
outcomes <- list(
  list(name = "BMI", strata = c("full")),
  list(name = "HbA1c", strata = c("full"))
)

# -----------------------------------------------------------------------------
# 3. Function to prepare exposure data
# -----------------------------------------------------------------------------
prepare_exposure <- function(pqtl_data, qtl_type) {
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
    exposure = paste0("LEP_", qtl_type),
    id.exposure = paste0("LEP_", qtl_type),
    stringsAsFactors = FALSE
  )
  
  exposure_data <- exposure_data %>% filter(!is.na(SNP) & SNP != "")
  
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
}

# -----------------------------------------------------------------------------
# 4. Function to load outcome data
# -----------------------------------------------------------------------------
load_outcome <- function(phenotype, stratum) {
  file_path <- sprintf("results/%s_%s_all_chr.regenie.gz", phenotype, stratum)
  
  if (!file.exists(file_path)) {
    cat(sprintf("  Warning: File not found: %s\n", file_path))
    return(NULL)
  }
  
  gwas <- fread(cmd = sprintf("zcat %s", file_path))
  
  outcome_name <- sprintf("%s_%s", phenotype, stratum)
  
  outcome_data <- data.frame(
    SNP = gwas$ID,
    chr = gwas$CHROM,
    pos = gwas$GENPOS,
    effect_allele = gwas$ALLELE1,
    other_allele = gwas$ALLELE0,
    eaf = gwas$A1FREQ,
    beta = gwas$BETA,
    se = gwas$SE,
    pval = 10^(-gwas$LOG10P),
    outcome = outcome_name,
    id.outcome = outcome_name,
    samplesize = gwas$N,
    stringsAsFactors = FALSE
  )
  
  return(outcome_data)
}

# -----------------------------------------------------------------------------
# 5. Run MR for each combination
# -----------------------------------------------------------------------------
all_results <- list()
all_harmonized <- list()

for (qtl_type in c("cis", "trans", "all")) {
  cat(sprintf("\n=== Processing %s pQTLs ===\n", toupper(qtl_type)))
  
  # Filter pQTLs
  if (qtl_type == "cis") {
    pqtl_subset <- lep_pqtls %>% filter(cis_trans == "cis")
  } else if (qtl_type == "trans") {
    pqtl_subset <- lep_pqtls %>% filter(cis_trans == "trans")
  } else {
    pqtl_subset <- lep_pqtls
  }
  
  if (nrow(pqtl_subset) == 0) {
    cat(sprintf("  No %s pQTLs available, skipping...\n", qtl_type))
    next
  }
  
  cat(sprintf("  Using %d SNPs\n", nrow(pqtl_subset)))
  
  # Prepare exposure
  exposure_dat <- prepare_exposure(pqtl_subset, qtl_type)
  snps_needed <- unique(exposure_dat$SNP)
  
  for (outcome_info in outcomes) {
    phenotype <- outcome_info$name
    
    for (stratum in outcome_info$strata) {
      cat(sprintf("  Processing %s_%s...\n", phenotype, stratum))
      
      # Load outcome
      outcome_raw <- load_outcome(phenotype, stratum)
      if (is.null(outcome_raw)) next
      
      # Filter for needed SNPs
      outcome_filtered <- outcome_raw %>% filter(SNP %in% snps_needed)
      
      if (nrow(outcome_filtered) == 0) {
        cat(sprintf("    No overlapping SNPs for %s_%s\n", phenotype, stratum))
        next
      }
      
      # Format outcome
      outcome_dat <- format_data(
        outcome_filtered,
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
      
      # Harmonize
      dat_harmonized <- harmonise_data(
        exposure_dat = exposure_dat,
        outcome_dat = outcome_dat,
        action = 2
      )
      
      if (nrow(dat_harmonized) == 0) {
        cat(sprintf("    No harmonized SNPs for %s_%s\n", phenotype, stratum))
        next
      }
      
      cat(sprintf("    Harmonized SNPs: %d\n", nrow(dat_harmonized)))
      
      # Store harmonized data
      key <- sprintf("%s_%s_%s", qtl_type, phenotype, stratum)
      all_harmonized[[key]] <- dat_harmonized
      
      # Run MR
      if (nrow(dat_harmonized) == 1) {
        # Single SNP: use Wald ratio
        mr_res <- mr(dat_harmonized, method_list = c("mr_wald_ratio"))
      } else {
        # Multiple SNPs: use all methods
        mr_res <- mr(dat_harmonized, method_list = c(
          "mr_ivw",
          "mr_egger_regression",
          "mr_weighted_median",
          "mr_weighted_mode"
        ))
      }
      
      mr_res$qtl_type <- qtl_type
      mr_res$phenotype <- phenotype
      mr_res$stratum <- stratum
      
      all_results[[key]] <- mr_res
    }
  }
}

# -----------------------------------------------------------------------------
# 6. Combine and save results
# -----------------------------------------------------------------------------
cat("\n=== Saving Results ===\n")

# Combine all MR results
combined_results <- do.call(rbind, all_results)

# Add OR columns
combined_results$OR <- exp(combined_results$b)
combined_results$OR_lci95 <- exp(combined_results$b - 1.96 * combined_results$se)
combined_results$OR_uci95 <- exp(combined_results$b + 1.96 * combined_results$se)

# Save combined results
write.csv(combined_results, file.path(output_dir, "MR_leptin_cis_trans_all.csv"), row.names = FALSE)
cat(sprintf("  Saved: %s\n", file.path(output_dir, "MR_leptin_cis_trans_all.csv")))

# Create summary table
summary_table <- combined_results %>%
  filter(method %in% c("Inverse variance weighted", "Wald ratio")) %>%
  select(qtl_type, phenotype, stratum, nsnp, b, se, pval) %>%
  arrange(phenotype, stratum, qtl_type)

cat("\n=== Summary (IVW/Wald Ratio) ===\n")
print(summary_table)

write.csv(summary_table, file.path(output_dir, "MR_leptin_cis_trans_summary.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 7. Generate Forest Plots
# -----------------------------------------------------------------------------
cat("\n=== Generating Forest Plots ===\n")

# Forest plot function
create_forest_plot <- function(results_df, title) {
  results_df <- results_df %>%
    mutate(
      label = sprintf("%s (%s)", stratum, qtl_type),
      ci_low = b - 1.96 * se,
      ci_high = b + 1.96 * se
    )
  
  p <- ggplot(results_df, aes(x = b, y = reorder(label, b))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.2) +
    geom_point(aes(color = qtl_type, size = nsnp)) +
    scale_color_manual(values = c("cis" = "#E41A1C", "trans" = "#377EB8", "all" = "#4DAF4A")) +
    labs(
      title = title,
      x = "Beta (95% CI)",
      y = "",
      color = "pQTL Type",
      size = "N SNPs"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
  
  return(p)
}

# Generate plots for BMI and HbA1c
for (phenotype in c("BMI", "HbA1c")) {
  plot_data <- combined_results %>%
    filter(phenotype == !!phenotype) %>%
    filter(method %in% c("Inverse variance weighted", "Wald ratio"))
  
  if (nrow(plot_data) > 0) {
    p <- create_forest_plot(plot_data, sprintf("LEP -> %s (Cis vs Trans)", phenotype))
    
    ggsave(
      file.path(output_dir, sprintf("forest_LEP_%s_cis_trans.png", phenotype)),
      p, width = 8, height = 6, dpi = 150
    )
    ggsave(
      file.path(output_dir, sprintf("forest_LEP_%s_cis_trans.pdf", phenotype)),
      p, width = 8, height = 6
    )
    cat(sprintf("  Saved forest plot: forest_LEP_%s_cis_trans.png/pdf\n", phenotype))
  }
}

# -----------------------------------------------------------------------------
# 8. TwoSampleMR native forest plots (for each harmonized dataset)
# -----------------------------------------------------------------------------
cat("\n=== Generating TwoSampleMR Forest Plots ===\n")

for (key in names(all_harmonized)) {
  dat <- all_harmonized[[key]]
  
  if (nrow(dat) > 1) {
    # Generate single SNP analysis for forest plot
    res_single <- mr_singlesnp(dat)
    
    # Forest plot
    p_forest <- mr_forest_plot(res_single)
    
    ggsave(
      file.path(output_dir, sprintf("forest_singlesnp_%s.png", key)),
      p_forest[[1]], width = 10, height = 6, dpi = 150
    )
    cat(sprintf("  Saved: forest_singlesnp_%s.png\n", key))
    
    # Leave-one-out plot
    res_loo <- mr_leaveoneout(dat)
    p_loo <- mr_leaveoneout_plot(res_loo)
    
    ggsave(
      file.path(output_dir, sprintf("leaveoneout_%s.png", key)),
      p_loo[[1]], width = 10, height = 6, dpi = 150
    )
    cat(sprintf("  Saved: leaveoneout_%s.png\n", key))
    
    # Scatter plot
    res_mr <- mr(dat)
    p_scatter <- mr_scatter_plot(res_mr, dat)
    
    ggsave(
      file.path(output_dir, sprintf("scatter_%s.png", key)),
      p_scatter[[1]], width = 8, height = 6, dpi = 150
    )
    cat(sprintf("  Saved: scatter_%s.png\n", key))
  }
}

cat("\n=== Analysis Complete ===\n")
cat(sprintf("Results saved to: %s\n", output_dir))
