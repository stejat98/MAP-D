# TwoSampleMR Analysis: Proteins -> Hallmarks (CHUNKED VERSION)
# Arguments: qtl_type (cis/trans), chunk_id, total_chunks

# Set library path
.libPaths("/n/groups/patel/sivateja/.R/library")

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript script.R <qtl_type: cis|trans> <chunk_id> <total_chunks>")
}

qtl_type <- args[1]  # "cis" or "trans"
chunk_id <- as.integer(args[2])
total_chunks <- as.integer(args[3])

if (!qtl_type %in% c("cis", "trans")) {
  stop("qtl_type must be 'cis' or 'trans'")
}

cat(sprintf("=== TwoSampleMR: %s pQTLs, Chunk %d/%d ===\n\n", 
            toupper(qtl_type), chunk_id, total_chunks))

# Load required packages
library(TwoSampleMR)
library(readxl)
library(dplyr)
library(data.table)

# Set working directory
setwd("/n/groups/patel/sivateja/regenie_pipeline")

# Create output directory
output_dir <- "results/twosampleMR/cis_trans_chunked"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Load validated proteins
# -----------------------------------------------------------------------------
cat("Step 1: Loading validated proteins...\n")

merged_data <- fread("merged_step1_2_ukb_T2D_coloc.csv")
unique_proteins <- sort(unique(merged_data$code_UKB))
cat(sprintf("  Total validated proteins: %d\n", length(unique_proteins)))

# -----------------------------------------------------------------------------
# 2. Load pQTL data from ST10 sheet
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

# Filter for validated proteins AND qtl_type
st10_filtered <- st10 %>% 
  filter(Gene_Symbol %in% unique_proteins) %>%
  filter(cis_trans == qtl_type)

cat(sprintf("  Total %s pQTLs for validated proteins: %d\n", qtl_type, nrow(st10_filtered)))

# Parse Variant ID
variant_parts <- strsplit(st10_filtered$`Variant_ID_CHROM_GENPOS_hg37_A0_A1_imp_v1`, ":")
st10_filtered$chr <- sapply(variant_parts, function(x) x[1])
st10_filtered$pos_hg37 <- sapply(variant_parts, function(x) x[2])
st10_filtered$other_allele <- sapply(variant_parts, function(x) x[3])
st10_filtered$effect_allele <- sapply(variant_parts, function(x) x[4])

# Get proteins with this qtl_type
proteins_with_qtl <- unique(st10_filtered$Gene_Symbol)
cat(sprintf("  Proteins with %s pQTLs: %d\n", qtl_type, length(proteins_with_qtl)))

# -----------------------------------------------------------------------------
# 3. Chunk proteins
# -----------------------------------------------------------------------------
cat("\nStep 3: Chunking proteins...\n")

# Split proteins into chunks
chunk_size <- ceiling(length(proteins_with_qtl) / total_chunks)
start_idx <- (chunk_id - 1) * chunk_size + 1
end_idx <- min(chunk_id * chunk_size, length(proteins_with_qtl))

if (start_idx > length(proteins_with_qtl)) {
  cat("  No proteins in this chunk. Exiting.\n")
  quit(save = "no", status = 0)
}

chunk_proteins <- proteins_with_qtl[start_idx:end_idx]
cat(sprintf("  This chunk: proteins %d-%d (%d proteins)\n", 
            start_idx, end_idx, length(chunk_proteins)))
cat(sprintf("  Proteins: %s\n", paste(head(chunk_proteins, 5), collapse = ", ")))
if (length(chunk_proteins) > 5) cat("  ...\n")

# Filter pQTLs for this chunk
st10_chunk <- st10_filtered %>% filter(Gene_Symbol %in% chunk_proteins)
cat(sprintf("  pQTLs in this chunk: %d\n", nrow(st10_chunk)))

# -----------------------------------------------------------------------------
# 4. Define outcomes
# -----------------------------------------------------------------------------
outcomes <- list(
  list(phenotype = "BMI", strata = c("full")),
  list(phenotype = "HbA1c", strata = c("full")),
  list(phenotype = "TRIG_HDL_RATIO", strata = c("full"))
)

# -----------------------------------------------------------------------------
# 5. Function to prepare exposure data
# -----------------------------------------------------------------------------
prepare_exposure <- function(pqtl_data, protein_name) {
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
# 6. Load and cache outcome data
# -----------------------------------------------------------------------------
cat("\nStep 4: Loading outcome GWAS data...\n")

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
# 7. Run MR for chunk proteins
# -----------------------------------------------------------------------------
cat("\nStep 5: Running MR analyses...\n")

all_results <- list()
total_proteins <- length(chunk_proteins)
counter <- 0

for (protein in chunk_proteins) {
  counter <- counter + 1
  cat(sprintf("  [%d/%d] %s\n", counter, total_proteins, protein))
  
  protein_pqtls <- st10_chunk %>% filter(Gene_Symbol == protein)
  
  if (nrow(protein_pqtls) == 0) next
  
  # Prepare exposure
  exposure_dat <- prepare_exposure(protein_pqtls, protein)
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
        suppressMessages(harmonise_data(
          exposure_dat = exposure_dat,
          outcome_dat = outcome_dat,
          action = 2
        ))
      }, error = function(e) NULL)
      
      if (is.null(dat_harmonized) || nrow(dat_harmonized) == 0) next
      
      # Run MR
      if (nrow(dat_harmonized) == 1) {
        mr_res <- tryCatch({
          suppressMessages(mr(dat_harmonized, method_list = c("mr_wald_ratio")))
        }, error = function(e) NULL)
      } else {
        mr_res <- tryCatch({
          suppressMessages(mr(dat_harmonized, method_list = c(
            "mr_ivw",
            "mr_egger_regression",
            "mr_weighted_median",
            "mr_weighted_mode"
          )))
        }, error = function(e) NULL)
      }
      
      if (!is.null(mr_res) && nrow(mr_res) > 0) {
        mr_res$qtl_type <- qtl_type
        mr_res$phenotype <- phenotype
        mr_res$stratum <- stratum
        mr_res$protein <- protein
        
        result_key <- sprintf("%s_%s_%s", protein, phenotype, stratum)
        all_results[[result_key]] <- mr_res
      }
    }
  }
}

# -----------------------------------------------------------------------------
# 8. Save results
# -----------------------------------------------------------------------------
cat("\nStep 6: Saving results...\n")

if (length(all_results) > 0) {
  combined_results <- do.call(rbind, all_results)
  
  # Add OR columns
  combined_results$OR <- exp(combined_results$b)
  combined_results$OR_lci95 <- exp(combined_results$b - 1.96 * combined_results$se)
  combined_results$OR_uci95 <- exp(combined_results$b + 1.96 * combined_results$se)
  
  # Save chunk results
  output_file <- sprintf("MR_%s_chunk%d_of_%d.csv", qtl_type, chunk_id, total_chunks)
  write.csv(combined_results, file.path(output_dir, output_file), row.names = FALSE)
  cat(sprintf("  Saved: %s (%d results)\n", output_file, nrow(combined_results)))
} else {
  cat("  No results to save.\n")
}

cat("\n=== Chunk Complete ===\n")
