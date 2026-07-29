# Combine chunked MR results for cis/trans analysis

.libPaths("/path/to/project/.R/library")
library(dplyr)
library(data.table)

cat("=== Combining Chunked MR Results ===\n\n")

setwd("/path/to/project/regenie_pipeline")
input_dir <- "results/twosampleMR/cis_trans_chunked"
output_dir <- "results/twosampleMR/cis_trans"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Find all chunk files
chunk_files <- list.files(input_dir, pattern = "MR_.*_chunk.*\\.csv$", full.names = TRUE)
cat(sprintf("Found %d chunk files\n", length(chunk_files)))

if (length(chunk_files) == 0) {
  stop("No chunk files found!")
}

# Read and combine all chunks
all_results <- lapply(chunk_files, function(f) {
  cat(sprintf("  Reading: %s\n", basename(f)))
  fread(f)
})

combined <- rbindlist(all_results, fill = TRUE)
cat(sprintf("\nTotal results: %d\n", nrow(combined)))

# Save combined results
write.csv(combined, file.path(output_dir, "MR_all_proteins_cis_trans.csv"), row.names = FALSE)
cat(sprintf("Saved: MR_all_proteins_cis_trans.csv\n"))

# Create IVW/Wald ratio summary
summary_results <- combined %>%
  filter(method %in% c("Inverse variance weighted", "Wald ratio")) %>%
  select(protein, qtl_type, phenotype, stratum, nsnp, b, se, pval, OR, OR_lci95, OR_uci95) %>%
  arrange(phenotype, stratum, protein, qtl_type)

write.csv(summary_results, file.path(output_dir, "MR_all_proteins_cis_trans_IVW_summary.csv"), row.names = FALSE)
cat(sprintf("Saved: MR_all_proteins_cis_trans_IVW_summary.csv\n"))

# Significant results
sig_cis <- summary_results %>% filter(qtl_type == "cis" & pval < 0.05)
sig_trans <- summary_results %>% filter(qtl_type == "trans" & pval < 0.05)

write.csv(sig_cis, file.path(output_dir, "MR_significant_cis_p05.csv"), row.names = FALSE)
write.csv(sig_trans, file.path(output_dir, "MR_significant_trans_p05.csv"), row.names = FALSE)

cat(sprintf("\n=== Summary ===\n"))
cat(sprintf("Unique proteins: %d\n", length(unique(combined$protein))))
cat(sprintf("Cis results: %d\n", sum(combined$qtl_type == "cis")))
cat(sprintf("Trans results: %d\n", sum(combined$qtl_type == "trans")))
cat(sprintf("Significant cis (p<0.05, IVW): %d\n", nrow(sig_cis)))
cat(sprintf("Significant trans (p<0.05, IVW): %d\n", nrow(sig_trans)))

cat("\n=== Done ===\n")
