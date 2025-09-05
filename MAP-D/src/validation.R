# =============================================================================
# MAP-D External Validation Module
# Metabolic Atlas of Progression to Diabetes (MAP-D)
# =============================================================================

# Source required files
source("config.R")
source("src/utils.R")

#' Load and process STEP1 validation data
#' 
#' @param step1_file Path to STEP1 Excel file
#' @param sheet_name Name of the sheet containing validation data
#' @return Data frame with STEP1 validation results
load_step1_validation <- function(step1_file = STEP1_VALIDATION_FILE, 
                                 sheet_name = "S2_tx_STEP1") {
  
  log_analysis_step("Loading STEP1 Validation Data")
  
  # Load required packages
  load_required_packages(c("readxl", "dplyr"))
  
  # Check if file exists
  if (!check_file_exists(step1_file, "STEP1 validation file")) {
    stop("STEP1 validation file not found")
  }
  
  cat("Loading STEP1 data from:", step1_file, "\n")
  cat("Sheet:", sheet_name, "\n")
  
  # Read validation data
  step1_data <- read_excel(step1_file, sheet = sheet_name)
  
  cat("STEP1 data dimensions:", nrow(step1_data), "rows,", ncol(step1_data), "columns\n")
  
  # Standardize column names
  if ("EntrezGeneSymbol" %in% names(step1_data)) {
    step1_data <- step1_data %>% rename(gene_symbol = EntrezGeneSymbol)
  }
  
  return(step1_data)
}

#' Load and process STEP2 validation data
#' 
#' @param step2_file Path to STEP2 Excel file
#' @param sheet_name Name of the sheet containing validation data
#' @return Data frame with STEP2 validation results
load_step2_validation <- function(step2_file = STEP1_VALIDATION_FILE, 
                                 sheet_name = "S3_tx_STEP2") {
  
  log_analysis_step("Loading STEP2 Validation Data")
  
  # Load required packages
  load_required_packages(c("readxl", "dplyr"))
  
  # Check if file exists
  if (!check_file_exists(step2_file, "STEP2 validation file")) {
    stop("STEP2 validation file not found")
  }
  
  cat("Loading STEP2 data from:", step2_file, "\n")
  cat("Sheet:", sheet_name, "\n")
  
  # Read validation data
  step2_data <- read_excel(step2_file, sheet = sheet_name)
  
  cat("STEP2 data dimensions:", nrow(step2_data), "rows,", ncol(step2_data), "columns\n")
  
  # Standardize column names
  if ("EntrezGeneSymbol" %in% names(step2_data)) {
    step2_data <- step2_data %>% rename(gene_symbol = EntrezGeneSymbol)
  }
  
  return(step2_data)
}

#' Validate PWAS results against external datasets
#' 
#' @param pwas_results MAP-D PWAS results data frame
#' @param validation_data External validation data frame
#' @param significance_threshold Significance threshold for validation
#' @return List with validation results
validate_pwas_results <- function(pwas_results, validation_data, 
                                 significance_threshold = BONFERRONI_THRESHOLD) {
  
  log_analysis_step("Validating PWAS Results Against External Data")
  
  # Extract gene symbols from protein names in PWAS results
  pwas_with_genes <- pwas_results %>%
    mutate(gene_symbol = str_extract(Protein, "^[^;]+")) %>%
    filter(!is.na(gene_symbol))
  
  # Filter for significant results
  significant_pwas <- pwas_with_genes %>%
    filter(Bonferroni < significance_threshold)
  
  cat("Significant PWAS results:", nrow(significant_pwas), "\n")
  cat("Validation dataset size:", nrow(validation_data), "\n")
  
  # Merge with validation data
  merged_results <- inner_join(
    significant_pwas, 
    validation_data, 
    by = "gene_symbol"
  )
  
  cat("Overlapping proteins:", nrow(merged_results), "\n")
  
  if (nrow(merged_results) == 0) {
    warning("No overlapping proteins found between PWAS and validation data")
    return(list(
      n_significant_pwas = nrow(significant_pwas),
      n_validation = nrow(validation_data),
      n_overlap = 0,
      validated_results = data.frame(),
      validation_stats = data.frame()
    ))
  }
  
  # Check for directional consistency
  if (all(c("estimate", "effect_size", "qvalue") %in% names(merged_results))) {
    
    validated_results <- merged_results %>%
      filter(
        qvalue < 0.05,  # Significant in validation dataset
        ((estimate > 0 & effect_size > 0) | (estimate < 0 & effect_size < 0))  # Same direction
      )
    
    cat("Directionally validated results:", nrow(validated_results), "\n")
    
  } else {
    warning("Required columns not found for directional validation")
    validated_results <- merged_results
  }
  
  # Calculate validation statistics by phenotype and subgroup
  if ("Phenotype" %in% names(merged_results) && "Subgroup" %in% names(merged_results)) {
    
    validation_stats <- merged_results %>%
      group_by(Phenotype, Subgroup) %>%
      summarise(
        n_total_significant = n(),
        n_validated = sum(qvalue < 0.05 & 
                         ((estimate > 0 & effect_size > 0) | 
                          (estimate < 0 & effect_size < 0)), na.rm = TRUE),
        validation_rate = n_validated / n_total_significant,
        median_pwas_p = median(p.value, na.rm = TRUE),
        median_validation_q = median(qvalue, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(validation_rate))
    
  } else {
    validation_stats <- data.frame(
      total_overlap = nrow(merged_results),
      validated = nrow(validated_results),
      validation_rate = nrow(validated_results) / nrow(merged_results)
    )
  }
  
  return(list(
    n_significant_pwas = nrow(significant_pwas),
    n_validation = nrow(validation_data),
    n_overlap = nrow(merged_results),
    merged_results = merged_results,
    validated_results = validated_results,
    validation_stats = validation_stats
  ))
}

#' Run comprehensive external validation
#' 
#' @param pwas_results MAP-D PWAS results
#' @param output_dir Directory to save validation results
#' @return List with validation results for both STEP1 and STEP2
run_external_validation <- function(pwas_results, output_dir = RESULTS_DIR) {
  
  log_analysis_step("Comprehensive External Validation", 
                   "Validating against STEP1 and STEP2 datasets")
  
  # Create validation output directory
  validation_dir <- file.path(output_dir, "validation")
  dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)
  
  validation_results <- list()
  
  # STEP1 validation
  tryCatch({
    step1_data <- load_step1_validation()
    step1_validation <- validate_pwas_results(pwas_results, step1_data)
    validation_results$step1 <- step1_validation
    
    # Save STEP1 results
    if (nrow(step1_validation$validated_results) > 0) {
      write_csv(step1_validation$validated_results, 
               file.path(validation_dir, "step1_validated_results.csv"))
    }
    
    if (nrow(step1_validation$validation_stats) > 0) {
      write_csv(step1_validation$validation_stats, 
               file.path(validation_dir, "step1_validation_stats.csv"))
    }
    
  }, error = function(e) {
    warning("STEP1 validation failed: ", e$message)
    validation_results$step1 <- NULL
  })
  
  # STEP2 validation
  tryCatch({
    step2_data <- load_step2_validation()
    step2_validation <- validate_pwas_results(pwas_results, step2_data)
    validation_results$step2 <- step2_validation
    
    # Save STEP2 results
    if (nrow(step2_validation$validated_results) > 0) {
      write_csv(step2_validation$validated_results, 
               file.path(validation_dir, "step2_validated_results.csv"))
    }
    
    if (nrow(step2_validation$validation_stats) > 0) {
      write_csv(step2_validation$validation_stats, 
               file.path(validation_dir, "step2_validation_stats.csv"))
    }
    
  }, error = function(e) {
    warning("STEP2 validation failed: ", e$message)
    validation_results$step2 <- NULL
  })
  
  # Create validation summary
  create_validation_summary(validation_results, validation_dir)
  
  # Save complete validation results
  saveRDS(validation_results, file.path(validation_dir, "complete_validation_results.rds"))
  
  return(validation_results)
}

#' Create validation summary report
#' 
#' @param validation_results Results from run_external_validation
#' @param output_dir Output directory
#' @return NULL (creates summary file)
create_validation_summary <- function(validation_results, output_dir) {
  
  summary_file <- file.path(output_dir, "validation_summary.txt")
  
  sink(summary_file)
  
  cat("=" , rep("=", 50), "=\n", sep = "")
  cat("MAP-D External Validation Summary\n")
  cat("Generated on:", as.character(Sys.time()), "\n")
  cat("=" , rep("=", 50), "=\n", sep = "")
  
  # STEP1 validation summary
  if (!is.null(validation_results$step1)) {
    cat("\nSTEP1 VALIDATION:\n")
    cat("-", rep("-", 30), "-\n", sep = "")
    
    step1 <- validation_results$step1
    cat("Significant PWAS results:", step1$n_significant_pwas, "\n")
    cat("STEP1 validation proteins:", step1$n_validation, "\n")
    cat("Overlapping proteins:", step1$n_overlap, "\n")
    cat("Validated results:", nrow(step1$validated_results), "\n")
    
    if (step1$n_overlap > 0) {
      validation_rate <- nrow(step1$validated_results) / step1$n_overlap
      cat("Overall validation rate:", sprintf("%.1f%%", validation_rate * 100), "\n")
    }
    
    # Top validated hits
    if (nrow(step1$validated_results) > 0) {
      top_hits <- step1$validated_results %>%
        arrange(Bonferroni) %>%
        head(5)
      
      cat("\nTop validated associations:\n")
      for (i in seq_len(nrow(top_hits))) {
        cat("  ", top_hits$gene_symbol[i], " (", top_hits$Phenotype[i], 
            ", FDR=", sprintf("%.2e", top_hits$FDR[i]), ")\n", sep = "")
      }
    }
  } else {
    cat("\nSTEP1 VALIDATION: Not available\n")
  }
  
  # STEP2 validation summary
  if (!is.null(validation_results$step2)) {
    cat("\nSTEP2 VALIDATION:\n")
    cat("-", rep("-", 30), "-\n", sep = "")
    
    step2 <- validation_results$step2
    cat("Significant PWAS results:", step2$n_significant_pwas, "\n")
    cat("STEP2 validation proteins:", step2$n_validation, "\n")
    cat("Overlapping proteins:", step2$n_overlap, "\n")
    cat("Validated results:", nrow(step2$validated_results), "\n")
    
    if (step2$n_overlap > 0) {
      validation_rate <- nrow(step2$validated_results) / step2$n_overlap
      cat("Overall validation rate:", sprintf("%.1f%%", validation_rate * 100), "\n")
    }
    
    # Top validated hits
    if (nrow(step2$validated_results) > 0) {
      top_hits <- step2$validated_results %>%
        arrange(Bonferroni) %>%
        head(5)
      
      cat("\nTop validated associations:\n")
      for (i in seq_len(nrow(top_hits))) {
        cat("  ", top_hits$gene_symbol[i], " (", top_hits$Phenotype[i], 
            ", FDR=", sprintf("%.2e", top_hits$FDR[i]), ")\n", sep = "")
      }
    }
  } else {
    cat("\nSTEP2 VALIDATION: Not available\n")
  }
  
  cat("\n", rep("=", 52), "\n", sep = "")
  cat("Validation analysis completed.\n")
  cat("Detailed results saved in:", output_dir, "\n")
  cat(rep("=", 52), "\n", sep = "")
  
  sink()
  
  cat("Validation summary saved to:", summary_file, "\n")
}

#' Create validation visualization
#' 
#' @param validation_results Results from run_external_validation
#' @param output_file Output file path (optional)
#' @return ggplot object
create_validation_plot <- function(validation_results, output_file = NULL) {
  
  load_required_packages(c("ggplot2", "dplyr"))
  
  # Combine validation statistics
  plot_data <- list()
  
  if (!is.null(validation_results$step1) && 
      nrow(validation_results$step1$validation_stats) > 0) {
    plot_data$step1 <- validation_results$step1$validation_stats %>%
      mutate(dataset = "STEP1")
  }
  
  if (!is.null(validation_results$step2) && 
      nrow(validation_results$step2$validation_stats) > 0) {
    plot_data$step2 <- validation_results$step2$validation_stats %>%
      mutate(dataset = "STEP2")
  }
  
  if (length(plot_data) == 0) {
    warning("No validation statistics available for plotting")
    return(NULL)
  }
  
  combined_data <- bind_rows(plot_data)
  
  # Create validation rate plot
  p <- ggplot(combined_data, aes(x = Phenotype, y = validation_rate, 
                                fill = dataset)) +
    geom_col(position = position_dodge(width = 0.8), alpha = 0.8) +
    facet_wrap(~Subgroup, scales = "free_x") +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
    scale_fill_manual(values = c("STEP1" = "#1f77b4", "STEP2" = "#ff7f0e")) +
    labs(
      title = "External Validation Rates",
      subtitle = "Proportion of significant associations validated in external datasets",
      x = "Phenotype",
      y = "Validation Rate",
      fill = "Validation Dataset"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom"
    )
  
  # Save plot if output file specified
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = 12, height = 8, dpi = PLOT_DPI, units = "in")
    cat("Validation plot saved to:", output_file, "\n")
  }
  
  return(p)
}

cat("External validation module loaded successfully.\n")
