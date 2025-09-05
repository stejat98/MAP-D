# =============================================================================
# MAP-D PWAS Analysis Module
# Metabolic Atlas of Progression to Diabetes (MAP-D)
# =============================================================================

# Source required files
source("config.R")
source("src/utils.R")

# Source existing PWAS helper function scripts
# Load required packages for PWAS helper functions
load_required_packages(c("tidyverse", "broom"))

if (file.exists(PWAS_LINEAR_FUNCTIONS)) {
  source(PWAS_LINEAR_FUNCTIONS)
  cat("Loaded linear PWAS functions from:", PWAS_LINEAR_FUNCTIONS, "\n")
} else {
  warning("Linear PWAS functions script not found at:", PWAS_LINEAR_FUNCTIONS)
}

if (file.exists(PWAS_LOGISTIC_FUNCTIONS)) {
  source(PWAS_LOGISTIC_FUNCTIONS)
  cat("Loaded logistic PWAS functions from:", PWAS_LOGISTIC_FUNCTIONS, "\n")
} else {
  warning("Logistic PWAS functions script not found at:", PWAS_LOGISTIC_FUNCTIONS)
}

#' Run Protein-Wide Association Study (PWAS) using existing helper functions
#' 
#' @param data Input data frame
#' @param outcome_var Outcome variable name
#' @param protein_vars Vector of protein variable names
#' @param adjustments Vector of adjustment variable names
#' @param analysis_type Type of analysis ("linear" or "logistic")
#' @param stratify_by Variable to stratify analysis by (optional)
#' @param output_prefix Prefix for output files
#' @return Data frame with PWAS results
run_pwas <- function(data, outcome_var, protein_vars, adjustments, 
                    analysis_type = "linear", stratify_by = NULL, 
                    output_prefix = "PWAS") {
  
  log_analysis_step(paste("Running PWAS for", outcome_var), 
                   paste("Analysis type:", analysis_type))
  
  # Validate inputs
  if (!outcome_var %in% names(data)) {
    stop("Outcome variable not found in data: ", outcome_var)
  }
  
  missing_proteins <- setdiff(protein_vars, names(data))
  if (length(missing_proteins) > 0) {
    warning("Some proteins not found in data: ", paste(missing_proteins, collapse = ", "))
    protein_vars <- intersect(protein_vars, names(data))
  }
  
  missing_adjustments <- setdiff(adjustments, names(data))
  if (length(missing_adjustments) > 0) {
    warning("Some adjustment variables not found: ", paste(missing_adjustments, collapse = ", "))
    adjustments <- intersect(adjustments, names(data))
  }
  
  # Run analysis
  if (is.null(stratify_by)) {
    # Single analysis using existing helper functions
    cat("Running PWAS for", length(protein_vars), "proteins using existing helper functions...\n")
    
    # Create output file name
    output_file <- file.path(RESULTS_DIR, paste0(output_prefix, "_", outcome_var, "_", 
                                                analysis_type, "_results"))
    
    if (analysis_type == "linear") {
      # Use existing EWAS function for linear regression (expects adjustments as vector)
      if (exists("EWAS")) {
        EWAS(data = data, 
             depvar = outcome_var, 
             adjustments = adjustments,  # Pass as vector, not string
             exposures = protein_vars, 
             outFileName = output_file)
      } else {
        stop("EWAS function not found. Please check linear PWAS functions script.")
      }
    } else if (analysis_type == "logistic") {
      # Use existing EWASLogistic function for logistic regression (expects adjustments as vector)
      if (exists("EWASLogistic")) {
        EWASLogistic(data = data, 
                    depvar = outcome_var, 
                    adjustments = adjustments,  # Pass as vector, not string
                    exposures = protein_vars, 
                    outFileName = output_file)
      } else {
        stop("EWASLogistic function not found. Please check logistic PWAS functions script.")
      }
    } else {
      stop("Unknown analysis type: ", analysis_type)
    }
    
    # Load and process results
    results_file <- paste0(output_file, ".RDS")
    if (file.exists(results_file)) {
      results <- readRDS(results_file)
      
      # Standardize column names to match expected format
      if ("Exposure" %in% names(results)) {
        results <- results %>% rename(Protein = Exposure)
      }
      
      # Add multiple testing corrections if not already present
      if (!"FDR" %in% names(results) && "p.value" %in% names(results)) {
        results <- results %>%
          mutate(
            FDR = p.adjust(p.value, method = "fdr"),
            Bonferroni = p.adjust(p.value, method = "bonferroni")
          )
      }
      
      # Save as CSV for easier access
      csv_file <- paste0(output_file, ".csv")
      write_csv(results, csv_file)
      
    } else {
      warning("Results file not found: ", results_file)
      return(NULL)
    }
    
  } else {
    # Stratified analysis
    cat("Running stratified PWAS by", stratify_by, "...\n")
    
    strata_levels <- unique(data[[stratify_by]])
    strata_levels <- strata_levels[!is.na(strata_levels)]
    
    results_by_strata <- list()
    
    for (stratum in strata_levels) {
      cat("  Processing stratum:", stratum, "\n")
      
      data_stratum <- data %>% filter(.data[[stratify_by]] == stratum)
      
      if (nrow(data_stratum) < MIN_SAMPLE_SIZE) {
        warning("Insufficient sample size for stratum: ", stratum)
        next
      }
      
      # Create stratum-specific output file
      stratum_output <- file.path(RESULTS_DIR, paste0(output_prefix, "_", outcome_var, "_", 
                                                     analysis_type, "_", stratum, "_results"))
      
      # Run PWAS for this stratum using existing helper functions
      if (analysis_type == "linear") {
        if (exists("EWAS")) {
          EWAS(data = data_stratum, 
               depvar = outcome_var, 
               adjustments = adjustments,  # Pass as vector, not string
               exposures = protein_vars, 
               outFileName = stratum_output)
        }
      } else if (analysis_type == "logistic") {
        if (exists("EWASLogistic")) {
          EWASLogistic(data = data_stratum, 
                      depvar = outcome_var, 
                      adjustments = adjustments,  # Pass as vector, not string
                      exposures = protein_vars, 
                      outFileName = stratum_output)
        }
      }
      
      # Load stratum results
      stratum_results_file <- paste0(stratum_output, ".RDS")
      if (file.exists(stratum_results_file)) {
        stratum_results <- readRDS(stratum_results_file)
        
        # Standardize column names to match expected format
        if ("Exposure" %in% names(stratum_results)) {
          stratum_results <- stratum_results %>% rename(Protein = Exposure)
        }
        
        stratum_results <- stratum_results %>% mutate(Subgroup = stratum)
        results_by_strata[[as.character(stratum)]] <- stratum_results
      }
    }
    
    # Combine results across strata
    if (length(results_by_strata) > 0) {
      results <- bind_rows(results_by_strata)
      
      # Add multiple testing corrections
      if ("p.value" %in% names(results)) {
        results <- results %>%
          mutate(
            FDR = p.adjust(p.value, method = "fdr"),
            Bonferroni = p.adjust(p.value, method = "bonferroni")
          ) %>%
          arrange(p.value)
      }
      
      # Save combined results
      combined_output <- file.path(RESULTS_DIR, paste0(output_prefix, "_", outcome_var, "_", 
                                                      analysis_type, "_combined_results.csv"))
      write_csv(results, combined_output)
      
    } else {
      warning("No valid strata results found")
      return(NULL)
    }
  }
  
  if (exists("results") && nrow(results) > 0) {
    cat("PWAS completed. Results processed successfully.\n")
    if ("FDR" %in% names(results)) {
      cat("Significant associations (FDR < 0.05):", sum(results$FDR < 0.05, na.rm = TRUE), "\n")
    }
    return(results)
  } else {
    warning("No results generated")
    return(NULL)
  }
}

#' Run metabolic phenotype PWAS across glycemic groups
#' 
#' @param data Input data frame with GlycemicStatus
#' @param protein_vars Vector of protein variable names
#' @param adjustments Base adjustment variables
#' @param add_fasting_adjustment Whether to add fasting time adjustment
#' @return List of PWAS results for each phenotype and group
run_metabolic_phenotype_pwas <- function(data, protein_vars, adjustments, 
                                        add_fasting_adjustment = TRUE) {
  
  log_analysis_step("Running Metabolic Phenotype PWAS", 
                   "Analyzing metabolic traits across glycemic groups")
  
  # Add fasting time adjustment if requested
  if (add_fasting_adjustment && FASTING_TIME_VAR %in% names(data)) {
    adjustments <- c(adjustments, FASTING_TIME_VAR)
    cat("Including fasting time adjustment\n")
  }
  
  # Check for required glycemic status variable
  if (!"GlycemicStatus" %in% names(data)) {
    stop("GlycemicStatus variable not found in data")
  }
  
  results_list <- list()
  
  # Run PWAS for each metabolic phenotype
  for (phenotype in METABOLIC_PHENOTYPES) {
    
    if (!phenotype %in% names(data)) {
      warning("Phenotype not found in data: ", phenotype)
      next
    }
    
    cat("\n", rep("-", 50), "\n")
    cat("Analyzing phenotype:", phenotype, "\n")
    cat(rep("-", 50), "\n")
    
    # Run stratified analysis by glycemic status
    phenotype_results <- run_pwas(
      data = data,
      outcome_var = phenotype,
      protein_vars = protein_vars,
      adjustments = adjustments,
      analysis_type = "linear",
      stratify_by = "GlycemicStatus",
      output_prefix = paste0("MAP_D_", phenotype)
    )
    
    results_list[[phenotype]] <- phenotype_results
  }
  
  # Combine all results
  combined_results <- bind_rows(results_list, .id = "phenotype_id")
  
  # Save combined results
  combined_output <- file.path(RESULTS_DIR, "MAP_D_metabolic_phenotypes_combined_results.csv")
  write_csv(combined_results, combined_output)
  
  cat("\nMetabolic phenotype PWAS completed.\n")
  cat("Combined results saved to:", combined_output, "\n")
  
  return(results_list)
}

#' Run complication outcome PWAS
#' 
#' @param data Input data frame
#' @param protein_vars Vector of protein variable names
#' @param adjustments Base adjustment variables
#' @param add_glycemic_adjustment Whether to add glycemic status and related adjustments
#' @return List of PWAS results for each complication
run_complication_pwas <- function(data, protein_vars, adjustments, 
                                 add_glycemic_adjustment = TRUE) {
  
  log_analysis_step("Running Complication Outcome PWAS", 
                   "Analyzing cardiovascular and metabolic complications")
  
  # Add glycemic adjustments if requested
  if (add_glycemic_adjustment) {
    glycemic_adjustments <- c(GLYCEMIC_STATUS_VAR, "BMI", "HbA1c")
    available_adjustments <- intersect(glycemic_adjustments, names(data))
    adjustments <- c(adjustments, available_adjustments)
    cat("Including glycemic adjustments:", paste(available_adjustments, collapse = ", "), "\n")
  }
  
  results_list <- list()
  
  # Run PWAS for each complication outcome
  for (outcome in COMPLICATION_OUTCOMES) {
    
    if (!outcome %in% names(data)) {
      warning("Outcome not found in data: ", outcome)
      next
    }
    
    # Check if outcome has sufficient cases
    outcome_table <- table(data[[outcome]], useNA = "ifany")
    n_cases <- sum(outcome_table[names(outcome_table) == "1"], na.rm = TRUE)
    
    if (n_cases < MIN_CASES_FOR_ANALYSIS) {
      warning("Insufficient cases for outcome ", outcome, ": ", n_cases)
      next
    }
    
    cat("\n", rep("-", 50), "\n")
    cat("Analyzing outcome:", outcome, "(", n_cases, "cases )\n")
    cat(rep("-", 50), "\n")
    
    # Run logistic regression PWAS
    outcome_results <- run_pwas(
      data = data,
      outcome_var = outcome,
      protein_vars = protein_vars,
      adjustments = adjustments,
      analysis_type = "logistic",
      stratify_by = NULL,
      output_prefix = paste0("MAP_D_", outcome)
    )
    
    results_list[[outcome]] <- outcome_results
  }
  
  # Combine all results
  combined_results <- bind_rows(results_list, .id = "outcome_id")
  
  # Save combined results
  combined_output <- file.path(RESULTS_DIR, "MAP_D_complications_combined_results.csv")
  write_csv(combined_results, combined_output)
  
  cat("\nComplication PWAS completed.\n")
  cat("Combined results saved to:", combined_output, "\n")
  
  return(results_list)
}

#' Generate PWAS summary statistics
#' 
#' @param pwas_results PWAS results data frame
#' @param significance_threshold FDR threshold for significance
#' @return Summary statistics data frame
generate_pwas_summary <- function(pwas_results, significance_threshold = FDR_THRESHOLD) {
  
  if (!"FDR" %in% names(pwas_results)) {
    warning("FDR column not found in results")
    return(NULL)
  }
  
  summary_stats <- pwas_results %>%
    group_by(Phenotype, Subgroup) %>%
    summarise(
      total_proteins = n(),
      significant_proteins = sum(FDR < significance_threshold, na.rm = TRUE),
      median_pvalue = median(p.value, na.rm = TRUE),
      min_pvalue = min(p.value, na.rm = TRUE),
      proportion_significant = significant_proteins / total_proteins,
      .groups = "drop"
    ) %>%
    arrange(desc(proportion_significant))
  
  return(summary_stats)
}

#' Create volcano plot for PWAS results
#' 
#' @param pwas_results PWAS results data frame
#' @param title Plot title
#' @param output_file Output file path (optional)
#' @return ggplot object
create_volcano_plot <- function(pwas_results, title = "PWAS Volcano Plot", 
                               output_file = NULL) {
  
  load_required_packages(c("ggplot2", "ggrepel"))
  
  # Prepare data for plotting
  plot_data <- pwas_results %>%
    mutate(
      neg_log10_p = -log10(p.value),
      significant = FDR < FDR_THRESHOLD,
      abs_estimate = abs(estimate)
    ) %>%
    filter(!is.na(neg_log10_p), !is.na(estimate))
  
  # Create plot
  p <- ggplot(plot_data, aes(x = estimate, y = neg_log10_p)) +
    geom_point(aes(color = significant), alpha = 0.6, size = 1) +
    scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "red")) +
    geom_hline(yintercept = -log10(FDR_THRESHOLD), linetype = "dashed", 
               color = "red", alpha = 0.7) +
    labs(
      title = title,
      x = "Effect Size (Beta)",
      y = "-log10(P-value)",
      color = paste("FDR <", FDR_THRESHOLD)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      legend.position = "bottom"
    )
  
  # Add labels for top significant hits
  top_hits <- plot_data %>%
    filter(significant) %>%
    arrange(p.value) %>%
    head(10)
  
  if (nrow(top_hits) > 0) {
    p <- p + geom_text_repel(
      data = top_hits,
      aes(label = Protein),
      size = 3,
      max.overlaps = 10
    )
  }
  
  # Save plot if output file specified
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = PLOT_WIDTH, height = PLOT_HEIGHT, 
           dpi = PLOT_DPI, units = "in")
    cat("Volcano plot saved to:", output_file, "\n")
  }
  
  return(p)
}

cat("PWAS analysis module loaded successfully.\n")
