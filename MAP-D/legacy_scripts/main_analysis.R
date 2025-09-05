#!/usr/bin/env Rscript

# =============================================================================
# MAP-D Main Analysis Pipeline
# Metabolic Atlas of Progression to Diabetes (MAP-D)
# 
# This script orchestrates the complete MAP-D analysis pipeline including:
# 1. Data preprocessing and quality control
# 2. Protein-wide association studies (PWAS)
# 3. LASSO regression modeling
# 4. Drug target analysis
# 5. Validation with external datasets
# 
# Usage:
#   Rscript main_analysis.R [--full-pipeline] [--skip-preprocessing] [--output-dir results/]
# =============================================================================

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
full_pipeline <- "--full-pipeline" %in% args
skip_preprocessing <- "--skip-preprocessing" %in% args
output_dir_arg <- grep("--output-dir", args, value = TRUE)
if (length(output_dir_arg) > 0) {
  RESULTS_DIR <- gsub("--output-dir=?", "", output_dir_arg)
}

# Source required modules
cat("Loading MAP-D analysis modules...\n")
source("config.R")
source("src/utils.R")
source("src/data_preprocessing.R")
source("src/pwas_analysis.R")
source("src/lasso_analysis.R")
source("src/drug_analysis.R")

# Load required packages
required_packages <- c(
  "tidyverse", "fst", "glmnet", "xml2", "tableone", 
  "ggplot2", "ggrepel", "parallel", "broom"
)

load_required_packages(required_packages)

#' Main analysis pipeline
#' 
#' @param run_full_pipeline Logical, whether to run complete pipeline from raw data
#' @param skip_data_preprocessing Logical, whether to skip data preprocessing
#' @param results_directory Directory for saving results
#' @return List with all analysis results
run_mapd_pipeline <- function(run_full_pipeline = FALSE, 
                             skip_data_preprocessing = FALSE,
                             results_directory = RESULTS_DIR) {
  
  log_analysis_step("MAP-D Analysis Pipeline", 
                   paste("Starting comprehensive analysis pipeline\n",
                        "Full pipeline:", run_full_pipeline, "\n",
                        "Skip preprocessing:", skip_data_preprocessing, "\n",
                        "Results directory:", results_directory))
  
  # Create results directory
  dir.create(results_directory, recursive = TRUE, showWarnings = FALSE)
  
  # Initialize results list
  pipeline_results <- list()
  
  # ==========================================================================
  # STEP 1: DATA PREPROCESSING
  # ==========================================================================
  
  if (!skip_data_preprocessing) {
    log_analysis_step("STEP 1: Data Preprocessing and Quality Control")
    
    # Run data preprocessing
    processed_data <- main_data_preprocessing(run_full_pipeline = run_full_pipeline)
    
    # Validate data quality
    validation_results <- validate_data_quality(processed_data)
    
    if (!validation_results$passes_validation) {
      warning("Data quality validation failed:")
      if (!validation_results$meets_min_size) {
        warning("  - Insufficient sample size: ", validation_results$sample_size)
      }
      if (!validation_results$meets_missing_threshold) {
        warning("  - High missing data in variables: ", 
                paste(validation_results$high_missing_vars, collapse = ", "))
      }
    }
    
    # Generate summary statistics
    if ("GlycemicStatus" %in% names(processed_data)) {
      summary_vars <- c(BASE_ADJUSTMENTS, METABOLIC_PHENOTYPES)
      available_vars <- intersect(summary_vars, names(processed_data))
      
      if (length(available_vars) > 0) {
        table1 <- generate_summary_stats(
          data = processed_data,
          group_var = "GlycemicStatus", 
          vars = available_vars,
          categorical_vars = c("x.sex", "x.738")
        )
        
        # Save table 1
        table1_df <- print(table1, showAllLevels = TRUE, quote = FALSE, 
                          noSpaces = TRUE, printToggle = FALSE)
        write_csv(as.data.frame(table1_df), 
                 file.path(results_directory, "MAP_D_baseline_characteristics.csv"))
      }
    }
    
    pipeline_results$processed_data <- processed_data
    pipeline_results$validation_results <- validation_results
    
  } else {
    log_analysis_step("STEP 1: Loading Existing Processed Data")
    processed_data <- main_data_preprocessing(run_full_pipeline = FALSE)
    pipeline_results$processed_data <- processed_data
  }
  
  # ==========================================================================
  # STEP 2: IDENTIFY PROTEIN VARIABLES AND ADJUSTMENTS
  # ==========================================================================
  
  log_analysis_step("STEP 2: Identifying Variables for Analysis")
  
  # Identify protein variables
  protein_vars <- grep("^protein_x\\.", names(processed_data), value = TRUE)
  cat("Identified", length(protein_vars), "protein variables\n")
  
  # Set up adjustment variables
  adjustments <- BASE_ADJUSTMENTS
  if (FASTING_TIME_VAR %in% names(processed_data)) {
    adjustments <- c(adjustments, FASTING_TIME_VAR)
  }
  
  cat("Base adjustments:", paste(adjustments, collapse = ", "), "\n")
  
  # Check for required variables
  missing_metabolic <- setdiff(METABOLIC_PHENOTYPES, names(processed_data))
  if (length(missing_metabolic) > 0) {
    warning("Missing metabolic phenotypes: ", paste(missing_metabolic, collapse = ", "))
  }
  
  missing_complications <- setdiff(COMPLICATION_OUTCOMES, names(processed_data))
  if (length(missing_complications) > 0) {
    warning("Missing complication outcomes: ", paste(missing_complications, collapse = ", "))
  }
  
  pipeline_results$protein_vars <- protein_vars
  pipeline_results$adjustments <- adjustments
  
  # ==========================================================================
  # STEP 3: METABOLIC PHENOTYPE PWAS
  # ==========================================================================
  
  log_analysis_step("STEP 3: Metabolic Phenotype PWAS Analysis")
  
  metabolic_pwas_results <- run_metabolic_phenotype_pwas(
    data = processed_data,
    protein_vars = protein_vars,
    adjustments = adjustments,
    add_fasting_adjustment = TRUE
  )
  
  pipeline_results$metabolic_pwas <- metabolic_pwas_results
  
  # Generate summary statistics
  for (phenotype in names(metabolic_pwas_results)) {
    if (!is.null(metabolic_pwas_results[[phenotype]])) {
      summary_stats <- generate_pwas_summary(metabolic_pwas_results[[phenotype]])
      if (!is.null(summary_stats)) {
        write_csv(summary_stats, 
                 file.path(results_directory, paste0("MAP_D_", phenotype, "_summary.csv")))
      }
    }
  }
  
  # ==========================================================================
  # STEP 4: COMPLICATION OUTCOME PWAS
  # ==========================================================================
  
  log_analysis_step("STEP 4: Complication Outcome PWAS Analysis")
  
  complication_pwas_results <- run_complication_pwas(
    data = processed_data,
    protein_vars = protein_vars,
    adjustments = adjustments,
    add_glycemic_adjustment = TRUE
  )
  
  pipeline_results$complication_pwas <- complication_pwas_results
  
  # ==========================================================================
  # STEP 5: LASSO REGRESSION ANALYSIS
  # ==========================================================================
  
  log_analysis_step("STEP 5: LASSO Regression Analysis")
  
  lasso_results <- run_metabolic_lasso_analysis(
    data = processed_data,
    protein_vars = protein_vars,
    baseline_vars = adjustments
  )
  
  pipeline_results$lasso_results <- lasso_results
  
  # Create model comparison plot
  if (!is.null(lasso_results$comparison)) {
    comparison_plot <- create_model_comparison_plot(
      lasso_results$comparison,
      output_file = file.path(results_directory, "figures", "MAP_D_model_comparison.pdf")
    )
    pipeline_results$comparison_plot <- comparison_plot
  }
  
  # ==========================================================================
  # STEP 6: DRUG TARGET ANALYSIS
  # ==========================================================================
  
  log_analysis_step("STEP 6: Drug Target Analysis")
  
  # Combine all PWAS results for drug analysis
  all_pwas_results <- bind_rows(
    bind_rows(metabolic_pwas_results, .id = "phenotype_source"),
    bind_rows(complication_pwas_results, .id = "phenotype_source")
  )
  
  if (nrow(all_pwas_results) > 0) {
    drug_analysis_results <- run_drug_target_analysis(
      pwas_results = all_pwas_results,
      output_dir = file.path(results_directory, "drug_analysis")
    )
    
    pipeline_results$drug_analysis <- drug_analysis_results
    
    # Create drug-target network if matches found
    if (!is.null(drug_analysis_results$drug_matches) && 
        nrow(drug_analysis_results$drug_matches) > 0) {
      
      network_plot <- create_drug_target_network(
        drug_analysis_results$drug_matches,
        output_file = file.path(results_directory, "figures", "MAP_D_drug_network.pdf")
      )
      pipeline_results$network_plot <- network_plot
    }
  } else {
    warning("No PWAS results available for drug analysis")
  }
  
  # ==========================================================================
  # STEP 7: EXTERNAL VALIDATION (if available)
  # ==========================================================================
  
  log_analysis_step("STEP 7: External Validation")
  
  # Check for validation data files
  step1_file <- STEP1_VALIDATION_FILE
  
  if (check_file_exists(step1_file, "STEP1 validation data")) {
    cat("External validation data found - running validation analysis\n")
    
    # Load validation functions (would need to be implemented)
    # validation_results <- run_external_validation(all_pwas_results, step1_file)
    # pipeline_results$validation <- validation_results
    
  } else {
    cat("External validation data not found - skipping validation step\n")
  }
  
  # ==========================================================================
  # STEP 8: GENERATE FINAL REPORTS AND VISUALIZATIONS
  # ==========================================================================
  
  log_analysis_step("STEP 8: Generating Final Reports")
  
  # Save session information
  save_session_info(file.path(results_directory, "session_info.txt"))
  
  # Create summary report
  create_analysis_summary(pipeline_results, results_directory)
  
  # Save complete pipeline results
  saveRDS(pipeline_results, file.path(results_directory, "MAP_D_complete_results.rds"))
  
  log_analysis_step("MAP-D Pipeline Completed Successfully", 
                   paste("Results saved to:", results_directory))
  
  return(pipeline_results)
}

#' Create analysis summary report
#' 
#' @param pipeline_results Results from run_mapd_pipeline
#' @param output_dir Output directory
#' @return NULL (creates summary files)
create_analysis_summary <- function(pipeline_results, output_dir) {
  
  summary_file <- file.path(output_dir, "MAP_D_analysis_summary.txt")
  
  sink(summary_file)
  
  cat("=" , rep("=", 60), "=\n", sep = "")
  cat("MAP-D Analysis Summary Report\n")
  cat("Generated on:", as.character(Sys.time()), "\n")
  cat("=" , rep("=", 60), "=\n", sep = "")
  
  # Data summary
  if (!is.null(pipeline_results$processed_data)) {
    data <- pipeline_results$processed_data
    cat("\nDATA SUMMARY:\n")
    cat("-", rep("-", 40), "-\n", sep = "")
    cat("Total samples:", nrow(data), "\n")
    
    if ("GlycemicStatus" %in% names(data)) {
      cat("Glycemic status distribution:\n")
      print(table(data$GlycemicStatus))
    }
    
    if (!is.null(pipeline_results$protein_vars)) {
      cat("Protein variables:", length(pipeline_results$protein_vars), "\n")
    }
  }
  
  # PWAS summary
  cat("\nPWAS RESULTS SUMMARY:\n")
  cat("-", rep("-", 40), "-\n", sep = "")
  
  if (!is.null(pipeline_results$metabolic_pwas)) {
    cat("Metabolic phenotype PWAS:\n")
    for (phenotype in names(pipeline_results$metabolic_pwas)) {
      results <- pipeline_results$metabolic_pwas[[phenotype]]
      if (!is.null(results) && nrow(results) > 0) {
        n_sig <- sum(results$FDR < FDR_THRESHOLD, na.rm = TRUE)
        cat("  ", phenotype, ": ", n_sig, " significant associations\n", sep = "")
      }
    }
  }
  
  if (!is.null(pipeline_results$complication_pwas)) {
    cat("Complication outcome PWAS:\n")
    for (outcome in names(pipeline_results$complication_pwas)) {
      results <- pipeline_results$complication_pwas[[outcome]]
      if (!is.null(results) && nrow(results) > 0) {
        n_sig <- sum(results$FDR < FDR_THRESHOLD, na.rm = TRUE)
        cat("  ", outcome, ": ", n_sig, " significant associations\n", sep = "")
      }
    }
  }
  
  # LASSO summary
  if (!is.null(pipeline_results$lasso_results)) {
    cat("\nLASSO MODELING SUMMARY:\n")
    cat("-", rep("-", 40), "-\n", sep = "")
    
    if (!is.null(pipeline_results$lasso_results$comparison)) {
      comp_results <- pipeline_results$lasso_results$comparison
      cat("Model performance (R² on test set):\n")
      
      for (outcome in unique(comp_results$Outcome)) {
        outcome_results <- comp_results %>% filter(Outcome == !!outcome)
        cat("  ", outcome, ":\n", sep = "")
        for (i in seq_len(nrow(outcome_results))) {
          cat("    ", outcome_results$ModelType[i], ": ", 
              sprintf("%.3f", outcome_results$R2[i]), "\n", sep = "")
        }
      }
    }
  }
  
  # Drug analysis summary
  if (!is.null(pipeline_results$drug_analysis)) {
    cat("\nDRUG TARGET ANALYSIS SUMMARY:\n")
    cat("-", rep("-", 40), "-\n", sep = "")
    
    drug_matches <- pipeline_results$drug_analysis$drug_matches
    if (!is.null(drug_matches) && nrow(drug_matches) > 0) {
      cat("Directionally consistent drug-target pairs:", nrow(drug_matches), "\n")
      cat("Unique drugs:", length(unique(drug_matches$drug_name)), "\n")
      cat("Unique protein targets:", length(unique(drug_matches$gene_symbol)), "\n")
    } else {
      cat("No directionally consistent drug-target pairs found\n")
    }
  }
  
  cat("\n", rep("=", 62), "\n", sep = "")
  cat("Analysis completed successfully.\n")
  cat("For detailed results, see individual output files in:", output_dir, "\n")
  cat(rep("=", 62), "\n", sep = "")
  
  sink()
  
  cat("Analysis summary saved to:", summary_file, "\n")
}

# =============================================================================
# MAIN EXECUTION
# =============================================================================

if (!interactive()) {
  cat("Starting MAP-D analysis pipeline...\n")
  
  # Run the complete pipeline
  results <- run_mapd_pipeline(
    run_full_pipeline = full_pipeline,
    skip_data_preprocessing = skip_preprocessing,
    results_directory = RESULTS_DIR
  )
  
  cat("MAP-D analysis pipeline completed successfully!\n")
  cat("Results saved to:", RESULTS_DIR, "\n")
  
} else {
  cat("MAP-D analysis pipeline loaded in interactive mode.\n")
  cat("Run: results <- run_mapd_pipeline() to execute the full pipeline.\n")
}
