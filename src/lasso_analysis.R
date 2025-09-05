# =============================================================================
# MAP-D LASSO Analysis Module
# Metabolic Atlas of Progression to Diabetes (MAP-D)
# =============================================================================

# Source required files
source("config.R")
source("src/utils.R")

#' Run LASSO regression analysis
#' 
#' @param data Input data frame
#' @param outcome_var Outcome variable name
#' @param predictor_vars Vector of predictor variable names
#' @param stratify_by Variable to stratify analysis by (optional)
#' @param alpha Elastic net mixing parameter (1 = LASSO, 0 = Ridge)
#' @param n_folds Number of cross-validation folds
#' @return List with model results and predictions
run_lasso_analysis <- function(data, outcome_var, predictor_vars, 
                              stratify_by = NULL, alpha = 1, n_folds = 10) {
  
  log_analysis_step(paste("Running LASSO Analysis for", outcome_var))
  
  # Load required packages
  load_required_packages(c("glmnet", "dplyr", "caret"))
  
  # Set seed for reproducibility
  set.seed(RANDOM_SEED)
  
  # Validate inputs
  if (!outcome_var %in% names(data)) {
    stop("Outcome variable not found in data: ", outcome_var)
  }
  
  missing_predictors <- setdiff(predictor_vars, names(data))
  if (length(missing_predictors) > 0) {
    warning("Some predictors not found in data: ", paste(missing_predictors, collapse = ", "))
    predictor_vars <- intersect(predictor_vars, names(data))
  }
  
  # Function to run LASSO on a single dataset
  run_single_lasso <- function(data_subset, outcome, predictors) {
    
    # Prepare matrices
    x_matrix <- as.matrix(data_subset[, predictors, drop = FALSE])
    y_vector <- data_subset[[outcome]]
    
    # Remove incomplete cases
    complete_cases <- complete.cases(x_matrix, y_vector)
    x_matrix <- x_matrix[complete_cases, , drop = FALSE]
    y_vector <- y_vector[complete_cases]
    
    # Check sample size
    if (length(y_vector) < MIN_SAMPLE_SIZE) {
      warning("Insufficient sample size: ", length(y_vector))
      return(NULL)
    }
    
    # Check for variance in outcome
    if (length(unique(y_vector)) < 2) {
      warning("Insufficient variance in outcome")
      return(NULL)
    }
    
    tryCatch({
      # Fit LASSO with cross-validation
      cv_fit <- cv.glmnet(x_matrix, y_vector, alpha = alpha, 
                         nfolds = n_folds, standardize = TRUE)
      
      # Extract coefficients at optimal lambda
      coefs <- coef(cv_fit, s = cv_fit$lambda.min)
      coef_df <- data.frame(
        variable = rownames(coefs),
        coefficient = as.numeric(coefs),
        stringsAsFactors = FALSE
      ) %>%
        filter(coefficient != 0, variable != "(Intercept)")
      
      # Calculate performance metrics
      predictions <- predict(cv_fit, s = cv_fit$lambda.min, newx = x_matrix)
      
      # Calculate R-squared
      ss_res <- sum((y_vector - predictions)^2)
      ss_tot <- sum((y_vector - mean(y_vector))^2)
      r_squared <- 1 - (ss_res / ss_tot)
      
      # Calculate cross-validated R-squared
      cv_r_squared <- 1 - (cv_fit$cvm[cv_fit$lambda == cv_fit$lambda.min] / 
                          var(y_vector))
      
      return(list(
        cv_fit = cv_fit,
        coefficients = coef_df,
        r_squared = r_squared,
        cv_r_squared = cv_r_squared,
        lambda_min = cv_fit$lambda.min,
        lambda_1se = cv_fit$lambda.1se,
        n_samples = length(y_vector),
        n_selected = nrow(coef_df)
      ))
      
    }, error = function(e) {
      warning("LASSO fitting failed: ", e$message)
      return(NULL)
    })
  }
  
  # Run analysis
  if (is.null(stratify_by)) {
    # Single analysis
    cat("Running LASSO for", length(predictor_vars), "predictors...\n")
    results <- run_single_lasso(data, outcome_var, predictor_vars)
    
    if (!is.null(results)) {
      results$stratum <- "All"
    }
    
    return(list(All = results))
    
  } else {
    # Stratified analysis
    cat("Running stratified LASSO by", stratify_by, "...\n")
    
    strata_levels <- unique(data[[stratify_by]])
    strata_levels <- strata_levels[!is.na(strata_levels)]
    
    results_by_strata <- list()
    
    for (stratum in strata_levels) {
      cat("  Processing stratum:", stratum, "\n")
      
      data_stratum <- data %>% filter(.data[[stratify_by]] == stratum)
      
      stratum_results <- run_single_lasso(data_stratum, outcome_var, predictor_vars)
      
      if (!is.null(stratum_results)) {
        stratum_results$stratum <- stratum
        results_by_strata[[as.character(stratum)]] <- stratum_results
      }
    }
    
    return(results_by_strata)
  }
}

#' Compare model performance across different predictor sets
#' 
#' @param data Input data frame
#' @param outcome_vars Vector of outcome variable names
#' @param predictor_sets Named list of predictor variable sets
#' @param train_proportion Proportion of data for training
#' @return Data frame with model comparison results
compare_model_performance <- function(data, outcome_vars, predictor_sets, 
                                    train_proportion = LASSO_TRAIN_PROPORTION) {
  
  log_analysis_step("Comparing Model Performance", 
                   paste("Testing", length(predictor_sets), "predictor sets on", 
                        length(outcome_vars), "outcomes"))
  
  # Load required packages
  load_required_packages(c("glmnet", "dplyr", "tidyr", "purrr"))
  
  # Set seed for reproducibility
  set.seed(RANDOM_SEED)
  
  # Create stratified train/test split by GlycemicStatus if available
  if ("GlycemicStatus" %in% names(data)) {
    split_data <- data %>%
      group_by(GlycemicStatus) %>%
      mutate(
        row_id = row_number(),
        n = n(),
        train_flag = row_id <= floor(train_proportion * n)
      ) %>%
      ungroup()
  } else {
    # Simple random split
    n_total <- nrow(data)
    train_indices <- sample(seq_len(n_total), size = floor(train_proportion * n_total))
    split_data <- data %>%
      mutate(
        train_flag = row_number() %in% train_indices
      )
  }
  
  train_data <- split_data %>% filter(train_flag)
  test_data <- split_data %>% filter(!train_flag)
  
  cat("Training samples:", nrow(train_data), "\n")
  cat("Testing samples:", nrow(test_data), "\n")
  
  # Function to evaluate single model
  evaluate_model <- function(outcome, predictors, train_df, test_df) {
    
    # Prepare training data
    x_train <- as.matrix(train_df[, predictors, drop = FALSE])
    y_train <- train_df[[outcome]]
    
    # Prepare test data
    x_test <- as.matrix(test_df[, predictors, drop = FALSE])
    y_test <- test_df[[outcome]]
    
    # Remove incomplete cases
    complete_train <- complete.cases(x_train, y_train)
    x_train <- x_train[complete_train, , drop = FALSE]
    y_train <- y_train[complete_train]
    
    complete_test <- complete.cases(x_test, y_test)
    x_test <- x_test[complete_test, , drop = FALSE]
    y_test <- y_test[complete_test]
    
    # Check sample sizes
    if (length(y_train) < MIN_SAMPLE_SIZE || length(y_test) < 10) {
      return(data.frame(R2 = NA, n_train = length(y_train), n_test = length(y_test)))
    }
    
    tryCatch({
      # Fit LASSO
      cv_fit <- cv.glmnet(x_train, y_train, alpha = 1, standardize = TRUE)
      
      # Predict on test set
      predictions <- predict(cv_fit, s = cv_fit$lambda.min, newx = x_test)
      
      # Calculate R-squared
      if (any(is.na(predictions)) || any(is.na(y_test)) || 
          length(predictions) != length(y_test)) {
        return(data.frame(R2 = NA, n_train = length(y_train), n_test = length(y_test)))
      }
      
      r_squared <- cor(y_test, predictions[, 1])^2
      
      return(data.frame(
        R2 = r_squared,
        n_train = length(y_train),
        n_test = length(y_test)
      ))
      
    }, error = function(e) {
      return(data.frame(R2 = NA, n_train = length(y_train), n_test = length(y_test)))
    })
  }
  
  # Run comparisons
  results <- expand_grid(
    Outcome = outcome_vars,
    ModelType = names(predictor_sets)
  ) %>%
    mutate(
      Predictors = map(ModelType, ~ predictor_sets[[.x]]),
      Results = map2(Outcome, Predictors, evaluate_model, 
                    train_df = train_data, test_df = test_data)
    ) %>%
    unnest(Results) %>%
    select(-Predictors)
  
  # Ensure proper factor ordering for plotting
  if ("BaselineOnly" %in% names(predictor_sets)) {
    model_order <- c("BaselineOnly", "ProteomicOnly", "BaselinePlusProteomic")
    available_models <- intersect(model_order, unique(results$ModelType))
    results$ModelType <- factor(results$ModelType, levels = available_models)
  }
  
  # Save results
  output_file <- file.path(RESULTS_DIR, "MAP_D_model_comparison_results.csv")
  write_csv(results, output_file)
  
  cat("Model comparison completed. Results saved to:", output_file, "\n")
  
  return(results)
}

#' Run comprehensive LASSO analysis for metabolic phenotypes
#' 
#' @param data Input data frame
#' @param protein_vars Vector of protein variable names
#' @param baseline_vars Vector of baseline covariate names
#' @return List with analysis results
run_metabolic_lasso_analysis <- function(data, protein_vars, baseline_vars) {
  
  log_analysis_step("Comprehensive LASSO Analysis for Metabolic Phenotypes")
  
  # Define predictor sets
  predictor_sets <- list(
    BaselineOnly = baseline_vars,
    ProteomicOnly = protein_vars,
    BaselinePlusProteomic = c(baseline_vars, protein_vars)
  )
  
  cat("Predictor sets defined:\n")
  for (set_name in names(predictor_sets)) {
    cat("  ", set_name, ":", length(predictor_sets[[set_name]]), "variables\n")
  }
  
  # Run model comparison
  comparison_results <- compare_model_performance(
    data = data,
    outcome_vars = METABOLIC_PHENOTYPES,
    predictor_sets = predictor_sets
  )
  
  # Run detailed LASSO analysis for each outcome with full model
  detailed_results <- list()
  
  for (outcome in METABOLIC_PHENOTYPES) {
    if (outcome %in% names(data)) {
      cat("\nRunning detailed LASSO for", outcome, "...\n")
      
      outcome_results <- run_lasso_analysis(
        data = data,
        outcome_var = outcome,
        predictor_vars = c(baseline_vars, protein_vars),
        stratify_by = "GlycemicStatus"
      )
      
      detailed_results[[outcome]] <- outcome_results
    }
  }
  
  # Save detailed results
  detailed_output <- file.path(RESULTS_DIR, "MAP_D_detailed_lasso_results.rds")
  saveRDS(detailed_results, detailed_output)
  
  return(list(
    comparison = comparison_results,
    detailed = detailed_results
  ))
}

#' Create model comparison forest plot
#' 
#' @param comparison_results Results from compare_model_performance
#' @param output_file Output file path (optional)
#' @return ggplot object
create_model_comparison_plot <- function(comparison_results, output_file = NULL) {
  
  load_required_packages(c("ggplot2", "dplyr"))
  
  # Order outcomes by R2 from the full model
  if ("BaselinePlusProteomic" %in% comparison_results$ModelType) {
    outcome_order <- comparison_results %>%
      filter(ModelType == "BaselinePlusProteomic") %>%
      arrange(R2) %>%
      pull(Outcome)
  } else {
    outcome_order <- unique(comparison_results$Outcome)
  }
  
  # Prepare data for plotting
  plot_data <- comparison_results %>%
    mutate(Outcome = factor(Outcome, levels = outcome_order)) %>%
    filter(!is.na(R2))
  
  # Create plot
  p <- ggplot(plot_data, aes(x = R2, y = Outcome, color = ModelType)) +
    geom_point(position = position_dodge(width = 0.6), size = 3) +
    scale_x_continuous(limits = c(0, 1), labels = scales::percent_format()) +
    scale_color_manual(values = c(
      "BaselineOnly" = "#1f77b4",
      "ProteomicOnly" = "#ff7f0e", 
      "BaselinePlusProteomic" = "#2ca02c"
    )) +
    labs(
      title = "Model Performance Comparison",
      subtitle = "R² on held-out test set",
      x = "R² (Test Set)",
      y = "Clinical Risk Factor",
      color = "Model Type"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      legend.position = "bottom",
      axis.title.y = element_blank()
    )
  
  # Save plot if output file specified
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = PLOT_WIDTH, height = PLOT_HEIGHT, 
           dpi = PLOT_DPI, units = "in")
    cat("Model comparison plot saved to:", output_file, "\n")
  }
  
  return(p)
}

#' Extract and summarize selected features from LASSO models
#' 
#' @param lasso_results Results from run_lasso_analysis
#' @param min_frequency Minimum frequency across strata for feature inclusion
#' @return Data frame with feature importance summary
summarize_lasso_features <- function(lasso_results, min_frequency = 0.5) {
  
  # Extract coefficients from all strata
  all_coefficients <- list()
  
  for (stratum_name in names(lasso_results)) {
    stratum_result <- lasso_results[[stratum_name]]
    
    if (!is.null(stratum_result) && !is.null(stratum_result$coefficients)) {
      coef_df <- stratum_result$coefficients %>%
        mutate(stratum = stratum_name)
      all_coefficients[[stratum_name]] <- coef_df
    }
  }
  
  if (length(all_coefficients) == 0) {
    return(NULL)
  }
  
  # Combine coefficients
  combined_coefs <- bind_rows(all_coefficients)
  
  # Summarize feature importance
  feature_summary <- combined_coefs %>%
    group_by(variable) %>%
    summarise(
      n_selected = n(),
      n_strata = length(unique(stratum)),
      frequency = n_selected / length(lasso_results),
      mean_coefficient = mean(coefficient),
      median_coefficient = median(coefficient),
      coefficient_range = max(coefficient) - min(coefficient),
      .groups = "drop"
    ) %>%
    filter(frequency >= min_frequency) %>%
    arrange(desc(frequency), desc(abs(mean_coefficient)))
  
  return(feature_summary)
}

cat("LASSO analysis module loaded successfully.\n")
