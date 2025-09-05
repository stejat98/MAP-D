# =============================================================================
# MAP-D Utility Functions
# Metabolic Atlas of Progression to Diabetes (MAP-D)
# =============================================================================

#' Load and validate required packages
#' 
#' @param packages Character vector of package names
#' @return NULL (packages are loaded)
load_required_packages <- function(packages) {
  cat("Loading required packages...\n")
  
  for (pkg in packages) {
    if (pkg == "tidyverse") {
      # Load core tidyverse packages individually since meta-package has issues
      core_packages <- c("dplyr", "tidyr", "readr", "purrr", "ggplot2", "tibble")
      for (core_pkg in core_packages) {
        if (requireNamespace(core_pkg, quietly = TRUE)) {
          library(core_pkg, character.only = TRUE)
          cat(paste("✓", core_pkg, "loaded\n"))
        } else {
          cat(paste("Note:", core_pkg, "not available, continuing...\n"))
        }
      }
    } else {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        cat("Installing missing package:", pkg, "\n")
        install.packages(pkg)
      }
      library(pkg, character.only = TRUE)
    }
  }
  
  cat("All packages loaded successfully.\n")
}

#' Create glycemic status variable
#' 
#' @param data Data frame containing diabetes status and HbA1c
#' @param diabetes_col Column name for diabetes status
#' @param hba1c_col Column name for HbA1c values
#' @param prediab_min Minimum HbA1c for prediabetes (default: 39)
#' @param prediab_max Maximum HbA1c for prediabetes (default: 48)
#' @return Data frame with GlycemicStatus column added
create_glycemic_status <- function(data, diabetes_col = "diabetes.y", 
                                  hba1c_col = "HbA1c", 
                                  prediab_min = 39, prediab_max = 48) {
  
  data$GlycemicStatus <- with(data, ifelse(
    get(diabetes_col) == 1, "Diabetes",
    ifelse(get(hba1c_col) >= prediab_min & get(hba1c_col) <= prediab_max, 
           "Prediabetes", "Normoglycemic")
  ))
  
  # Convert to factor with preferred order
  data$GlycemicStatus <- factor(data$GlycemicStatus, 
                               levels = c("Normoglycemic", "Prediabetes", "Diabetes"))
  
  return(data)
}

#' Apply strict normoglycemic filter
#' 
#' @param data Data frame with GlycemicStatus and HbA1c
#' @param hba1c_col Column name for HbA1c values
#' @param max_normo_hba1c Maximum HbA1c for strict normoglycemic (default: 39)
#' @return Filtered data frame
apply_strict_normoglycemic_filter <- function(data, hba1c_col = "HbA1c", 
                                            max_normo_hba1c = 39) {
  
  # Remove normoglycemic individuals with HbA1c > threshold
  excluded_indices <- which(data$GlycemicStatus == "Normoglycemic" & 
                           data[[hba1c_col]] > max_normo_hba1c)
  
  if (length(excluded_indices) > 0) {
    cat("Excluding", length(excluded_indices), 
        "normoglycemic individuals with HbA1c >", max_normo_hba1c, "\n")
    data <- data[-excluded_indices, ]
  }
  
  return(data)
}

#' Calculate ALST (Appendicular Lean Soft Tissue)
#' 
#' @param data Data frame containing appendicular FFM and sex
#' @param ffm_col Column name for appendicular FFM
#' @param sex_col Column name for sex (0 = female, 1 = male)
#' @return Data frame with ALST column added
calculate_alst <- function(data, ffm_col = "f.23129.0.0", sex_col = "x.sex") {
  
  # ALST (kg) = (0.958 × appendicular FFM) − (0.166 × S) − 0.308
  # where S = 0 for female, 1 for male
  data$ALST <- (0.958 * data[[ffm_col]]) - (0.166 * data[[sex_col]]) - 0.308
  
  return(data)
}

#' Validate data quality
#' 
#' @param data Data frame to validate
#' @param min_sample_size Minimum required sample size
#' @param max_missing_prop Maximum allowed missing proportion
#' @return List with validation results
validate_data_quality <- function(data, min_sample_size = 100, max_missing_prop = 0.1) {
  
  validation_results <- list()
  
  # Check sample size
  validation_results$sample_size <- nrow(data)
  validation_results$meets_min_size <- validation_results$sample_size >= min_sample_size
  
  # Check missing data
  missing_props <- sapply(data, function(x) sum(is.na(x)) / length(x))
  validation_results$missing_proportions <- missing_props
  validation_results$high_missing_vars <- names(missing_props)[missing_props > max_missing_prop]
  validation_results$meets_missing_threshold <- length(validation_results$high_missing_vars) == 0
  
  # Overall validation
  validation_results$passes_validation <- validation_results$meets_min_size && 
                                        validation_results$meets_missing_threshold
  
  return(validation_results)
}

#' Generate summary statistics by group
#' 
#' @param data Data frame
#' @param group_var Grouping variable name
#' @param vars Variables to summarize
#' @param categorical_vars Variables to treat as categorical
#' @return Summary table
generate_summary_stats <- function(data, group_var, vars, categorical_vars = NULL) {
  
  library(tableone)
  
  # Convert categorical variables to factors
  if (!is.null(categorical_vars)) {
    for (var in categorical_vars) {
      if (var %in% names(data)) {
        data[[var]] <- as.factor(as.character(data[[var]]))
      }
    }
  }
  
  # Create table one
  table1 <- CreateTableOne(vars = vars, strata = group_var, 
                          data = data, test = FALSE)
  
  return(table1)
}

#' Extract protein information from Olink format
#' 
#' @param protein_string Protein string in format "CODE;NAME"
#' @return Named list with code and name
extract_protein_info <- function(protein_string) {
  
  if (is.na(protein_string) || protein_string == "") {
    return(list(code = NA, name = NA))
  }
  
  parts <- strsplit(protein_string, ";", fixed = TRUE)[[1]]
  
  return(list(
    code = trimws(parts[1]),
    name = if (length(parts) > 1) trimws(parts[2]) else NA
  ))
}

#' Clean and standardize variable names
#' 
#' @param data Data frame
#' @param rename_map Named vector for specific renames
#' @return Data frame with cleaned names
clean_variable_names <- function(data, rename_map = NULL) {
  
  # Apply specific renames if provided
  if (!is.null(rename_map)) {
    old_names <- names(rename_map)
    new_names <- unname(rename_map)
    
    for (i in seq_along(old_names)) {
      if (old_names[i] %in% names(data)) {
        names(data)[names(data) == old_names[i]] <- new_names[i]
      }
    }
  }
  
  return(data)
}

#' Log analysis step
#' 
#' @param step_name Name of the analysis step
#' @param details Additional details to log
#' @return NULL (prints to console)
log_analysis_step <- function(step_name, details = NULL) {
  
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("[", timestamp, "] ", step_name, "\n", sep = "")
  
  if (!is.null(details)) {
    cat(details, "\n")
  }
  
  cat(rep("=", 60), "\n", sep = "")
}

#' Save session information for reproducibility
#' 
#' @param output_file Path to save session info
#' @return NULL (saves session info to file)
save_session_info <- function(output_file = "results/session_info.txt") {
  
  # Ensure directory exists
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  
  # Capture session information
  sink(output_file)
  cat("MAP-D Analysis Session Information\n")
  cat("Generated on:", as.character(Sys.time()), "\n\n")
  
  cat("R Version:\n")
  print(R.version.string)
  cat("\n")
  
  cat("Platform:\n")
  print(Sys.info())
  cat("\n")
  
  cat("Loaded Packages:\n")
  print(sessionInfo())
  cat("\n")
  
  sink()
  
  cat("Session information saved to:", output_file, "\n")
}

#' Check file existence and warn if missing
#' 
#' @param file_path Path to file
#' @param file_description Description of file for error message
#' @return Logical indicating if file exists
check_file_exists <- function(file_path, file_description = "File") {
  
  if (!file.exists(file_path)) {
    warning(file_description, " not found at:", file_path)
    return(FALSE)
  }
  
  return(TRUE)
}

cat("MAP-D utility functions loaded successfully.\n")
