# =============================================================================
# MAP-D Data Preprocessing Module
# Metabolic Atlas of Progression to Diabetes (MAP-D)
# =============================================================================

# Source required files
source("config.R")
source("src/utils.R")

#' Load and preprocess UK Biobank proteomics data
#' 
#' @param olink_file Path to Olink proteomics data file
#' @param output_file Path to save processed data
#' @return Wide-format proteomics data frame
load_and_process_proteomics <- function(olink_file = OLINK_DATA, 
                                       output_file = NULL) {
  
  log_analysis_step("Loading and Processing Proteomics Data")
  
  # Load required packages
  load_required_packages(c("tidyverse", "reshape2", "fst"))
  
  # Check if file exists
  if (!check_file_exists(olink_file, "Olink proteomics data")) {
    stop("Cannot proceed without proteomics data")
  }
  
  cat("Loading Olink proteomics data from:", olink_file, "\n")
  olink_data <- read_tsv(olink_file, show_col_types = FALSE)
  
  cat("Original data dimensions:", nrow(olink_data), "rows,", ncol(olink_data), "columns\n")
  
  # Reshape to wide format
  cat("Reshaping data to wide format...\n")
  olink_data_wide <- spread(olink_data, key = protein_id, value = result)
  
  # Add protein prefix to column names (except eid)
  protein_cols <- setdiff(names(olink_data_wide), "eid")
  names(olink_data_wide)[names(olink_data_wide) %in% protein_cols] <- 
    sprintf("protein_x.%s", protein_cols)
  
  cat("Wide format dimensions:", nrow(olink_data_wide), "rows,", ncol(olink_data_wide), "columns\n")
  
  # Save if output file specified
  if (!is.null(output_file)) {
    cat("Saving processed proteomics data to:", output_file, "\n")
    saveRDS(olink_data_wide, output_file)
  }
  
  return(olink_data_wide)
}

#' Load and merge UK Biobank main dataset with proteomics
#' 
#' @param ukb_file Path to UK Biobank main data file
#' @param proteomics_data Proteomics data frame (wide format)
#' @param output_file Path to save merged data
#' @return Merged data frame
merge_ukb_proteomics <- function(ukb_file = UKB_MAIN_DATA, 
                                proteomics_data = NULL,
                                output_file = NULL) {
  
  log_analysis_step("Merging UK Biobank Data with Proteomics")
  
  # Load required packages
  load_required_packages(c("fst", "dplyr"))
  
  # Load UKB main data if proteomics not provided
  if (is.null(proteomics_data)) {
    proteomics_data <- readRDS(file.path(BASE_DATA_DIR, "olink_data_reshaped.RDS"))
  }
  
  # Check if UKB file exists
  if (!check_file_exists(ukb_file, "UK Biobank main data")) {
    stop("Cannot proceed without UK Biobank main data")
  }
  
  cat("Loading UK Biobank main data...\n")
  # Note: This assumes the base data is already available
  # In practice, you would load the full UKB dataset here
  base_data_file <- file.path(BASE_DATA_DIR, "all_cause_mort_data_exposures_full.fst")
  
  if (check_file_exists(base_data_file, "Base UK Biobank data")) {
    data <- read.fst(base_data_file)
    cat("Base data dimensions:", nrow(data), "rows,", ncol(data), "columns\n")
  } else {
    stop("Base UK Biobank data file not found. Please ensure data preprocessing is complete.")
  }
  
  # Merge with proteomics
  cat("Merging with proteomics data...\n")
  merged_data <- left_join(data, proteomics_data, by = c("eid" = "protein_x.eid"))
  
  cat("Merged data dimensions:", nrow(merged_data), "rows,", ncol(merged_data), "columns\n")
  
  # Save merged data
  if (!is.null(output_file)) {
    cat("Saving merged data to:", output_file, "\n")
    write.fst(merged_data, output_file)
  }
  
  return(merged_data)
}

#' Apply sample filters and create derived variables
#' 
#' @param data Input data frame
#' @param ancestry_filter Vector of ancestry groups to include
#' @return Filtered and processed data frame
apply_sample_filters <- function(data, ancestry_filter = ANCESTRY_FILTER) {
  
  log_analysis_step("Applying Sample Filters and Creating Derived Variables")
  
  cat("Original sample size:", nrow(data), "\n")
  
  # Apply ancestry filter
  if ("f.21000.0.0" %in% names(data)) {
    data <- data %>% filter(f.21000.0.0 %in% ancestry_filter)
    cat("After ancestry filter:", nrow(data), "samples\n")
  } else {
    warning("Ancestry column not found, skipping ancestry filter")
  }
  
  # Standardize variable names
  rename_map <- c(
    "21001_0.0.y" = "BMI",
    "30760_0.0.y" = "HDL", 
    "30780_0.0.y" = "LDL",
    "4080_0.0.y" = "systolic_BP",
    "4079_0.0.y" = "diastolic_BP",
    "30620_0.0.y" = "ALT",
    "30750_0.0.y" = "HbA1c",
    "30870_0.0.y" = "triglycerides"
  )
  
  data <- clean_variable_names(data, rename_map)
  
  # Create derived variables
  cat("Creating derived variables...\n")
  
  # Triglyceride to HDL ratio
  if (all(c("triglycerides", "HDL") %in% names(data))) {
    data <- data %>% mutate(TRIG_HDL_RATIO = triglycerides / HDL)
  }
  
  # Add ALST if appendicular FFM data available
  if (all(c("f.23129.0.0", "x.sex") %in% names(data))) {
    data <- calculate_alst(data)
  }
  
  # Create glycemic status
  if (all(c("diabetes.y", "HbA1c") %in% names(data))) {
    data <- create_glycemic_status(data)
    
    # Apply strict normoglycemic filter
    data <- apply_strict_normoglycemic_filter(data)
    
    cat("Glycemic status distribution:\n")
    print(table(data$GlycemicStatus))
  }
  
  cat("Final processed sample size:", nrow(data), "\n")
  
  return(data)
}

#' Process complication outcomes from ICD-10 codes
#' 
#' @param ukb_data UK Biobank main dataset
#' @param output_file Path to save outcomes data
#' @return Data frame with binary complication outcomes
process_complication_outcomes <- function(ukb_data, output_file = NULL) {
  
  log_analysis_step("Processing Complication Outcomes from ICD-10 Codes")
  
  # Load required packages
  load_required_packages(c("dplyr", "tidyr", "lubridate", "purrr"))
  
  # Identify ICD-10 diagnosis and date columns
  dx_cols <- grep("^f.41270.0\\.", colnames(ukb_data), value = TRUE)
  dt_cols <- gsub("41270", "41280", dx_cols)
  
  cat("Found", length(dx_cols), "diagnosis columns\n")
  
  # Reshape to long format
  dx_long <- ukb_data %>%
    select(f.eid, all_of(dx_cols)) %>%
    pivot_longer(cols = -f.eid, names_to = "col", values_to = "icd10") %>%
    mutate(index = gsub(".*\\.", "", col)) %>%
    select(-col)
  
  dt_long <- ukb_data %>%
    select(f.eid, all_of(dt_cols)) %>%
    pivot_longer(cols = -f.eid, names_to = "col", values_to = "dx_date") %>%
    mutate(index = gsub(".*\\.", "", col)) %>%
    select(-col)
  
  # Join diagnosis and date
  dx_combined <- left_join(dx_long, dt_long, by = c("f.eid", "index")) %>%
    mutate(dx_date = ymd(dx_date))
  
  # Get baseline assessment date
  baseline_df <- ukb_data %>%
    select(f.eid, baseline_date = f.53.0.0) %>%
    mutate(baseline_date = ymd(baseline_date))
  
  # Function to get incident cases
  get_incident_cases <- function(codes, condition_name) {
    dx_combined %>%
      filter(icd10 %in% codes) %>%
      left_join(baseline_df, by = "f.eid") %>%
      filter(!is.na(dx_date) & dx_date > baseline_date) %>%
      group_by(f.eid) %>%
      summarise(
        event_date = min(dx_date),
        condition = condition_name,
        .groups = "drop"
      )
  }
  
  # Function to get prevalent cases (for exclusion)
  get_prevalent_cases <- function(codes) {
    dx_combined %>%
      filter(icd10 %in% codes) %>%
      left_join(baseline_df, by = "f.eid") %>%
      filter(!is.na(dx_date) & !is.na(baseline_date)) %>%
      filter(dx_date <= baseline_date) %>%
      distinct(f.eid)
  }
  
  # Process each complication
  complications <- list(
    ckd = CKD_CODES,
    nafld = NAFLD_CODES,
    hf = HEART_FAILURE_CODES
  )
  
  # Get all unique IDs
  all_ids <- dx_combined %>% select(f.eid) %>% distinct()
  
  outcome_dfs <- list()
  
  for (comp_name in names(complications)) {
    cat("Processing", comp_name, "outcomes...\n")
    
    codes <- complications[[comp_name]]
    
    # Get incident cases
    incident_cases <- get_incident_cases(codes, comp_name)
    
    # Get prevalent cases for exclusion
    prevalent_cases <- get_prevalent_cases(codes)
    
    # Create binary outcome
    binary_outcome <- all_ids %>%
      left_join(incident_cases %>% select(f.eid) %>% mutate(outcome = 1), 
                by = "f.eid") %>%
      mutate(outcome = ifelse(is.na(outcome), 0, outcome)) %>%
      anti_join(prevalent_cases, by = "f.eid") %>%  # Exclude prevalent cases
      rename(!!paste0("y_", comp_name) := outcome)
    
    outcome_dfs[[comp_name]] <- binary_outcome
    
    cat("  Incident cases:", sum(binary_outcome[[paste0("y_", comp_name)]]), "\n")
    cat("  Controls:", sum(binary_outcome[[paste0("y_", comp_name)]] == 0), "\n")
  }
  
  # Combine all outcomes
  outcomes_combined <- reduce(outcome_dfs, full_join, by = "f.eid")
  
  # Save outcomes
  if (!is.null(output_file)) {
    cat("Saving complication outcomes to:", output_file, "\n")
    write_csv(outcomes_combined, output_file)
  }
  
  return(outcomes_combined)
}

#' Perform quality control on proteomics data
#' 
#' @param data Data frame with proteomics data
#' @param min_detection_rate Minimum detection rate for proteins
#' @param max_missing_prop Maximum missing proportion
#' @return List with cleaned data and QC metrics
proteomics_quality_control <- function(data, 
                                     min_detection_rate = MIN_PROTEIN_DETECTION_RATE,
                                     max_missing_prop = MAX_MISSING_PROPORTION) {
  
  log_analysis_step("Proteomics Quality Control")
  
  # Identify protein columns
  protein_cols <- grep("^protein_x\\.", names(data), value = TRUE)
  cat("Found", length(protein_cols), "protein columns\n")
  
  if (length(protein_cols) == 0) {
    warning("No protein columns found")
    return(list(data = data, qc_metrics = NULL))
  }
  
  # Calculate missing proportions for each protein
  missing_props <- data[protein_cols] %>%
    summarise_all(~ sum(is.na(.)) / length(.)) %>%
    pivot_longer(everything(), names_to = "protein", values_to = "missing_prop")
  
  # Identify proteins to exclude
  proteins_to_exclude <- missing_props %>%
    filter(missing_prop > max_missing_prop) %>%
    pull(protein)
  
  cat("Excluding", length(proteins_to_exclude), "proteins with >", 
      max_missing_prop * 100, "% missing data\n")
  
  # Remove problematic proteins
  if (length(proteins_to_exclude) > 0) {
    data <- data %>% select(-all_of(proteins_to_exclude))
    protein_cols <- setdiff(protein_cols, proteins_to_exclude)
  }
  
  # Filter samples with complete proteomics data
  complete_protein_data <- data %>%
    filter(complete.cases(select(., all_of(protein_cols))))
  
  cat("Samples with complete proteomics data:", nrow(complete_protein_data), "\n")
  
  # QC metrics
  qc_metrics <- list(
    original_proteins = length(grep("^protein_x\\.", names(data), value = TRUE)) + length(proteins_to_exclude),
    final_proteins = length(protein_cols),
    excluded_proteins = proteins_to_exclude,
    original_samples = nrow(data),
    final_samples = nrow(complete_protein_data),
    missing_proportions = missing_props
  )
  
  return(list(
    data = complete_protein_data,
    qc_metrics = qc_metrics
  ))
}

#' Main data preprocessing pipeline
#' 
#' @param run_full_pipeline Logical, whether to run full pipeline or load existing data
#' @return Processed data frame ready for analysis
main_data_preprocessing <- function(run_full_pipeline = FALSE) {
  
  log_analysis_step("MAP-D Data Preprocessing Pipeline", 
                   "Starting comprehensive data preprocessing")
  
  if (run_full_pipeline) {
    # Step 1: Load and process proteomics data
    proteomics_data <- load_and_process_proteomics()
    
    # Step 2: Merge with UK Biobank data
    merged_data <- merge_ukb_proteomics(proteomics_data = proteomics_data)
    
    # Step 3: Process complications data
    ukb_main <- read.fst(UKB_MAIN_DATA)
    complications <- process_complication_outcomes(ukb_main)
    
    # Step 4: Merge complications
    final_data <- left_join(merged_data, complications, by = c("eid" = "f.eid"))
    
    # Step 5: Apply sample filters and create derived variables
    processed_data <- apply_sample_filters(final_data)
    
    # Step 6: Quality control
    qc_results <- proteomics_quality_control(processed_data)
    
    # Save final processed data
    output_file <- file.path(BASE_DATA_DIR, "MAP_D_processed_data.fst")
    write.fst(qc_results$data, output_file)
    
    cat("Data preprocessing complete. Final data saved to:", output_file, "\n")
    
    return(qc_results$data)
    
  } else {
    # Load existing processed data
    processed_file <- PROCESSED_DATA_FILE
    
    if (check_file_exists(processed_file, "Processed data")) {
      cat("Loading existing processed data from:", processed_file, "\n")
      data <- read.fst(processed_file)
      
      # Apply standard filters
      data <- apply_sample_filters(data)
      
      return(data)
    } else {
      stop("Processed data file not found. Set run_full_pipeline = TRUE to process from scratch.")
    }
  }
}

cat("Data preprocessing module loaded successfully.\n")
