#!/usr/bin/env Rscript

# Test MAP-D refactored pipeline with synthetic data
cat("=== MAP-D Synthetic Data Test ===\n")

# Load modules
cat("Loading modules...\n")
source("config.R")
source("src/utils.R")
source("src/pwas_analysis.R")

# Create synthetic data that mimics your real data structure
cat("Creating synthetic data...\n")
set.seed(RANDOM_SEED)

n_samples <- 500
n_proteins <- 10

# Create synthetic data
synthetic_data <- data.frame(
  eid = 1:n_samples,
  
  # Demographics (matching your adjustments)
  x.age = rnorm(n_samples, 50, 15),
  x.sex = sample(c(0, 1), n_samples, replace = TRUE),
  x.738 = sample(1:5, n_samples, replace = TRUE),
  f.74.0.0 = rnorm(n_samples, 8, 3),  # fasting time
  
  # Diabetes status for glycemic classification
  diabetes.y = sample(c(0, 1), n_samples, replace = TRUE, prob = c(0.8, 0.2)),
  HbA1c = rnorm(n_samples, 40, 8),
  
  # Metabolic phenotypes
  BMI = rnorm(n_samples, 27, 5),
  HDL = rnorm(n_samples, 1.3, 0.3),
  LDL = rnorm(n_samples, 3.2, 0.8),
  systolic_BP = rnorm(n_samples, 130, 20),
  diastolic_BP = rnorm(n_samples, 80, 10),
  triglycerides = rnorm(n_samples, 1.7, 0.8),
  
  # Complication outcomes
  y_ckd = sample(c(0, 1), n_samples, replace = TRUE, prob = c(0.9, 0.1)),
  y_nafld = sample(c(0, 1), n_samples, replace = TRUE, prob = c(0.85, 0.15)),
  y_hf = sample(c(0, 1), n_samples, replace = TRUE, prob = c(0.95, 0.05)),
  cad = sample(c(0, 1), n_samples, replace = TRUE, prob = c(0.9, 0.1))
)

# Add protein columns
for (i in 1:n_proteins) {
  synthetic_data[[paste0("protein_x.", i)]] <- rnorm(n_samples, 0, 1)
}

# Create derived variables
synthetic_data$TRIG_HDL_RATIO <- synthetic_data$triglycerides / synthetic_data$HDL

# Create glycemic status
synthetic_data <- create_glycemic_status(synthetic_data)

cat("Synthetic data created:\n")
cat("- Samples:", nrow(synthetic_data), "\n")
cat("- Variables:", ncol(synthetic_data), "\n")
cat("- Proteins:", n_proteins, "\n")

# Check glycemic status distribution
cat("- Glycemic status distribution:\n")
print(table(synthetic_data$GlycemicStatus))

# Test 1: Linear PWAS (metabolic phenotypes)
cat("\n=== Testing Linear PWAS ===\n")
protein_vars <- grep("^protein_x\\.", names(synthetic_data), value = TRUE)
adjustments <- c("x.age", "x.sex", "x.738")

tryCatch({
  linear_results <- run_pwas(
    data = synthetic_data,
    outcome_var = "BMI",
    protein_vars = protein_vars[1:3],  # Test with first 3 proteins
    adjustments = adjustments,
    analysis_type = "linear",
    output_prefix = "test_synthetic_linear"
  )
  
  if (!is.null(linear_results) && nrow(linear_results) > 0) {
    cat("✓ Linear PWAS test successful!\n")
    cat("- Results shape:", nrow(linear_results), "x", ncol(linear_results), "\n")
    cat("- Columns:", paste(names(linear_results), collapse = ", "), "\n")
    
    if ("p.value" %in% names(linear_results)) {
      cat("- P-value range:", min(linear_results$p.value, na.rm = TRUE), 
          "to", max(linear_results$p.value, na.rm = TRUE), "\n")
    }
  } else {
    cat("✗ Linear PWAS test failed - no results\n")
  }
  
}, error = function(e) {
  cat("✗ Linear PWAS test failed:", e$message, "\n")
})

# Test 2: Logistic PWAS (complications)
cat("\n=== Testing Logistic PWAS ===\n")
tryCatch({
  logistic_results <- run_pwas(
    data = synthetic_data,
    outcome_var = "y_ckd",
    protein_vars = protein_vars[1:3],  # Test with first 3 proteins
    adjustments = adjustments,
    analysis_type = "logistic",
    output_prefix = "test_synthetic_logistic"
  )
  
  if (!is.null(logistic_results) && nrow(logistic_results) > 0) {
    cat("✓ Logistic PWAS test successful!\n")
    cat("- Results shape:", nrow(logistic_results), "x", ncol(logistic_results), "\n")
    cat("- Columns:", paste(names(logistic_results), collapse = ", "), "\n")
    
    if ("p.value" %in% names(logistic_results)) {
      cat("- P-value range:", min(logistic_results$p.value, na.rm = TRUE), 
          "to", max(logistic_results$p.value, na.rm = TRUE), "\n")
    }
  } else {
    cat("✗ Logistic PWAS test failed - no results\n")
  }
  
}, error = function(e) {
  cat("✗ Logistic PWAS test failed:", e$message, "\n")
})

# Test 3: Stratified Analysis
cat("\n=== Testing Stratified PWAS ===\n")
tryCatch({
  stratified_results <- run_pwas(
    data = synthetic_data,
    outcome_var = "BMI",
    protein_vars = protein_vars[1:2],  # Test with 2 proteins
    adjustments = adjustments,
    analysis_type = "linear",
    stratify_by = "GlycemicStatus",
    output_prefix = "test_synthetic_stratified"
  )
  
  if (!is.null(stratified_results) && nrow(stratified_results) > 0) {
    cat("✓ Stratified PWAS test successful!\n")
    cat("- Results shape:", nrow(stratified_results), "x", ncol(stratified_results), "\n")
    
    if ("Subgroup" %in% names(stratified_results)) {
      cat("- Subgroups:", paste(unique(stratified_results$Subgroup), collapse = ", "), "\n")
    }
  } else {
    cat("✗ Stratified PWAS test failed - no results\n")
  }
  
}, error = function(e) {
  cat("✗ Stratified PWAS test failed:", e$message, "\n")
})

# Test 4: Utility Functions
cat("\n=== Testing Utility Functions ===\n")
tryCatch({
  # Test data validation
  validation <- validate_data_quality(synthetic_data)
  cat("✓ Data validation works - passes:", validation$passes_validation, "\n")
  
  # Test summary stats
  if (exists("generate_summary_stats")) {
    table1 <- generate_summary_stats(
      data = synthetic_data,
      group_var = "GlycemicStatus",
      vars = c("x.age", "BMI", "HDL"),
      categorical_vars = "x.sex"
    )
    cat("✓ Summary statistics generation works\n")
  }
  
}, error = function(e) {
  cat("✗ Utility function test failed:", e$message, "\n")
})

# Cleanup test files
cat("\n=== Cleaning up test files ===\n")
test_files <- list.files(pattern = "^test_synthetic.*\\.(RDS|csv)$")
if (length(test_files) > 0) {
  file.remove(test_files)
  cat("Removed", length(test_files), "test files\n")
}

cat("\n=== Synthetic Data Test Complete ===\n")
cat("If you see ✓ for all major tests, it is working correctly!\n")
cat("\nNext steps:\n")
cat("1. Update config.R with your actual data paths\n")
cat("2. Test with your real data\n")
cat("3. Run the full pipeline: source('main_analysis.R'); run_mapd_pipeline()\n")
