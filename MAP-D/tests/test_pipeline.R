# =============================================================================
# MAP-D Pipeline Test Suite
# Metabolic Atlas of Progression to Diabetes (MAP-D)
# =============================================================================

# Load required packages for testing
library(testthat)

# Source the main modules
source("../config.R")
source("../src/utils.R")

#' Test utility functions
test_that("Utility functions work correctly", {
  
  # Test create_glycemic_status function
  test_data <- data.frame(
    diabetes.y = c(0, 0, 1, 0, 0),
    HbA1c = c(35, 42, 55, 30, 45),
    stringsAsFactors = FALSE
  )
  
  result <- create_glycemic_status(test_data)
  
  expect_true("GlycemicStatus" %in% names(result))
  expect_equal(levels(result$GlycemicStatus), c("Normoglycemic", "Prediabetes", "Diabetes"))
  expect_equal(as.character(result$GlycemicStatus), 
               c("Normoglycemic", "Prediabetes", "Diabetes", "Normoglycemic", "Prediabetes"))
})

test_that("ALST calculation works correctly", {
  
  test_data <- data.frame(
    f.23129.0.0 = c(20, 25, 30),  # appendicular FFM
    x.sex = c(0, 1, 0),           # 0 = female, 1 = male
    stringsAsFactors = FALSE
  )
  
  result <- calculate_alst(test_data)
  
  expect_true("ALST" %in% names(result))
  
  # Check calculation: ALST = (0.958 × FFM) − (0.166 × sex) − 0.308
  expected_alst_1 <- (0.958 * 20) - (0.166 * 0) - 0.308  # Female
  expected_alst_2 <- (0.958 * 25) - (0.166 * 1) - 0.308  # Male
  expected_alst_3 <- (0.958 * 30) - (0.166 * 0) - 0.308  # Female
  
  expect_equal(result$ALST[1], expected_alst_1, tolerance = 0.001)
  expect_equal(result$ALST[2], expected_alst_2, tolerance = 0.001)
  expect_equal(result$ALST[3], expected_alst_3, tolerance = 0.001)
})

test_that("Data validation works correctly", {
  
  # Test with good data
  good_data <- data.frame(
    var1 = rnorm(200),
    var2 = rnorm(200),
    var3 = sample(c(0, 1), 200, replace = TRUE)
  )
  
  validation <- validate_data_quality(good_data, min_sample_size = 100, max_missing_prop = 0.1)
  
  expect_true(validation$meets_min_size)
  expect_true(validation$meets_missing_threshold)
  expect_true(validation$passes_validation)
  
  # Test with insufficient sample size
  small_data <- data.frame(
    var1 = rnorm(50),
    var2 = rnorm(50)
  )
  
  validation_small <- validate_data_quality(small_data, min_sample_size = 100)
  expect_false(validation_small$meets_min_size)
  expect_false(validation_small$passes_validation)
  
  # Test with high missing data
  missing_data <- data.frame(
    var1 = c(rnorm(50), rep(NA, 150)),  # 75% missing
    var2 = rnorm(200)
  )
  
  validation_missing <- validate_data_quality(missing_data, max_missing_prop = 0.1)
  expect_false(validation_missing$meets_missing_threshold)
  expect_false(validation_missing$passes_validation)
})

#' Test configuration settings
test_that("Configuration is properly loaded", {
  
  expect_true(exists("PROJECT_NAME"))
  expect_true(exists("METABOLIC_PHENOTYPES"))
  expect_true(exists("FDR_THRESHOLD"))
  expect_true(exists("RANDOM_SEED"))
  
  expect_equal(PROJECT_NAME, "MAP-D")
  expect_true(is.numeric(FDR_THRESHOLD))
  expect_true(FDR_THRESHOLD > 0 && FDR_THRESHOLD < 1)
  expect_true(is.numeric(RANDOM_SEED))
})

#' Test file checking utility
test_that("File existence checking works", {
  
  # Test with existing file
  temp_file <- tempfile()
  writeLines("test", temp_file)
  
  expect_true(check_file_exists(temp_file))
  
  # Test with non-existing file
  expect_false(check_file_exists("non_existent_file.txt"))
  
  # Clean up
  unlink(temp_file)
})

#' Test protein information extraction
test_that("Protein information extraction works", {
  
  # Test normal protein string
  protein_string <- "ENSG00000123456;Protein Name"
  result <- extract_protein_info(protein_string)
  
  expect_equal(result$code, "ENSG00000123456")
  expect_equal(result$name, "Protein Name")
  
  # Test protein string without name
  protein_string_no_name <- "ENSG00000123456"
  result_no_name <- extract_protein_info(protein_string_no_name)
  
  expect_equal(result_no_name$code, "ENSG00000123456")
  expect_true(is.na(result_no_name$name))
  
  # Test empty string
  result_empty <- extract_protein_info("")
  expect_true(is.na(result_empty$code))
  expect_true(is.na(result_empty$name))
  
  # Test NA input
  result_na <- extract_protein_info(NA)
  expect_true(is.na(result_na$code))
  expect_true(is.na(result_na$name))
})

#' Integration test with sample data
test_that("Pipeline components integrate correctly", {
  
  # Create sample data
  set.seed(123)
  n_samples <- 500
  n_proteins <- 10
  
  sample_data <- data.frame(
    eid = 1:n_samples,
    x.age = rnorm(n_samples, 50, 10),
    x.sex = sample(c(0, 1), n_samples, replace = TRUE),
    x.738 = sample(1:5, n_samples, replace = TRUE),
    diabetes.y = sample(c(0, 1), n_samples, replace = TRUE, prob = c(0.8, 0.2)),
    HbA1c = rnorm(n_samples, 40, 8),
    BMI = rnorm(n_samples, 27, 5),
    HDL = rnorm(n_samples, 1.3, 0.3),
    stringsAsFactors = FALSE
  )
  
  # Add protein columns
  for (i in 1:n_proteins) {
    sample_data[[paste0("protein_x.", i)]] <- rnorm(n_samples, 0, 1)
  }
  
  # Test glycemic status creation
  sample_data_with_status <- create_glycemic_status(sample_data)
  expect_true("GlycemicStatus" %in% names(sample_data_with_status))
  
  # Test data validation
  validation <- validate_data_quality(sample_data_with_status)
  expect_true(validation$passes_validation)
  
  # Test protein variable identification
  protein_vars <- grep("^protein_x\\.", names(sample_data_with_status), value = TRUE)
  expect_equal(length(protein_vars), n_proteins)
  expect_true(all(grepl("^protein_x\\.", protein_vars)))
})

# Run tests if script is executed directly
if (!interactive()) {
  cat("Running MAP-D pipeline tests...\n")
  
  # Set working directory to tests folder
  if (basename(getwd()) != "tests") {
    if (dir.exists("tests")) {
      setwd("tests")
    }
  }
  
  # Run all tests
  test_results <- test_dir(".", reporter = "summary")
  
  cat("\nTest results summary:\n")
  print(test_results)
  
  if (any(test_results$failed > 0)) {
    cat("\nSome tests failed. Please review and fix issues before running the pipeline.\n")
    quit(status = 1)
  } else {
    cat("\nAll tests passed! Pipeline is ready to run.\n")
  }
}
