#!/usr/bin/env Rscript

# Quick test script to validate refactored code
cat("=== MAP-D Refactoring Quick Test ===\n")

# Test 1: Basic module loading
cat("\n1. Testing module loading...\n")
tryCatch({
  source('config.R')
  cat("✓ Config loaded successfully\n")
  
  source('src/utils.R') 
  cat("✓ Utils loaded successfully\n")
  
  source('src/pwas_analysis.R')
  cat("✓ PWAS analysis module loaded successfully\n")
  
  # Check if PWAS helper functions are available
  if (exists("EWAS")) {
    cat("✓ Linear PWAS function (EWAS) available\n")
  } else {
    cat("✗ Linear PWAS function (EWAS) not found\n")
  }
  
  if (exists("EWASLogistic")) {
    cat("✓ Logistic PWAS function (EWASLogistic) available\n")
  } else {
    cat("✗ Logistic PWAS function (EWASLogistic) not found\n")
  }
  
}, error = function(e) {
  cat("✗ Module loading failed:", e$message, "\n")
  quit(status = 1)
})

# Test 2: Configuration validation
cat("\n2. Testing configuration...\n")
tryCatch({
  cat("Project:", PROJECT_NAME, "\n")
  cat("Random seed:", RANDOM_SEED, "\n")
  cat("FDR threshold:", FDR_THRESHOLD, "\n")
  cat("✓ Configuration variables accessible\n")
}, error = function(e) {
  cat("✗ Configuration test failed:", e$message, "\n")
})

# Test 3: Utility functions
cat("\n3. Testing utility functions...\n")
tryCatch({
  # Test glycemic status creation
  test_data <- data.frame(
    diabetes.y = c(0, 0, 1),
    HbA1c = c(35, 42, 55)
  )
  
  result <- create_glycemic_status(test_data)
  if ("GlycemicStatus" %in% names(result)) {
    cat("✓ create_glycemic_status function works\n")
  } else {
    cat("✗ create_glycemic_status function failed\n")
  }
  
  # Test data validation
  validation <- validate_data_quality(test_data, min_sample_size = 2)
  if (!is.null(validation$passes_validation)) {
    cat("✓ validate_data_quality function works\n")
  } else {
    cat("✗ validate_data_quality function failed\n")
  }
  
}, error = function(e) {
  cat("✗ Utility function test failed:", e$message, "\n")
})

# Test 4: File structure
cat("\n4. Testing file structure...\n")
required_dirs <- c("src", "data", "results", "docs", "tests", "legacy_scripts")
for (dir in required_dirs) {
  if (dir.exists(dir)) {
    cat("✓", dir, "directory exists\n")
  } else {
    cat("✗", dir, "directory missing\n")
  }
}

required_files <- c("main_analysis.R", "config.R", "requirements.R", "README.md")
for (file in required_files) {
  if (file.exists(file)) {
    cat("✓", file, "exists\n")
  } else {
    cat("✗", file, "missing\n")
  }
}

cat("\n=== Quick Test Complete ===\n")
cat("If all tests show ✓, the refactoring is working correctly!\n")
