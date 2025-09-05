# MAP-D Testing Guide

This guide provides a systematic approach to test the refactored MAP-D codebase and ensure nothing was broken during the refactoring process.

## 🧪 **Testing Strategy**

### **Phase 1: Quick Validation Tests**

#### Test 1: Basic Module Loading
```bash
# Run the quick test script
Rscript test_quick.R
```

This tests:
- ✅ All modules load without syntax errors
- ✅ PWAS helper functions are accessible
- ✅ Configuration variables are set correctly
- ✅ Utility functions work
- ✅ File structure is complete

#### Test 2: Package Dependencies
```r
# In R console
source("requirements.R")
check_package_versions()
```

### **Phase 2: Function-Level Testing**

#### Test 3: PWAS Helper Functions Direct Test
```r
# Test with minimal synthetic data
source("config.R")
source("src/pwas_analysis.R")

# Create minimal test data
set.seed(123)
test_data <- data.frame(
  eid = 1:100,
  x.age = rnorm(100, 50, 10),
  x.sex = sample(c(0, 1), 100, replace = TRUE),
  x.738 = sample(1:5, 100, replace = TRUE),
  protein_x.1 = rnorm(100, 0, 1),
  protein_x.2 = rnorm(100, 0, 1),
  BMI = rnorm(100, 27, 5),
  y_test = sample(c(0, 1), 100, replace = TRUE)
)

# Test linear PWAS function
if (exists("EWAS")) {
  cat("Testing EWAS function...\n")
  EWAS(data = test_data, 
       depvar = "BMI", 
       adjustments = "x.age + x.sex",
       exposures = c("protein_x.1", "protein_x.2"),
       outFileName = "test_linear")
  
  if (file.exists("test_linear.RDS")) {
    cat("✓ Linear PWAS function works!\n")
    results <- readRDS("test_linear.RDS")
    print(head(results))
  }
}

# Test logistic PWAS function  
if (exists("EWASLogistic")) {
  cat("Testing EWASLogistic function...\n")
  EWASLogistic(data = test_data,
              depvar = "y_test",
              adjustments = "x.age + x.sex", 
              exposures = c("protein_x.1", "protein_x.2"),
              outFileName = "test_logistic")
  
  if (file.exists("test_logistic.RDS")) {
    cat("✓ Logistic PWAS function works!\n")
    results <- readRDS("test_logistic.RDS")
    print(head(results))
  }
}
```

### **Phase 3: Integration Testing with Real Paths**

#### Test 4: Update Config for Your System
```r
# Edit config.R with your actual data paths
BASE_DATA_DIR <- "/your/actual/data/path"  # Update this
RAW_DATA_DIR <- "/your/actual/raw/data/path"  # Update this
RESULTS_DIR <- "/your/actual/results/path"  # Update this
```

#### Test 5: Test Data Loading (if data available)
```r
# Test data preprocessing functions
source("main_analysis.R")

# Test with your actual data paths
if (file.exists(PROCESSED_DATA_FILE)) {
  cat("Testing data loading...\n")
  data <- main_data_preprocessing(run_full_pipeline = FALSE)
  cat("✓ Data loaded successfully!\n")
  cat("Data dimensions:", nrow(data), "x", ncol(data), "\n")
  
  # Check protein variables
  protein_vars <- grep("^protein_x\\.", names(data), value = TRUE)
  cat("Found", length(protein_vars), "protein variables\n")
}
```

### **Phase 4: End-to-End Pipeline Test**

#### Test 6: Mini Pipeline Run
```r
# Run a small subset analysis to test the full pipeline
source("main_analysis.R")

# Create test subset if you have real data
if (exists("data") && nrow(data) > 1000) {
  cat("Running mini pipeline test...\n")
  
  # Take small random sample
  test_subset <- data %>% 
    sample_n(min(1000, nrow(data))) %>%
    select(eid, x.age, x.sex, x.738, BMI, 
           matches("^protein_x\\."), GlycemicStatus) %>%
    na.omit()
  
  # Test PWAS on subset
  protein_test <- grep("^protein_x\\.", names(test_subset), value = TRUE)[1:5]  # Just 5 proteins
  
  results <- run_pwas(
    data = test_subset,
    outcome_var = "BMI", 
    protein_vars = protein_test,
    adjustments = c("x.age", "x.sex"),
    analysis_type = "linear",
    output_prefix = "test_mini"
  )
  
  if (!is.null(results)) {
    cat("✓ Mini pipeline test successful!\n")
    cat("Results shape:", nrow(results), "x", ncol(results), "\n")
  }
}
```

## 🚀 **Quick Start Testing (Recommended)**

### **Option A: If you have your data readily accessible**
1. Update paths in `config.R`
2. Run: `Rscript test_quick.R`
3. If successful, run a mini analysis on a subset

### **Option B: If you need to set up data paths first**
1. Check what data files you have available
2. Update the paths in `config.R` to point to your actual data
3. Test with the synthetic data approach first

### **Option C: Compare with Legacy Results**
```r
# Run the same analysis with both old and new code
# Compare outputs to ensure consistency

# Old way (from legacy_scripts)
source("legacy_scripts/run_metabolic_hallmark_phenotypes_PWAS.R")

# New way (refactored)
source("main_analysis.R")
results_new <- run_mapd_pipeline()

# Compare key results
# (You'll need to adapt this based on your specific comparisons)
```

## 🔍 **What to Look For**

### ✅ **Success Indicators:**
- All modules load without errors
- PWAS helper functions are accessible (`EWAS`, `EWASLogistic`)
- Test data produces expected output format
- File paths resolve correctly
- Results have expected columns: `Protein`, `Phenotype`, `estimate`, `p.value`, `FDR`

### ❌ **Failure Indicators:**
- Module loading errors
- Missing function errors
- File path errors
- Unexpected result formats
- Missing output files

## 🛠 **Troubleshooting**

### Common Issues:
1. **R not found**: Ensure R is installed and in PATH
2. **Package missing**: Run `source("requirements.R"); install_mapd_packages()`
3. **File paths wrong**: Update paths in `config.R`
4. **Helper functions not found**: Check that helper scripts are in `src/` directory

### Quick Fixes:
```bash
# If R/Rscript not found, try:
/usr/local/bin/R --version  # or wherever R is installed

# If packages missing:
R -e "install.packages(c('tidyverse', 'fst', 'broom'))"

# If paths wrong, check what files you actually have:
ls -la /your/data/directory/
```

## 📋 **Testing Checklist**

- [ ] Quick test script runs successfully
- [ ] All required packages installed
- [ ] PWAS helper functions accessible
- [ ] Config paths updated for your system
- [ ] Synthetic data test works
- [ ] Real data loads (if available)
- [ ] Mini pipeline runs successfully
- [ ] Results match expected format

Once all tests pass, your refactored MAP-D pipeline is ready for production use! 🎉
