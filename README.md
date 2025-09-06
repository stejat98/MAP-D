# MAP-D: Metabolic Atlas of Progression to Diabetes

[![R](https://img.shields.io/badge/R-%3E%3D4.0-blue)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

## Overview

The **Metabolic Atlas of Progression to Diabetes (MAP-D)** is a comprehensive proteomics analysis pipeline that maps patient progression from a healthy metabolic state to Type 2 Diabetes (T2D). This repository contains the complete computational framework used in our manuscript.

MAP-D leverages proteomic data on 2,923 proteins measured in a median of 47,963 UK Biobank participants to compute associations with hallmarks of metabolic disease across normoglycemia, prediabetes, and type 2 diabetes populations.

## Key Features

- **Comprehensive PWAS Analysis**: Protein-wide association studies across metabolic phenotypes and glycemic groups
- **Machine Learning Integration**: LASSO regression modeling for predictive biomarker discovery
- **Drug Target Discovery**: Integration with DrugBank for therapeutic target identification
- **External Validation**: Cross-validation with independent datasets
- **Reproducible Pipeline**: Fully automated analysis with comprehensive documentation


MAP-D provides:

1. **Proteomic signatures** that discriminate between patient subpopulations (e.g., obesity with normoglycemia vs. obesity with T2D)
2. **Predictive models** with R² up to 0.8 for metabolic traits in T2D
3. **Therapeutic insights** from GLP-1 agonist intervention trials
4. **Drug repurposing opportunities** through "therapeutically intransient" protein identification

## Installation

### Prerequisites

- R (≥ 4.0.0)
- Required R packages (automatically installed by the pipeline):
  - `tidyverse` (≥ 1.3.0)
  - `fst` (≥ 0.9.4)
  - `glmnet` (≥ 4.1)
  - `xml2` (≥ 1.3.0)
  - `tableone` (≥ 0.13.0)
  - `ggplot2` (≥ 3.3.0)
  - `ggrepel` (≥ 0.9.0)
  - `parallel`
  - `broom` (≥ 0.7.0)

### Setup

1. **Clone the repository:**
   ```bash
   git clone https://github.com/your-repo/MAP-D.git
   cd MAP-D
   ```

2. **Configure data paths:**
   Edit `config.R` to set your data directory paths:
   ```r
   # Modify these paths according to your system
   BASE_DATA_DIR <- "/path/to/your/data"
   RAW_DATA_DIR <- "/path/to/uk_biobank"
   ```

3. **Download required external data:**
   - UK Biobank proteomics data (Olink platform)
   - DrugBank XML database (https://go.drugbank.com/releases/latest)
   - External validation datasets (see Data Requirements)

## Usage

### Quick Start

#### Test with Synthetic Data (Recommended First Step)

Before running with your own data, test the pipeline with synthetic data to verify everything works:

```bash
# Test the pipeline with synthetic data
Rscript test_with_synthetic_data.R
```

This test script will:
- Generate synthetic proteomic and clinical data (500 samples, 10 proteins)
- Run linear PWAS analysis (metabolic phenotypes like BMI)
- Run logistic PWAS analysis (complications like CKD)
- Run stratified analysis by glycemic status
- Test utility functions and data validation
- Clean up test files automatically

**Expected output**: You should see ✓ marks for successful tests. If any test fails with ✗, check your R environment and dependencies.

#### Run Full Pipeline

Once the synthetic test passes, run the complete MAP-D analysis pipeline:

```bash
# Run full pipeline from raw data
Rscript main_analysis.R --full-pipeline

# Run with existing processed data
Rscript main_analysis.R

# Specify custom output directory
Rscript main_analysis.R --output-dir results/my_analysis/
```

### Interactive Analysis

For interactive analysis or custom modifications:

```r
# Load the pipeline
source("main_analysis.R")

# Run complete pipeline
results <- run_mapd_pipeline(
  run_full_pipeline = FALSE,
  skip_data_preprocessing = FALSE,
  results_directory = "results/"
)

# Access specific results
metabolic_pwas <- results$metabolic_pwas
lasso_models <- results$lasso_results
drug_targets <- results$drug_analysis
```

### Step-by-Step Analysis

Run individual components:

```r
# Data preprocessing
source("src/data_preprocessing.R")
processed_data <- main_data_preprocessing()

# PWAS analysis
source("src/pwas_analysis.R")
pwas_results <- run_metabolic_phenotype_pwas(processed_data, protein_vars, adjustments)

# LASSO modeling
source("src/lasso_analysis.R")
lasso_results <- run_metabolic_lasso_analysis(processed_data, protein_vars, baseline_vars)

# Drug target analysis
source("src/drug_analysis.R")
drug_results <- run_drug_target_analysis(pwas_results)
```

## Data Requirements

### Primary Data Sources

1. **UK Biobank Proteomics Data**
   - Olink proteomics platform data (Field 30901)
   - Main dataset (ukb34521.fst or similar)
   - Baseline assessment data

2. **Clinical Phenotypes**
   - BMI, HDL, LDL cholesterol
   - Triglyceride/HDL ratio
   - Systolic and diastolic blood pressure
   - HbA1c levels
   - Glycaemic status (Normoglycemic, Prediabetes, T2D)

3. **Complication Outcomes**
   - Coronary artery disease (CAD)
   - Chronic kidney disease (CKD)
   - Non-alcoholic fatty liver disease (NAFLD)

### External Validation Data

- STEP1 & STEP2 validation dataset (https://www.nature.com/articles/s41591-024-03355-2#data-availability)

### Drug Database

- DrugBank XML database (full version)
- Available from: https://go.drugbank.com/releases/latest
- Requires free academic license

## Directory Structure

```
MAP-D/
├── main_analysis.R                                                    # Main pipeline script
├── config.R                                                          # Configuration file
├── requirements.R                                                    # R package dependencies
├── test_with_synthetic_data.R                                       # Test script with synthetic data
├── README.md                                                         # This file
├── LICENSE                                                           # License file
├── src/                                                              # Source code modules
│   ├── utils.R                                                       # Utility functions
│   ├── data_preprocessing.R                                          # Data preprocessing
│   ├── pwas_analysis.R                                              # PWAS analysis
│   ├── lasso_analysis.R                                             # LASSO modeling
│   ├── drug_analysis.R                                              # Drug target analysis
│   ├── validation.R                                                 # Validation functions
│   ├── Baseline_PEWAS_Logistic_Functions_script.R                   # Logistic PWAS functions
│   └── cardiometabolic_Linear_reg_INRT_Baseline_PEWAS_Functions_Script.R  # Linear PWAS functions
├── data/                                                             # Data directory
│   ├── raw/                                                          # Raw data files
│   ├── processed/                                                    # Processed data
│   └── external/                                                     # External datasets
└── results/                                                          # Analysis results
    ├── figures/                                                      # Generated plots
    ├── tables/                                                       # Summary tables
    └── *.csv, *.RDS                                                  # Analysis output files
```

## Analysis Pipeline

### 1. Data Preprocessing (`src/data_preprocessing.R`)

- **Proteomics QC**: Missing data filtering, outlier detection
- **Sample Filtering**: Ancestry, completeness criteria
- **Variable Creation**: Glycemic status, derived metabolic measures
- **Complication Outcomes**: ICD-10 code processing for incident cases

### 2. PWAS Analysis (`src/pwas_analysis.R`)

- **Existing Helper Functions**: Uses validated PWAS functions from original analysis
  - Linear regression: `cardiometabolic_Linear_reg_INRT_Baseline_PEWAS_Functions_Script.R`
  - Logistic regression: `Baseline_PEWAS_Logistic_Functions_script.R`
- **Metabolic Phenotypes**: BMI, lipids, blood pressure, HbA1c
- **Stratified Analysis**: By glycemic status (normoglycemic, prediabetes, T2D)
- **Complication Analysis**: Cardiovascular, kidney, liver outcomes
- **Multiple Testing Correction**: FDR and Bonferroni methods

### 3. LASSO Modeling (`src/lasso_analysis.R`)

- **Model Comparison**: Baseline-only vs. protein-only vs. combined models
- **Cross-Validation**: Stratified train/test splits
- **Performance Metrics**: R² on held-out test sets
- **Feature Selection**: Identification of predictive proteins

### 4. Drug Target Analysis (`src/drug_analysis.R`)

- **DrugBank Integration**: Approved drug-target relationships
- **Directional Matching**: Therapeutic action consistency
- **Repurposing Opportunities**: Multi-target and multi-phenotype analysis
- **Network Visualization**: Drug-target interaction networks

## Output Files

### Main Results

- `MAP_D_complete_results.rds`: Complete pipeline results
- `MAP_D_analysis_summary.txt`: Human-readable summary
- `MAP_D_baseline_characteristics.csv`: Study population characteristics

### PWAS Results

- `MAP_D_[phenotype]_linear_results.csv`: Individual phenotype associations
- `MAP_D_metabolic_phenotypes_combined_results.csv`: Combined metabolic results
- `MAP_D_complications_combined_results.csv`: Complication outcome results

### LASSO Results

- `MAP_D_model_comparison_results.csv`: Model performance comparison
- `MAP_D_detailed_lasso_results.rds`: Detailed LASSO models
- `MAP_D_model_comparison.pdf`: Performance visualization

### Drug Analysis

- `drugbank_processed.csv`: Processed DrugBank data
- `drug_target_matches.csv`: Directionally consistent matches
- `drug_repurposing_by_[outcome/drug/protein].csv`: Repurposing summaries
- `MAP_D_drug_network.pdf`: Network visualization

## Reproducibility

### Session Information

The pipeline automatically captures:
- R version and platform information
- Package versions and dependencies
- Analysis timestamps and parameters
- Random seeds for reproducible results

### Configuration Management

All analysis parameters are centralized in `config.R`:
- Statistical thresholds (FDR, Bonferroni)
- Sample inclusion criteria
- ICD-10 code definitions
- File paths and directories

### Quality Control

Built-in validation includes:
- Sample size requirements
- Missing data thresholds
- Model convergence checks
- Result consistency validation

## Advanced Usage

### Custom Analysis

Modify `config.R` for custom analyses:

```r
# Custom significance thresholds
FDR_THRESHOLD <- 0.01
BONFERRONI_THRESHOLD <- 0.001

# Custom phenotype sets
CUSTOM_PHENOTYPES <- c("BMI", "waist_circumference", "body_fat_percentage")

# Custom ICD-10 codes
CUSTOM_DISEASE_CODES <- c("E11", "E111", "E112")  # T2D complications
```

### Parallel Processing

Enable parallel processing for large datasets:

```r
# In config.R
N_CORES <- 8  # Set to desired number of cores

# The pipeline automatically uses parallel processing where beneficial
```

### Memory Optimization

For large datasets, consider:

```r
# Process data in chunks
options(fst.use.index = TRUE)  # Enable FST indexing

# Reduce memory usage
gc()  # Garbage collection between major steps
```


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Data Availability
- **UK Biobank Data**: The data from the UK Biobank that support the findings of this study are available upon application (https://www.ukbiobank.ac.uk/register-apply/).
- **MAP-D Web Resource**: https://btierneyshiny.shinyapps.io/mapd-visualizer/
- **Figshare Repository**: https://doi.org/10.6084/m9.figshare.30007306.v1
- **GitHub Repository**: https://github.com/stejat98/MAP-D

## Acknowledgments

- UK Biobank participants and investigators
- Olink Proteomics for platform development
- DrugBank database maintainers
- R Core Team and package developers


---

