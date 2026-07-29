# MAP-D: Metabolic Atlas of the Proteome in Diabetes

[![R](https://img.shields.io/badge/R-%3E%3D4.0-blue)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

## Overview

The **Metabolic Atlas of the Proteome in Diabetes (MAP-D)** is a comprehensive proteomics analysis pipeline that maps patient progression from a healthy metabolic state to Type 2 Diabetes (T2D). This repository contains the complete computational framework used in our manuscript.

MAP-D leverages proteomic data on 2,923 proteins measured in a median of 47,963 UK Biobank participants to compute associations with hallmarks of metabolic disease across normoglycemia, prediabetes, and type 2 diabetes populations.

## Key Features

- **Comprehensive PWAS Analysis**: Protein-wide association studies across metabolic phenotypes and glycemic groups
- **Predicitive analysis**: LASSO regression modeling for predictive biomarker discovery
- **Drug Target Discovery**: Integration with DrugBank for therapeutic target identification
- **Split-sample GWAS and Mendelian randomization**: REGENIE hallmark GWAS in a held-out sample plus forward, reverse, and deCODE-replication MR under `gwas_regenie/`, `mr_twosample/`, and `figures_tables/` (see [below](#split-sample-gwas-and-mendelian-randomization))
- **Interactive web atlas**: A Shiny companion app (`app/`) for exploring every result layer — atlas, trials, MR, evidence funnel, CAD, and drug targeting (see [below](#interactive-web-application-app))


MAP-D provides:

1. **Proteomic signatures** that discriminate between patient subpopulations (e.g., obesity with normoglycemia vs. obesity with T2D)
2. **Predictive models** with R² up to 0.8 for metabolic traits in T2D
3. **Therapeutic insights** from GLP-1 agonist intervention trials
4. **Drug repurposing opportunities** through "therapeutically intransigent" protein identification

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
├── app/                                                              # Interactive Shiny web app (manuscript companion)
│   ├── app.R                                                         # Shiny UI + server (bslib navbar, 9 sections)
│   ├── plots.R                                                       # Figure builders (volcano, triangulation, MR, CAD, ...)
│   ├── build_data.R                                                  # Builds app/data/mapd.rds from the supplementary workbook
│   ├── data/mapd.rds                                                 # Pre-built app data bundle (single source of truth)
│   ├── SUPP_TABLE_1_05_20_2026.xlsx                                  # Supplementary workbook the app data is built from
│   ├── manifest.json                                                 # Posit Connect Cloud deployment manifest
│   ├── deploy.R                                                      # shinyapps.io deploy helper
│   └── README.md                                                     # App-specific documentation
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
├── gwas_regenie/                                                     # Split-sample REGENIE GWAS
│   ├── scripts/                                                      # Input generation, submission, post-processing
│   └── docs/                                                         # Run plans and summary-statistic paths
├── mr_twosample/                                                     # Two-sample and bidirectional MR
│   └── scripts/                                                      # Forward, reverse, and deCODE replication MR
├── figures_tables/                                                   # Manuscript figures and summary tables
│   └── scripts/
├── docs/genetics_pipeline/                                           # Genetics pipeline documentation
│   ├── WALKTHROUGH_SPLIT_SAMPLE_MR.md                                # End-to-end split-sample GWAS + MR guide
│   ├── FILE_MANIFEST.md                                              # Provenance of every staged script
│   └── build_github_staging.R                                        # Reproducible staging script
└── results/                                                          # Analysis results
    ├── figures/                                                      # Generated plots
    ├── tables/                                                       # Summary tables
    └── *.csv, *.RDS                                                  # Analysis output files
```

## Interactive Web Application (`app/`)

The `app/` directory contains the interactive **MAP-D atlas** — a Shiny application that
accompanies the manuscript and lets readers explore every result layer. It is organised as a
guided tour through the paper's triangulation narrative (Figures 1–5):

- **Overview** — study design and triangulation framework
- **Atlas** — proteome-wide hallmark associations (volcano, top proteins, cross-hallmark correlation)
- **Triangulation & Trials** — UK Biobank vs. STEP 1/2 semaglutide effects (reversion vs. persistent)
- **Bidirectional MR** — protein-first and hallmark-first Mendelian-randomization architecture
- **Evidence Funnel** — winnowing of proteins across evidence layers
- **Persistent → CAD** — incident coronary-artery-disease associations
- **Drug Targeting** — DrugBank-matched approved drugs
- **Protein Explorer** — every line of evidence for any single protein

### Data

The app reads a single pre-built bundle, `app/data/mapd.rds`, produced from the supplementary
workbook (`SUPP_TABLE_1_05_20_2026.xlsx`) by `app/build_data.R` — the single source of truth.
Rebuild it whenever the workbook changes:

```bash
cd app
Rscript build_data.R      # writes data/mapd.rds and asserts the manuscript's headline numbers
```

### Run locally

```bash
cd app
Rscript -e "shiny::runApp('.', launch.browser=TRUE)"
```

Requires `shiny`, `bslib`, `plotly`, `DT`, `dplyr`, `tidyr`, `stringr`, and `ggplot2`.

### Deploy

The app is deployed as a public web resource (see [Data Availability](#data-availability)).
`app/manifest.json` (generated with `rsconnect::writeManifest()`) makes the `app/` directory
directly publishable on **Posit Connect Cloud** from this repository; `app/deploy.R` is a helper
for shinyapps.io. Only `app.R`, `plots.R`, and `data/mapd.rds` are required at runtime.

## Split-Sample GWAS and Mendelian Randomization

The genetics arm of MAP-D tests causal relationships between plasma proteins and three
cardiometabolic hallmarks: **BMI**, **HbA1c**, and the **triglyceride/HDL ratio
(`TRIG_HDL_RATIO`)**. The GWAS and MR analyses reported in the manuscript were run in the
**full cohort only**; the glycemic subgroups used elsewhere in MAP-D are not part of this arm.

To avoid sample overlap between the exposure and outcome GWAS — a key assumption of
two-sample MR — UK Biobank is split into two non-overlapping subsets:

- **Proteomics sample**: participants with Olink measurements, used for the protein GWAS
- **Held-out sample**: the remaining participants, used for the hallmark trait GWAS

The split is performed in `gwas_regenie/scripts/generate_inputs.R` and verified by
`gwas_regenie/scripts/preflight_checks.R`, which asserts that the two EID sets are disjoint.

### MR directions

| Direction | Instruments from | Outcome looked up in | Key scripts |
|-----------|------------------|----------------------|-------------|
| Forward | UKB-PPP *cis*-pQTLs | Held-out hallmark GWAS | `mr_twosample/scripts/twosampleMR_cis_trans_chunked.R`, `twosampleMR_proteins_{bmi,hba1c,trig_hdl_ratio_full}.R` |
| Reverse | Held-out hallmark GWAS | deCODE Iceland pQTLs | `mr_twosample/scripts/reverse_mr_bmi_to_all_proteins.R`, `reverse_mr_phenotype_to_all_proteins.R` |
| Forward replication | deCODE *cis*-pQTLs | Held-out hallmark GWAS | `mr_twosample/scripts/bidirectional_mr_decode_trig_hdl.R` |

Using deCODE as the reverse-MR outcome gives zero sample overlap with UK Biobank, and
restricting the forward analysis to *cis*-pQTLs limits horizontal pleiotropy.

### Sensitivity analysis

Reverse-MR instruments are distance-pruned at 500 kb in the primary analysis. A sensitivity
analysis repeats this with stricter LD clumping at r² < 0.01 against the 1000G EUR panel:

```bash
sbatch mr_twosample/scripts/submit_reverse_mr_r2clump.sh   # BMI and TRIG_HDL_RATIO, 16 chunks each
Rscript mr_twosample/scripts/combine_reverse_mr_r2clump.R  # combine chunks, compare against 500 kb results
```

A full end-to-end walkthrough — input generation, REGENIE submission, every MR direction,
and the resulting summary-statistic paths — is in
[`docs/genetics_pipeline/WALKTHROUGH_SPLIT_SAMPLE_MR.md`](docs/genetics_pipeline/WALKTHROUGH_SPLIT_SAMPLE_MR.md).

### Notes on running these scripts

These scripts were written for an HPC cluster with SLURM and contain absolute paths to
UK Biobank genotype and phenotype files; adapt the paths at the top of each script before
running. Individual-level UK Biobank data cannot be redistributed and is not included here.
Some scripts additionally support the glycemic strata, but those runs are exploratory and
were not used for the manuscript. LD score regression uses the official
[LDSC](https://github.com/bulik/ldsc) tool via
`gwas_regenie/scripts/regenie_to_ldsc_munge_input.py`. Scripts that query the OpenGWAS API
read a JWT from `$OPENGWAS_JWT` or a gitignored `.secrets/opengwas_jwt` file.

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


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Data Availability
- **UK Biobank Data**: The data from the UK Biobank that support the findings of this study are available upon application (https://www.ukbiobank.ac.uk/register-apply/).
- **MAP-D atlas summary statistics (Supplementary Data S1)**: https://doi.org/10.6084/m9.figshare.30007306 — the integrated atlas workbook (bidirectional observational associations, bidirectional Mendelian randomization, STEP 1/2 semaglutide comparisons, incident-CAD associations, DrugBank matches, GTEx tissue enrichment, and LDSC genetic correlations). figshare is the canonical citable copy; the same workbook is mirrored in this repository at [`app/SUPP_TABLE_1_05_20_2026.xlsx`](app/SUPP_TABLE_1_05_20_2026.xlsx) (the source the interactive atlas is built from) and its worksheets are documented in [`docs/supplementary_data_S1.md`](docs/supplementary_data_S1.md).
- **Hallmark GWAS summary statistics (BMI, HbA1c, TG/HDL)**: https://doi.org/10.6084/m9.figshare.32048550
- **MAP-D Web Resource (interactive atlas)**: https://stejat98-map-d.share.connect.posit.cloud/ (source in [`app/`](app/))
- **Code archive (Zenodo)**: https://doi.org/10.5281/zenodo.17071087
- **GitHub Repository**: https://github.com/stejat98/MAP-D

## Acknowledgments

- UK Biobank participants and investigators
- Olink Proteomics for platform development
- DrugBank database maintainers
- R Core Team and package developers


---

