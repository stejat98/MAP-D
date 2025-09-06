# =============================================================================
# MAP-D Analysis Configuration File
# Metabolic Atlas of Progression to Diabetes (MAP-D)
# =============================================================================

# Project Information
PROJECT_NAME <- "MAP-D"
PROJECT_VERSION <- "1.0.0"
ANALYSIS_DATE <- Sys.Date()

# =============================================================================
# DATA PATHS
# =============================================================================

# Base data directory (modify these paths according to your system)
BASE_DATA_DIR <- ""
RAW_DATA_DIR <- ""
RESULTS_DIR <- "results"

# UK Biobank data files
UKB_MAIN_DATA <- file.path(RAW_DATA_DIR, "main_data_34521/ukb34521.fst")
OLINK_DATA <- file.path(RAW_DATA_DIR, "olink_22881_52887/olink_data_22881.txt")

# Processed data files
PROCESSED_DATA_FILE <- file.path(BASE_DATA_DIR, "all_cause_mort_data_exposures_full_plus_proteomics_plus_GLP_complications.fst")
PROTEOMIC_PCA_FILE <- file.path(BASE_DATA_DIR, "proteomic_all_overall_cad.RDS")
PCA_RISK_FACTORS_FILE <- file.path(BASE_DATA_DIR, "all_pca_risk_factors_raw_abundances_cad_df.RDS")

# External validation data
STEP1_VALIDATION_FILE <- "data/external/41591_2024_3355_MOESM3_ESM.xlsx"
DRUGBANK_XML_FILE <- "data/external/drugbank/full database.xml"

# PWAS Helper Function Scripts (existing validated functions)
PWAS_LINEAR_FUNCTIONS <- "src/cardiometabolic_Linear_reg_INRT_Baseline_PEWAS_Functions_Script.R"
PWAS_LOGISTIC_FUNCTIONS <- "src/Baseline_PEWAS_Logistic_Functions_script.R"

# =============================================================================
# ANALYSIS PARAMETERS
# =============================================================================

# Sample inclusion criteria
ANCESTRY_FILTER <- c("British", "Irish", "Any other white background")

# Glycemic status definitions
PREDIABETES_HBA1C_MIN <- 39  # mmol/mol
PREDIABETES_HBA1C_MAX <- 48  # mmol/mol
NORMOGLYCEMIC_HBA1C_MAX <- 39  # mmol/mol (for strict normoglycemic definition)

# Statistical thresholds
FDR_THRESHOLD <- 0.05
BONFERRONI_THRESHOLD <- 0.05
LASSO_TRAIN_PROPORTION <- 0.7

# ICD-10 codes for complications
HEART_FAILURE_CODES <- c("I50", "I500", "I501", "I509", "I110", "I130", "I132")
NAFLD_CODES <- c("K760", "K758")
CKD_CODES <- c("N18", "N181", "N182", "N183", "N184", "N185", "N189",
               "I12", "I120", "I121", "I13", "I131", "I132")

# =============================================================================
# COVARIATES AND ADJUSTMENTS
# =============================================================================

# Base demographic adjustments
BASE_ADJUSTMENTS <- c("x.age", "x.sex", "x.738")  # age, sex, assessment center

# Additional adjustments for specific analyses
FASTING_TIME_VAR <- "f.74.0.0"
GLYCEMIC_STATUS_VAR <- "GlycemicStatus"

# Clinical phenotypes for analysis
METABOLIC_PHENOTYPES <- c("BMI", "TRIG_HDL_RATIO", "HDL", "LDL", 
                         "systolic_BP", "diastolic_BP", "HbA1c")

# Complication outcomes
COMPLICATION_OUTCOMES <- c("y_ckd", "y_nafld", "y_hf", "cad")

# =============================================================================
# DRUG ANALYSIS PARAMETERS
# =============================================================================

# DrugBank action terms for drug-target matching
INHIBIT_TERMS <- c("inhibitor", "antagonist", "blocker", "suppressor")
ACTIVATE_TERMS <- c("agonist", "activator", "inducer", "stimulator")

# =============================================================================
# OUTPUT SETTINGS
# =============================================================================

# File naming conventions
RESULTS_PREFIX <- "MAP_D"
DATE_SUFFIX <- format(Sys.Date(), "%Y%m%d")

# Plot settings
PLOT_WIDTH <- 10
PLOT_HEIGHT <- 8
PLOT_DPI <- 300

# =============================================================================
# COMPUTATIONAL SETTINGS
# =============================================================================

# Random seed for reproducibility
RANDOM_SEED <- 123

# Parallel processing
N_CORES <- parallel::detectCores() - 1

# Memory settings
options(scipen = 999)  # Disable scientific notation
options(stringsAsFactors = FALSE)

# =============================================================================
# VALIDATION SETTINGS
# =============================================================================

# Minimum sample size requirements
MIN_SAMPLE_SIZE <- 100
MIN_CASES_FOR_ANALYSIS <- 50

# Quality control thresholds
MAX_MISSING_PROPORTION <- 0.1  # 10% maximum missing data
MIN_PROTEIN_DETECTION_RATE <- 0.8  # 80% minimum detection rate

cat("Configuration loaded successfully for MAP-D analysis\n")
cat("Project:", PROJECT_NAME, "Version:", PROJECT_VERSION, "\n")
cat("Analysis date:", as.character(ANALYSIS_DATE), "\n")
