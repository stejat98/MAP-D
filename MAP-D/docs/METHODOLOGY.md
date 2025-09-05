# MAP-D Methodology

## Overview

The Metabolic Atlas of Progression to Diabetes (MAP-D) employs a comprehensive proteomics analysis framework to map the molecular progression from healthy metabolic states to Type 2 Diabetes (T2D). This document details the statistical methods, analytical approaches, and computational procedures used in the MAP-D study.

## Study Design

### Population

- **Source**: UK Biobank participants with available proteomics data
- **Sample Size**: Median 47,963 participants across analyses
- **Inclusion Criteria**: 
  - European ancestry (British, Irish, Any other white background)
  - Complete proteomics measurements
  - Baseline clinical measurements
- **Exclusion Criteria**:
  - Prevalent disease at baseline (for incident outcome analyses)
  - Missing key covariates or outcomes

### Proteomics Platform

- **Technology**: Olink Proximity Extension Assay (PEA)
- **Proteins**: 2,923 proteins measured
- **Quality Control**: 
  - Minimum detection rate: 80%
  - Maximum missing data: 10% per protein
  - Outlier detection and removal

## Analytical Framework

### 1. Data Preprocessing

#### Glycemic Status Classification

Participants were classified into three glycemic groups:

- **Normoglycemic**: No diabetes diagnosis AND HbA1c < 39 mmol/mol
- **Prediabetes**: No diabetes diagnosis AND HbA1c 39-48 mmol/mol  
- **Type 2 Diabetes**: Physician diagnosis of diabetes OR HbA1c ≥ 48 mmol/mol

#### Quality Control Procedures

1. **Sample-level QC**:
   - Minimum 80% protein detection rate per sample
   - Ancestry verification using genetic data
   - Outlier detection using Mahalanobis distance

2. **Protein-level QC**:
   - Minimum 80% detection rate across samples
   - Coefficient of variation assessment
   - Batch effect correction

3. **Phenotype QC**:
   - Biologically plausible ranges
   - Consistency checks across time points
   - Missing data patterns assessment

### 2. Protein-Wide Association Studies (PWAS)

#### Statistical Models

**Linear Models** (for continuous outcomes):
```
Outcome ~ Protein + Age + Sex + Assessment_Center + Fasting_Time + ε
```

**Logistic Models** (for binary outcomes):
```
logit(P(Outcome=1)) = β₀ + β₁×Protein + β₂×Age + β₃×Sex + β₄×Assessment_Center + β₅×Glycemic_Status + β₆×BMI + β₇×HbA1c
```

#### Stratified Analyses

Analyses were performed:
1. **Overall**: All participants combined
2. **By Glycemic Status**: Normoglycemic, Prediabetes, T2D separately
3. **By BMI Categories**: Normal weight, overweight, obese

#### Multiple Testing Correction

- **False Discovery Rate (FDR)**: Benjamini-Hochberg procedure
- **Family-wise Error Rate (FWER)**: Bonferroni correction
- **Significance Thresholds**: 
  - FDR < 0.05 for discovery
  - Bonferroni < 0.05 for stringent significance

### 3. Machine Learning Approaches

#### LASSO Regression

**Objective**: Identify parsimonious sets of predictive proteins

**Model Specification**:
```
min(β) { (1/2n) ||y - Xβ||²₂ + λ||β||₁ }
```

Where:
- `y`: Outcome vector
- `X`: Predictor matrix (proteins + covariates)
- `λ`: Regularization parameter (selected via cross-validation)
- `β`: Coefficient vector

**Cross-Validation Strategy**:
- 10-fold cross-validation for λ selection
- Stratified 70/30 train-test split
- Glycemic status-stratified sampling

**Model Comparison**:
1. **Baseline Model**: Demographics + clinical variables only
2. **Protein Model**: Proteins only
3. **Combined Model**: Demographics + clinical + proteins

**Performance Metrics**:
- R² on held-out test set
- Root Mean Square Error (RMSE)
- Mean Absolute Error (MAE)

### 4. Drug Target Analysis

#### DrugBank Integration

**Data Source**: DrugBank database (version 5.1.8)
**Inclusion Criteria**:
- Approved drugs only
- Known mechanism of action
- Human protein targets

#### Directional Matching Algorithm

For each significant protein association:

1. **Effect Direction**: Determine if protein is positively or negatively associated with outcome
2. **Therapeutic Action**: Classify drug action as inhibitory or activating
3. **Directional Consistency**: Match proteins requiring inhibition with inhibitory drugs, and proteins requiring activation with activating drugs

**Action Classification**:
- **Inhibitory**: antagonist, inhibitor, blocker, suppressor
- **Activating**: agonist, activator, inducer, stimulator

#### Repurposing Opportunities

**Multi-target Analysis**: Drugs targeting multiple significant proteins
**Multi-phenotype Analysis**: Proteins associated with multiple metabolic outcomes
**Therapeutic Intransience**: Proteins associated with outcomes but not affected by existing therapies

### 5. External Validation

#### Validation Datasets

**STEP1 Dataset**: 
- Source: Semaglutide intervention trials
- Proteins: Overlapping set with MAP-D
- Validation criterion: Consistent effect direction and significance

**STEP2 Dataset**:
- Source: Independent proteomics cohorts
- Design: Cross-sectional and longitudinal
- Validation criterion: Replication of associations

#### Validation Metrics

- **Replication Rate**: Proportion of significant associations validated
- **Effect Concordance**: Correlation of effect sizes
- **Directional Consistency**: Agreement in effect direction

### 6. Statistical Considerations

#### Power Calculations

Sample size requirements were calculated assuming:
- Type I error rate: α = 0.05
- Type II error rate: β = 0.20 (80% power)
- Effect sizes: Cohen's d ≥ 0.1 for continuous outcomes, OR ≥ 1.2 for binary outcomes

#### Missing Data Handling

- **Complete Case Analysis**: Primary approach for main analyses
- **Multiple Imputation**: Sensitivity analyses for key findings
- **Missing Data Patterns**: Assessed for systematic bias

#### Sensitivity Analyses

1. **Ancestry Stratification**: European ancestry subgroups
2. **Medication Effects**: Exclusion of participants on relevant medications
3. **Fasting Status**: Stratification by fasting time
4. **Batch Effects**: Technical replicate analysis

## Computational Implementation

### Software and Packages

- **R Version**: ≥ 4.0.0
- **Key Packages**:
  - `glmnet`: LASSO regression
  - `tidyverse`: Data manipulation
  - `fst`: Fast data I/O
  - `xml2`: DrugBank parsing
  - `parallel`: Parallel processing

### Reproducibility Measures

- **Random Seeds**: Fixed for all stochastic procedures
- **Session Information**: Captured for all analyses
- **Version Control**: Git-based tracking of all code changes
- **Container Environment**: Docker images for computational environment

### Quality Assurance

- **Unit Testing**: Automated tests for all functions
- **Integration Testing**: End-to-end pipeline validation
- **Code Review**: Peer review of all analytical code
- **Result Validation**: Independent verification of key findings

## Limitations and Assumptions

### Study Limitations

1. **Cross-sectional Design**: Cannot establish causality for most associations
2. **Single Time Point**: Proteomics measured at baseline only
3. **Population Homogeneity**: Primarily European ancestry
4. **Platform Specificity**: Olink-specific protein measurements

### Statistical Assumptions

1. **Linear Relationships**: Assumed for continuous outcomes
2. **Independence**: Observations assumed independent
3. **Normality**: Residuals assumed normally distributed
4. **Homoscedasticity**: Constant variance assumed

### Biological Assumptions

1. **Protein Stability**: Plasma proteins reflect tissue-level biology
2. **Causal Relationships**: Protein changes drive phenotypic changes
3. **Therapeutic Targets**: Druggable proteins represent viable targets

## Future Directions

### Methodological Enhancements

1. **Longitudinal Analysis**: Multi-time point protein measurements
2. **Causal Inference**: Mendelian randomization approaches
3. **Multi-omics Integration**: Genomics, metabolomics, transcriptomics
4. **Machine Learning**: Advanced ML approaches (random forests, neural networks)

### Clinical Translation

1. **Biomarker Validation**: Prospective validation studies
2. **Risk Stratification**: Clinical risk prediction models
3. **Therapeutic Development**: Drug target prioritization
4. **Precision Medicine**: Personalized treatment strategies

---

This methodology represents the comprehensive analytical framework underlying the MAP-D study, designed to ensure reproducible, robust, and clinically relevant insights into the proteomics of diabetes progression.
