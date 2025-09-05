##run metabolic hallmark phenotypes PWAS


## BMI

# stratify by diabetes

## BMI
## run EWAS for specified phenotype (Baseline covariates) [extract p-values, etc.]
EWAS(data=data_non_T2D, depvar = "BMI",
     adjustments = adjustments,exposures=exposures,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/BMI_Linear_regression_proteomic_lm_results_non_T2D"))


## run EWAS for specified phenotype (Baseline covariates) [extract p-values, etc.]
EWAS(data=data_T2D, depvar = "BMI",
     adjustments = adjustments,exposures=protein_vars,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/BMI_Linear_regression_proteomic_lm_results_T2D"))

## run EWAS for specified phenotype (Baseline covariates) [extract p-values, etc.]
EWAS(data=data_prediabetes, depvar = "BMI",
     adjustments = adjustments,exposures=protein_vars,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/BMI_Linear_regression_proteomic_lm_results_prediabetes"))


# stratify by diabetes (and adjust for fasting time)

adjustments <- c(adjustments, "f.74.0.0")
## BMI
## run EWAS for specified phenotype (Baseline covariates) [extract p-values, etc.]
EWAS(data=data_non_T2D, depvar = "BMI",
     adjustments = adjustments,exposures=exposures,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/BMI_Linear_regression_proteomic_lm_results_non_T2D_adj_fasting_time"))


## run EWAS for specified phenotype (Baseline covariates) [extract p-values, etc.]
EWAS(data=data_T2D, depvar = "BMI",
     adjustments = adjustments,exposures=protein_vars,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/BMI_Linear_regression_proteomic_lm_results_T2D_adj_fasting_time"))

## run EWAS for specified phenotype (Baseline covariates) [extract p-values, etc.]
EWAS(data=data_prediabetes, depvar = "BMI",
     adjustments = adjustments,exposures=protein_vars,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/BMI_Linear_regression_proteomic_lm_results_prediabetes_adj_fasting_time"))



## T2D

## non T2D


## Trig/HDL
## run EWAS for specified phenotype (Baseline covariates) [extract p-values, etc.]
EWAS(data=data_non_T2D, depvar = "TRIG_HDL_RATIO",
     adjustments = adjustments,exposures=exposures,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/TRIG_HDL_RATIO_Linear_regression_proteomic_lm_results_non_T2D"))

## run EWAS for specified phenotype (Baseline covariates) [extract p-values, etc.]
EWAS(data=data_T2D, depvar = "TRIG_HDL_RATIO",
     adjustments = adjustments,exposures=protein_vars,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/TRIG_HDL_RATIO_Linear_regression_proteomic_lm_results_T2D"))

## run EWAS for specified phenotype (Baseline covariates) [extract p-values, etc.]
EWAS(data=data_prediabetes, depvar = "TRIG_HDL_RATIO",
     adjustments = adjustments,exposures=protein_vars,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/TRIG_HDL_RATIO_Linear_regression_proteomic_lm_results_prediabetes"))


# stratify by diabetes (and adjust for fasting time)

## Trig/HDL
## run EWAS for specified phenotype (Baseline covariates) [extract p-values, etc.]
EWAS(data=data_non_T2D, depvar = "TRIG_HDL_RATIO",
     adjustments = adjustments,exposures=exposures,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/TRIG_HDL_RATIO_Linear_regression_proteomic_lm_results_non_T2D_adj_fasting_time"))

## run EWAS for specified phenotype (Baseline covariates) [extract p-values, etc.]
EWAS(data=data_T2D, depvar = "TRIG_HDL_RATIO",
     adjustments = adjustments,exposures=protein_vars,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/TRIG_HDL_RATIO_Linear_regression_proteomic_lm_results_T2D_adj_fasting_time"))

## run EWAS for specified phenotype (Baseline covariates) [extract p-values, etc.]
EWAS(data=data_prediabetes, depvar = "TRIG_HDL_RATIO",
     adjustments = adjustments,exposures=protein_vars,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/TRIG_HDL_RATIO_Linear_regression_proteomic_lm_results_prediabetes_adj_fasting_time"))





## HDL

EWAS(data=data_non_T2D, depvar = "HDL",
     adjustments = adjustments,exposures=exposures,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/HDL_Linear_regression_proteomic_lm_results_non_T2D"))

EWAS(data=data_T2D, depvar = "HDL",
     adjustments = adjustments,exposures=protein_vars,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/HDL_Linear_regression_proteomic_lm_results_T2D"))

## run EWAS for specified phenotype (Baseline covariates) [extract p-values, etc.]
EWAS(data=data_prediabetes, depvar = "HDL",
     adjustments = adjustments,exposures=protein_vars,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/HDL_Linear_regression_proteomic_lm_results_prediabetes"))


# adjust for fasting time

## HDL

EWAS(data=data_non_T2D, depvar = "HDL",
     adjustments = adjustments,exposures=exposures,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/HDL_Linear_regression_proteomic_lm_results_non_T2D_adj_fasting_time"))

EWAS(data=data_T2D, depvar = "HDL",
     adjustments = adjustments,exposures=protein_vars,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/HDL_Linear_regression_proteomic_lm_results_T2D_adj_fasting_time"))

## run EWAS for specified phenotype (Baseline covariates) [extract p-values, etc.]
EWAS(data=data_prediabetes, depvar = "HDL",
     adjustments = adjustments,exposures=protein_vars,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/HDL_Linear_regression_proteomic_lm_results_prediabetes_adj_fasting_time"))



## LDL


EWAS(data=data_non_T2D, depvar = "LDL",
     adjustments = adjustments,exposures=exposures,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/LDL_Linear_regression_proteomic_lm_results_non_T2D"))

EWAS(data=data_T2D, depvar = "LDL",
     adjustments = adjustments,exposures=protein_vars,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/LDL_Linear_regression_proteomic_lm_results_T2D"))

## run EWAS for specified phenotype (Baseline covariates) [extract p-values, etc.]
EWAS(data=data_prediabetes, depvar = "LDL",
     adjustments = adjustments,exposures=protein_vars,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/LDL_Linear_regression_proteomic_lm_results_prediabetes"))


# adjust for fasting time

## LDL

EWAS(data=data_non_T2D, depvar = "LDL",
     adjustments = adjustments,exposures=exposures,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/LDL_Linear_regression_proteomic_lm_results_non_T2D_adj_fasting_time"))

EWAS(data=data_T2D, depvar = "LDL",
     adjustments = adjustments,exposures=protein_vars,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/LDL_Linear_regression_proteomic_lm_results_T2D_adj_fasting_time"))

## run EWAS for specified phenotype (Baseline covariates) [extract p-values, etc.]
EWAS(data=data_prediabetes, depvar = "LDL",
     adjustments = adjustments,exposures=protein_vars,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/LDL_Linear_regression_proteomic_lm_results_prediabetes_adj_fasting_time"))



## Systolic BP

EWAS(data=data_non_T2D, depvar = "systolic_BP",
     adjustments = adjustments,exposures=exposures,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/systolic_BP_Linear_regression_proteomic_lm_results_non_T2D"))

EWAS(data=data_T2D, depvar = "systolic_BP",
     adjustments = adjustments,exposures=protein_vars,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/systolic_BP_Linear_regression_proteomic_lm_results_T2D"))

## run EWAS for specified phenotype (Baseline covariates) [extract p-values, etc.]
EWAS(data=data_prediabetes, depvar = "systolic_BP",
     adjustments = adjustments,exposures=protein_vars,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/systolic_BP_Linear_regression_proteomic_lm_results_prediabetes"))


# adjust for fasting time

## systolic_BP

EWAS(data=data_non_T2D, depvar = "systolic_BP",
     adjustments = adjustments,exposures=exposures,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/systolic_BP_Linear_regression_proteomic_lm_results_non_T2D_adj_fasting_time"))

EWAS(data=data_T2D, depvar = "systolic_BP",
     adjustments = adjustments,exposures=protein_vars,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/systolic_BP_Linear_regression_proteomic_lm_results_T2D_adj_fasting_time"))

## run EWAS for specified phenotype (Baseline covariates) [extract p-values, etc.]
EWAS(data=data_prediabetes, depvar = "systolic_BP",
     adjustments = adjustments,exposures=protein_vars,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/systolic_BP_Linear_regression_proteomic_lm_results_prediabetes_adj_fasting_time"))



## Diastolic BP

EWAS(data=data_non_T2D, depvar = "diastolic_BP",
     adjustments = adjustments,exposures=exposures,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/diastolic_BP_Linear_regression_proteomic_lm_results_non_T2D"))

EWAS(data=data_T2D, depvar = "diastolic_BP",
     adjustments = adjustments,exposures=protein_vars,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/diastolic_BP_Linear_regression_proteomic_lm_results_T2D"))

## run EWAS for specified phenotype (Baseline covariates) [extract p-values, etc.]
EWAS(data=data_prediabetes, depvar = "diastolic_BP",
     adjustments = adjustments,exposures=protein_vars,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/diastolic_BP_Linear_regression_proteomic_lm_results_prediabetes"))


# adjust for fasting time

## diastolic_BP

EWAS(data=data_non_T2D, depvar = "diastolic_BP",
     adjustments = adjustments,exposures=exposures,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/diastolic_BP_Linear_regression_proteomic_lm_results_non_T2D_adj_fasting_time"))

EWAS(data=data_T2D, depvar = "diastolic_BP",
     adjustments = adjustments,exposures=protein_vars,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/diastolic_BP_Linear_regression_proteomic_lm_results_T2D_adj_fasting_time"))

## run EWAS for specified phenotype (Baseline covariates) [extract p-values, etc.]
EWAS(data=data_prediabetes, depvar = "diastolic_BP",
     adjustments = adjustments,exposures=protein_vars,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/diastolic_BP_Linear_regression_proteomic_lm_results_prediabetes_adj_fasting_time"))


## HbA1c

EWAS(data=data_non_T2D, depvar = "HbA1c",
     adjustments = adjustments,exposures=exposures,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/HbA1c_Linear_regression_proteomic_lm_results_non_T2D"))

EWAS(data=data_T2D, depvar = "HbA1c",
     adjustments = adjustments,exposures=protein_vars,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/HbA1c_Linear_regression_proteomic_lm_results_T2D"))

EWAS(data=data_prediabetes, depvar = "HbA1c",
     adjustments = adjustments,exposures=protein_vars,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/HbA1c_Linear_regression_proteomic_lm_results_prediabetes"))

# adjust for fasting time

## HbA1c

EWAS(data=data_non_T2D, depvar = "HbA1c",
     adjustments = adjustments,exposures=exposures,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/HbA1c_Linear_regression_proteomic_lm_results_non_T2D_adj_fasting_time"))

EWAS(data=data_T2D, depvar = "HbA1c",
     adjustments = adjustments,exposures=protein_vars,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/HbA1c_Linear_regression_proteomic_lm_results_T2D_adj_fasting_time"))

## run EWAS for specified phenotype (Baseline covariates) [extract p-values, etc.]
EWAS(data=data_prediabetes, depvar = "HbA1c",
     adjustments = adjustments,exposures=protein_vars,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/HbA1c_Linear_regression_proteomic_lm_results_prediabetes_adj_fasting_time"))

