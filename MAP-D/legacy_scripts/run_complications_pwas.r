## run complications  script 

## NAFLD

library(tidyverse)


merged_validated_proteins <- read_rds("/n/groups/patel/sivateja/UKB/merged_validated_proteins_2.csv")
validated_proteins_unique_all_step1_step_2 <- unique(merged_validated_proteins$Exposure_UKB)

source("Baseline_PEWAS_Logistic_Functions_script.R")

#  adjust for fasting time
adjustments <- c(adjustments, "f.74.0.0")
ukb34521 <- read.fst("/n/no_backup2/patel/uk_biobank/main_data_34521/ukb34521.fst")
ukb34521 <- ukb34521 %>% select(c("f.eid","f.74.0.0"))

data_plus_GLP_complications <- left_join(data_plus_GLP_complications, ukb34521, by =c("eid"="f.eid"))


## run EWAS for specified phenotype (Baseline covariates) [extract p-values, etc.]
EWASLogistic(data=data_plus_GLP_complications, depvar = "y_nafld",
     adjustments = adjustments,exposures= validated_proteins_unique_all_step1_step_2,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/NAFLD_logistic_regression_proteomic_glm_results_all_glycemic_group"))


## run EWAS for specified phenotype (Baseline covariates) [extract p-values, etc.]
EWASLogistic(data=data_plus_GLP_complications, depvar = "y_ckd",
             adjustments = adjustments,exposures= validated_proteins_unique_all_step1_step_2,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/CKD_logistic_regression_proteomic_glm_results_all_glycemic_group"))

## run EWAS for specified phenotype (Baseline covariates) [extract p-values, etc.]
EWASLogistic(data=data_plus_GLP_complications, depvar = "y_hf",
             adjustments = adjustments,exposures= validated_proteins_unique_all_step1_step_2,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/HF_logistic_regression_proteomic_glm_results_all_glycemic_group"))


data_plus_GLP_complications <- data_plus_GLP_complications %>% rename(BMI = `21001_0.0.y`,HbA1c = `30750_0.0.y`)

# Create GlycemicStatus variable
data_plus_GLP_complications$GlycemicStatus <- with(data_plus_GLP_complications, ifelse(
  diabetes.y == 1, "Diabetes",
  ifelse(HbA1c >= 39 & HbA1c <= 48, "Prediabetes", "Normoglycemic")
))

# Convert it to a factor with preferred order
data_plus_GLP_complications$GlycemicStatus <- factor(data_plus_GLP_complications$GlycemicStatus, levels = c("Normoglycemic", "Prediabetes", "Diabetes"))

data_plus_GLP_complications_glycemic_status_HbA1c_adjusted <- data_plus_GLP_complications[-which(data_plus_GLP_complications$GlycemicStatus == "Normoglycemic" & data_plus_GLP_complications$HbA1c > 39),]

adjustments <- c(adjustments, "f.74.0.0", "GlycemicStatus", "BMI", "HbA1c")

## run EWAS for specified phenotype (Baseline covariates) [extract p-values, etc.]
EWASLogistic(data=data_plus_GLP_complications_glycemic_status_HbA1c_adjusted, depvar = "y_nafld",
             adjustments = adjustments,exposures= validated_proteins_unique_all_step1_step_2,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/NAFLD_logistic_regression_proteomic_glm_results_all_glycemic_group_adj_glyc_BMI_HbA1c"))

## run EWAS for specified phenotype (Baseline covariates) [extract p-values, etc.]
EWASLogistic(data=data_plus_GLP_complications_glycemic_status_HbA1c_adjusted, depvar = "cad",
             adjustments = adjustments,exposures= validated_proteins_unique_all_step1_step_2,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/CAD_logistic_regression_proteomic_glm_results_all_glycemic_group_adj_glyc_BMI_HbA1c"))


## run EWAS for specified phenotype (Baseline covariates) [extract p-values, etc.]
EWASLogistic(data=data_plus_GLP_complications_glycemic_status_HbA1c_adjusted, depvar = "y_ckd",
             adjustments = adjustments,exposures= validated_proteins_unique_all_step1_step_2,outFileName = sprintf("/n/groups/patel/sivateja/UKB/PEWAS_results/CKD_logistic_regression_proteomic_glm_results_all_glycemic_group_adj_glyc_BMI_HbA1c"))

