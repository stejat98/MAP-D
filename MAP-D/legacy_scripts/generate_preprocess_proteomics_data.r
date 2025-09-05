
## generate/process proteomics data



library(tidyverse)


olink_data_22881 <- read_tsv("/n/groups/patel/uk_biobank/olink_22881_52887/olink_data_22881.txt")


library(reshape2)

library(fst)

olink_data_22881_wide <- spread(olink_data_22881 , key =protein_id, value = result)

colnames(olink_data_22881_wide) <- sprintf("protein_x.%s", colnames(olink_data_22881_wide))
#olink_data_22881_wide <- reshape(olink_data_22881, idvar = "eid", timevar = "protein_id", direction = "wide")
saveRDS(olink_data_22881_wide,"/n/groups/patel/sivateja/olink_data_reshaped.RDS")

olink_data_22881_wide <- readRDS("/n/groups/patel/sivateja/olink_data_reshaped.RDS")

data <- read.fst("/n/groups/patel/sivateja/UKB/all_cause_mort_data_exposures_full.fst")

data <- left_join(data, olink_data_22881_wide, by = c("eid" = "protein_x.eid"))
#olink_data_22881_wide <- olink_data_22881 %>% dcast(protein_id ~ result)

write.fst(data, "/n/groups/patel/sivateja/UKB/all_cause_mort_data_exposures_full_plus_proteomics.fst")


#protein_vars <- colnames(olink_data_22881_wide)[-c(1,2)]

library(fst)
source("/home/st320/cardiometabolic_Linear_reg_INRT_Baseline_PEWAS_Functions_Script.R") 


load("/home/st320/UKB_PEWAS/adjustments_survival_analysis.Rdata")

data <- read.fst("/n/groups/patel/sivateja/UKB/all_cause_mort_data_exposures_full_plus_proteomics.fst")
data <- data %>% filter(f.21000.0.0 %in% c("British","Irish", "Any other white background"))
adjustments <- c(adjustments)
protein_vars <- colnames(data)[grep("protein_x.",colnames(data))][-1]




exposures <- protein_vars

#protein_vars <- gsub("protein_x.x.","protein_x.", protein_vars)

#protein_vars <- gsub("protein_x.protein_","protein_x.", protein_vars)

data <- data %>% rename(BMI = `21001_0.0.y`, HDL = `30760_0.0.y`,  LDL = `30780_0.0.y`, systolic_BP = `4080_0.0.y`, diastolic_BP = `4079_0.0.y`, ALT = `30620_0.0.y` )

data <- data %>% mutate(TRIG_HDL_RATIO = `30870_0.0.y`/HDL)

data <- data %>% rename(HbA1c = `30750_0.0.y`)

ukb34521 <- read.fst("/n/no_backup2/patel/uk_biobank/main_data_34521/ukb34521.fst")
ukb34521 <- ukb34521 %>% select(c("f.eid","f.74.0.0","f.23129.0.0", "f.23101.0.0"))

# ALST (kg) = (0.958 × [appendicular FFM (kg) [23129]]) − (0.166 × S) − 0.308, with S taking the value 0 if female and 1 if male.

## 




data <- left_join(data, ukb34521, by =c("eid"="f.eid"))

data <- data %>% mutate(ALST = (0.958 * `f.23129.0.0`) - (0.166 * `x.sex`) - 0.308)
  
# Create GlycemicStatus variable
data$GlycemicStatus <- with(data, ifelse(
  diabetes.y == 1, "Diabetes",
  ifelse(HbA1c >= 39 & HbA1c <= 48, "Prediabetes", "Normoglycemic")
))

# Convert it to a factor with preferred order
data$GlycemicStatus <- factor(data$GlycemicStatus, levels = c("Normoglycemic", "Prediabetes", "Diabetes"))

# Check how many samples in each category
table(data$GlycemicStatus)

data_Normo <- data %>% filter(GlycemicStatus == "Normoglycemic")
summary(data_Normo$HbA1c)

## filtered accounting for number of norm

data_glycemic_status_HbA1c_adjusted <- data[-which(data$GlycemicStatus == "Normoglycemic" & data$HbA1c > 39),]

#data_Normo_HbA1c_adjusted <- data_Normo %>% filter(HbA1c < 39)


demo_HbA1c_adjusted_proteomic_values_only <- data_glycemic_status_HbA1c_adjusted[, names(data_glycemic_status_HbA1c_adjusted) %in% c(protein_vars,"x.age","x.sex","GlycemicStatus","eid","BMI", "TRIG_HDL_RATIO", "HDL", "LDL", "systolic_BP", "diastolic_BP", "HbA1c", "cad", "x.738")]

demo_HbA1c_adjusted_proteomic_values_only <- as.data.frame(demo_HbA1c_adjusted_proteomic_values_only)

demo_HbA1c_adjusted_proteomic_values_only_proteomics_filtered <- demo_HbA1c_adjusted_proteomic_values_only %>% filter(!is.na(protein_x.1))

table(demo_HbA1c_adjusted_proteomic_values_only_proteomics_filtered$GlycemicStatus)

table(demo_HbA1c_adjusted_proteomic_values_only_proteomics_filtered$GlycemicStatus)/nrow(demo_HbA1c_adjusted_proteomic_values_only_proteomics_filtered)

table(demo_HbA1c_adjusted_proteomic_values_only_proteomics_filtered$x.sex)/nrow(demo_HbA1c_adjusted_proteomic_values_only_proteomics_filtered)


demo_HbA1c_adjusted_proteomic_values_only_proteomics_filtered %>% group_by(GlycemicStatus) %>% summarise(Median_BMI=median(BMI,na.rm=T),TRIG_HDL_RATIO=median(TRIG_HDL_RATIO,na.rm=T), HbA1c = median(HbA1c, na.rm=T), systolic_BP = median(systolic_BP, na.rm=T), HDL = median(HDL, na.rm=T), LDL = median(LDL, na.rm=T), diastolic_BP = median(diastolic_BP, na.rm=T)) %>% as.data.frame()
demo_HbA1c_adjusted_proteomic_values_only_proteomics_filtered %>% group_by(GlycemicStatus, x.738) %>% summarise(n = n(), .groups = "drop") %>%
  group_by(GlycemicStatus) %>%
  mutate(proportion = n / sum(n)) %>%
  arrange(GlycemicStatus, desc(proportion)) %>% as.data.frame()


demo_HbA1c_adjusted_proteomic_values_only_proteomics_filtered %>% group_by(GlycemicStatus, cad) %>% summarise(n = n(), .groups = "drop") %>%
  group_by(GlycemicStatus) %>%
  mutate(proportion = n / sum(n)) %>%
  arrange(GlycemicStatus, desc(proportion)) %>% as.data.frame()

# Keep only proteomic data for PCA
proteomic_values_only <- data[, names(data) %in% c(protein_vars,"GlycemicStatus","eid","BMI", "TRIG_HDL_RATIO", "HDL", "LDL", "systolic_BP", "diastolic_BP", "HbA1c", "cad")]

proteomic_values_only <- as.data.frame(proteomic_values_only)

proteomic_values_only_proteomics_filtered <- proteomic_values_only %>% filter(!is.na(protein_x.1))

#is_na_matrix <- is.na(olink_data_22881_wide)


# Count how many rows are incomplete *because* of each column
#col_na_impact <- colSums(is_na_matrix)

# Sort columns by number of incomplete cases they introduce
sort(col_na_impact, decreasing = TRUE)

proteomic_values_clean <- proteomic_values_only %>% select(-c("protein_x.1173", "protein_x.1889", "protein_x.1991"))


proteomic_values_complete <- na.omit(proteomic_values_clean)

Glycemic_status_complete_cases <- proteomic_values_complete %>% select(c("eid", "GlycemicStatus","BMI", "TRIG_HDL_RATIO", "HDL", "LDL", "systolic_BP", "diastolic_BP", "HbA1c","cad"))

proteomic_values_complete <- proteomic_values_complete %>% select(-c("eid", "GlycemicStatus","BMI", "TRIG_HDL_RATIO", "HDL", "LDL", "systolic_BP", "diastolic_BP", "HbA1c","cad"))
# PCA with scaling
pca_result <- prcomp(proteomic_values_complete, scale. = TRUE)

saveRDS(pca_result, "/n/groups/patel/sivateja/UKB/proteomic_all_overall_cad.RDS")
# Create a data frame with PCA results for plotting
pca_df <- data.frame(
  pca_result$x,
  proteomic_values_complete,
  Glycemic_status_complete_cases 
)

saveRDS(pca_df, "/n/groups/patel/sivateja/UKB/all_pca_risk_factors_raw_abundances_cad_df.RDS")
