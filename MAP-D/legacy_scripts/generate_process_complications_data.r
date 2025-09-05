## generate/process complication data script



## CKD processing

# also eGFR
# Load required libraries
library(dplyr)
library(tidyr)
library(readr)
library(lubridate)
library(fst)

ukb34521 <- read.fst("/n/no_backup2/patel/uk_biobank/main_data_34521/ukb34521.fst")



library(dplyr)
library(tidyr)
library(lubridate)

# Define ICD-10 disease codes
#hf_codes <- c("I50", "I110", "I130", "I132")
hf_codes <- c("I50", "I500", "I501", "I509",  # Heart failure
              "I110",                        # Hypertensive heart disease with HF
              "I130", "I132")                # HTN heart and renal disease with HF
nafld_codes <- c("K760", "K758")
# ckd_icd10_codes <- c("N180", "N181", "N182", "N183", "N184", "N185", "N189",
#                      "N18",   # Base CKD category
#                      "I120", "I131")  # Hypertensive kidney disease codes
ckd_codes <- c(
  "N18", "N181", "N182", "N183", "N184", "N185", "N189",
  "I12", "I120", "I121",
  "I13", "I131", "I132")

# Identify all column names for ICD-10 diagnoses and dates
dx_cols <- grep("^f.41270.0\\.", colnames(ukb34521), value = TRUE)
dt_cols <- gsub("41270", "41280", dx_cols)

# Reshape long format
dx_long <- ukb34521 %>%
  select(f.eid, all_of(dx_cols)) %>%
  pivot_longer(cols = -`f.eid`, names_to = "col", values_to = "icd10") %>%
  mutate(index = gsub(".*\\.", "", col)) %>%
  select(-col)

dt_long <- ukb34521 %>%
  select(f.eid, all_of(dt_cols)) %>%
  pivot_longer(cols = -`f.eid`, names_to = "col", values_to = "dx_date") %>%
  mutate(index = gsub(".*\\.", "", col)) %>%
  select(-col)

# Join diagnosis and date
dx_combined <- left_join(dx_long, dt_long, by = c("f.eid", "index")) %>%
  mutate(dx_date = ymd(dx_date))

# Get baseline assessment date
baseline_df <- ukb34521 %>%
  select(f.eid, baseline_date = `f.53.0.0`) %>%
  mutate(baseline_date = ymd(baseline_date))

# Function to filter incident cases
get_incident <- function(df, codes, condition) {
  df %>%
    filter(icd10 %in% codes) %>%
    left_join(baseline_df, by = "f.eid") %>%
    filter(!is.na(dx_date) & dx_date > baseline_date) %>%
    group_by(f.eid) %>%
    summarise(first_date = min(dx_date), condition = condition, source = "ICD10_EHR", .groups = "drop")
}

get_incident_ckd <- function(df, baseline_df, codes) {
  df %>%
    filter(icd10 %in% codes) %>%
    left_join(baseline_df, by = "f.eid") %>%
    filter(!is.na(dx_date) & dx_date > baseline_date) %>%
    group_by(f.eid) %>%
    summarise(
      event_date = min(dx_date),
      condition = "ckd",
      source = "ICD10_EHR",
      .groups = "drop"
    )
}

# Incident cases
hf_incident <- get_incident(dx_combined, hf_codes, "heart_failure")
nafld_incident <- get_incident(dx_combined, nafld_codes, "nafld")
ckd_incident <- get_incident_ckd(dx_combined, baseline_df, ckd_codes)

# Combine all
incident_dx <- bind_rows(hf_incident, nafld_incident, ckd_incident)


##########


library(dplyr)
library(tidyr)
library(lubridate)
library(purrr)

# Your ICD-10 definitions
hf_codes <- c("I50", "I500", "I501", "I509", "I110", "I130", "I132")
nafld_codes <- c("K760", "K758")
ckd_codes <- c("N18", "N181", "N182", "N183", "N184", "N185", "N189",
               "I12", "I120", "I121", "I13", "I131", "I132")

# Step 1: Master list of IDs from full dataset
all_ids <- dx_combined %>% select(f.eid) %>% distinct()

# Step 2: Get prevalent case function
get_prevalent <- function(dx_df, baseline_df, codes) {
  dx_df %>%
    filter(icd10 %in% codes) %>%
    left_join(baseline_df, by = "f.eid") %>%
    filter(!is.na(dx_date) & !is.na(baseline_date)) %>%
    filter(dx_date <= baseline_date) %>%
    distinct(f.eid)
}

# Step 3: Build binary variable (1 = incident, 0 = non-case, excluding prevalent)
build_binary_outcome <- function(incident_df, dx_df, baseline_df, codes, disease_label) {
  prevalent_ids <- get_prevalent(dx_df, baseline_df, codes)
  
  binary_df <- all_ids %>%
    left_join(incident_df %>% select(f.eid) %>% mutate(y = 1), by = "f.eid") %>%
    mutate(y = ifelse(is.na(y), 0, y)) %>%
    anti_join(prevalent_ids, by = "f.eid") %>%  # Remove prevalent cases
    rename(!!paste0("y_", disease_label) := y)
  
  return(binary_df)
}

# Step 4: Run for each condition
ckd_y    <- build_binary_outcome(ckd_incident,    dx_combined, baseline_df, ckd_codes,    "ckd")
hf_y     <- build_binary_outcome(hf_incident,     dx_combined, baseline_df, hf_codes,     "hf")
nafld_y  <- build_binary_outcome(nafld_incident,  dx_combined, baseline_df, nafld_codes,  "nafld")



# Step 5: Combine outcomes into one wide dataframe
outcomes_df <- reduce(list(ckd_y, hf_y, nafld_y), full_join, by = "f.eid")



# Step 6: Inspect
head(outcomes_df)
table(outcomes_df$y_ckd)
table(outcomes_df$y_hf)
table(outcomes_df$y_nafld)

library(fst)
data <- read.fst("/n/groups/patel/sivateja/UKB/all_cause_mort_data_exposures_full_plus_proteomics.fst")

data_plus_GLP_complications <- left_join(data,outcomes_df, by = c("eid" = "f.eid"))

write.fst(data_plus_GLP_complications, "/n/groups/patel/sivateja/UKB/all_cause_mort_data_exposures_full_plus_proteomics_plus_GLP_complications.fst")

library(fst)

library(tidyverse)


data_plus_GLP_complications <- read.fst("/n/groups/patel/sivateja/UKB/all_cause_mort_data_exposures_full_plus_proteomics_plus_GLP_complications.fst")

load("/home/st320/UKB_PEWAS/adjustments_survival_analysis.Rdata")

library(dplyr)
data_plus_GLP_complications  <- data_plus_GLP_complications  %>% filter(f.21000.0.0 %in% c("British","Irish", "Any other white background"))
adjustments <- c(adjustments[1:3], "x.738")


## 

library(tableone)

data_plus_GLP_complications <- data_plus_GLP_complications %>% rename(BMI = `21001_0.0.y`, HDL = `30760_0.0.y`,  LDL = `30780_0.0.y`, systolic_BP = `4080_0.0.y`, diastolic_BP = `4079_0.0.y`, ALT = `30620_0.0.y` )

data_plus_GLP_complications <- data_plus_GLP_complications %>% mutate(TRIG_HDL_RATIO = `30870_0.0.y`/HDL)

data_plus_GLP_complications <- data_plus_GLP_complications %>% rename(HbA1c = `30750_0.0.y`)

data <- data_plus_GLP_complications
# Create GlycemicStatus variable
data$GlycemicStatus <- with(data, ifelse(
  diabetes.y == 1, "Diabetes",
  ifelse(HbA1c >= 39 & HbA1c <= 48, "Prediabetes", "Normoglycemic")
))

# Convert it to a factor with preferred order
data$GlycemicStatus <- factor(data$GlycemicStatus, levels = c("Normoglycemic", "Prediabetes", "Diabetes"))


data_glycemic_status_HbA1c_adjusted <- data[-which(data$GlycemicStatus == "Normoglycemic" & data$HbA1c > 39),]

vars <- c(adjustments, "BMI", "HbA1c", "HDL", "LDL", "TRIG_HDL_RATIO", "systolic_BP", "diastolic_BP",
          "cad", "y_ckd", "y_nafld")

data_glycemic_status_HbA1c_adjusted$x.sex <- as.factor(as.character(data_glycemic_status_HbA1c_adjusted$x.sex))
data_glycemic_status_HbA1c_adjusted$x.738 <- as.factor(as.character(data_glycemic_status_HbA1c_adjusted$x.738))
data_glycemic_status_HbA1c_adjusted$cad <- as.factor(as.character(data_glycemic_status_HbA1c_adjusted$cad))
data_glycemic_status_HbA1c_adjusted$y_ckd <- as.factor(as.character(data_glycemic_status_HbA1c_adjusted$y_ckd))
data_glycemic_status_HbA1c_adjusted$y_nafld <- as.factor(as.character(data_glycemic_status_HbA1c_adjusted$y_nafld))

table1 <- CreateTableOne(vars = vars, strata = "GlycemicStatus", data = data_glycemic_status_HbA1c_adjusted, test = FALSE)

print(table1, showAllLevels = TRUE)

# Convert to data frame
table1_df <- print(table1, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)

# Save as CSV
write_csv(as.data.frame(table1_df), "/n/groups/patel/sivateja/UKB/pwas_metabolic_disease_hallmarks_table_one.csv")


library(dplyr)

# Assuming your dataframe is called baseline_df
median_ldl_t2d <-data_glycemic_status_HbA1c_adjusted %>%
  filter(GlycemicStatus == "Diabetes") %>%
  summarise(median_LDL = median(LDL, na.rm = TRUE))

print(median_ldl_t2d)

# Assuming your dataframe is called baseline_df
median_ldl_t2d <-data_glycemic_status_HbA1c_adjusted %>%
  filter(GlycemicStatus == "Prediabetes") %>%
  summarise(median_LDL = median(LDL, na.rm = TRUE))

print(median_ldl_t2d)

# Assuming your dataframe is called baseline_df
median_ldl_t2d <-data_glycemic_status_HbA1c_adjusted %>%
  filter(GlycemicStatus == "Normoglycemic") %>%
  summarise(median_LDL = median(LDL, na.rm = TRUE))

print(median_ldl_t2d)



