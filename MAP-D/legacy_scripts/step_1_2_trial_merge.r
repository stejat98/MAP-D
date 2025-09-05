## step 1 and step 2 integration script


# Load required packages
library(readxl)
library(dplyr)
library(purrr)

# Define the primary file and sheet
primary_file <- "/Users/sivatejatang/Downloads/41591_2024_3355_MOESM3_ESM.xlsx"  # Change this to your filename
primary_sheet <- "S2_tx_STEP1"  # Specify the relevant sheet name
primary_column <- "EntrezGeneSymbol"  # Adjust to match your dataset


# Read the primary protein list
primary_data <- read_excel(primary_file, sheet = primary_sheet)

# Define a list of PWAS files (modify with actual filenames)
pwas_files <- list.files("/Users/sivatejatang/HMS Dropbox/Sivateja Tangirala/arpa/data_packet/PWAS_data_raw", full.names=T)[-7]

# Function to read a specific sheet and extract protein column
read_pwas_data <- function(file, col_name="Protein") {
  pwas_data <- readr::read_csv(file) 
  
  pwas_data_formatted <- tidyr::separate(pwas_data, Protein, into = c("code", "name"), sep = ";", fill = "right", extra = "drop")
  
  return(pwas_data_formatted)
}

# Compute intersections
results <- map(pwas_files, ~ {
  pwas_data <- read_pwas_data(.x)
  pwas_data$FDR <- p.adjust(pwas_data$p.value, method="fdr")
  pwas_data$Bonferroni <- p.adjust(pwas_data$p.value, method="bonferroni")
  merged_step1_pwas_results <- left_join(pwas_data, primary_data, by = c("code" = "EntrezGeneSymbol"))
  
  return(merged_step1_pwas_results)
})

# Combine results into a single dataframe
STEP1_merged_results <- bind_rows(results)

# Define the primary file and sheet
primary_STEP2_file <- "/Users/sivatejatang/Downloads/41591_2024_3355_MOESM3_ESM.xlsx"  # Change this to your filename
primary_STEP2_sheet <- "S3_tx_STEP2"  # Specify the relevant sheet name
primary_STEP2_column <- "EntrezGeneSymbol"  # Adjust to match your dataset

# Read the primary protein list
primary_STEP2_data <- read_excel(primary_STEP2_file, sheet = primary_STEP2_sheet)

# Compute intersections
STEP2_results <- map(pwas_files, ~ {
  pwas_data <- read_pwas_data(.x)
  pwas_data$FDR <- p.adjust(pwas_data$p.value, method="fdr")
  pwas_data$Bonferroni <- p.adjust(pwas_data$p.value, method="bonferroni")
  merged_step1_pwas_results <- left_join(pwas_data, primary_STEP2_data, by = c("code" = "EntrezGeneSymbol"))
  
  return(merged_step1_pwas_results)
})

# Combine results into a single dataframe
STEP2_merged_results <- bind_rows(STEP2_results)


STEP1_validated_stats <- STEP1_merged_results %>% filter(Bonferroni < 0.05 & qvalue < 0.05 & ((estimate < 0 & effect_size < 0) | (estimate > 0 & effect_size > 0))) %>% group_by(Phenotype,Subgroup) %>% summarise(NumValidated=n())


STEP1_Total_stats <- STEP1_merged_results %>% filter(Bonferroni < 0.05) %>% group_by(Phenotype, Subgroup) %>% summarise(NumTotal=n())


STEP1_validated_stats <- STEP1_validated_stats %>% mutate(Validated_Proportion = NumValidated/STEP1_Total_stats$NumTotal)





STEP2_validated_stats <- STEP2_merged_results %>% filter(Bonferroni < 0.05 & qvalue < 0.05 & ((estimate < 0 & effect_size < 0) | (estimate > 0 & effect_size > 0))) %>% group_by(Phenotype,Subgroup) %>% summarise(NumValidated=n())


STEP2_Total_stats <- STEP2_merged_results %>% filter(Bonferroni < 0.05) %>% group_by(Phenotype, Subgroup) %>% summarise(NumTotal=n())


STEP2_validated_stats <- STEP2_validated_stats %>% mutate(Validated_Proportion = NumValidated/STEP2_Total_stats$NumTotal)



readr::write_csv(STEP1_merged_results, "/Users/sivatejatang/HMS Dropbox/Sivateja Tangirala/arpa/data_packet/STEP1_merged_results.csv")

readr::write_csv(STEP2_merged_results, "/Users/sivatejatang/HMS Dropbox/Sivateja Tangirala/arpa/data_packet/STEP2_merged_results.csv")

## treating T2D as a subgroup
## interaction of protein with T2D status ??
## significant protein interactions


STEP1_merged_results %>% filter(Bonferroni < 0.05 & qvalue < 0.05 & ((estimate < 0 & effect_size < 0) | (estimate > 0 & effect_size > 0))) %>% arrange(Bonferroni) %>% group_by(Phenotype,Subgroup) %>% slice(1:2)




STEP1_merged_results %>% filter(Bonferroni < 0.05 & qvalue < 0.05 & ((estimate < 0 & effect_size < 0) | (estimate > 0 & effect_size > 0))) %>% arrange(Bonferroni) %>% group_by(Phenotype,Subgroup) %>% slice(1:2) %>% as.data.frame() %>% 
  readr::write_csv("/Users/sivatejatang/HMS Dropbox/Sivateja Tangirala/arpa/data_packet/STEP1_top_val_hits.csv")
  

STEP2_merged_results %>% filter(Bonferroni < 0.05 & qvalue < 0.05 & ((estimate < 0 & effect_size < 0) | (estimate > 0 & effect_size > 0))) %>% arrange(Bonferroni) %>% group_by(Phenotype,Subgroup) %>% slice(1:2) %>% as.data.frame() %>% 
  readr::write_csv("/Users/sivatejatang/HMS Dropbox/Sivateja Tangirala/arpa/data_packet/STEP2_top_val_hits.csv")
