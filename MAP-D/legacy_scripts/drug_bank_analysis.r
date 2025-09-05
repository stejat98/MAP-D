### updated with BMI HbA1c adjusted outcomes pwas 


library(xml2)
library(dplyr)
library(purrr)
library(stringr)
library(readr)

# Load incident outcome files
cad_df <- read_csv("/Users/sivatejatang/HMS Dropbox/Sivateja Tangirala/arpa/data_packet/PWAS_data_raw/CAD_logistic_regression_proteomic_glm_results_all_glycemic_group_adj_glyc_BMI_HbA1c_filtered.csv")
ckd_df <- read_csv("/Users/sivatejatang/HMS Dropbox/Sivateja Tangirala/arpa/data_packet/PWAS_data_raw/CKD_logistic_regression_proteomic_glm_results_all_glycemic_group_adj_glyc_BMI_HbA1c_filtered.csv")
nafld_df <- read_csv("/Users/sivatejatang/HMS Dropbox/Sivateja Tangirala/arpa/data_packet/PWAS_data_raw/NAFLD_logistic_regression_proteomic_glm_results_all_glycemic_group_adj_glyc_BMI_HbA1c_filtered.csv")

# Combine outcomes with gene_symbol extracted
incident_outcomes <- bind_rows(
  cad_df %>% mutate(outcome = "CAD"),
  ckd_df %>% mutate(outcome = "CKD"),
  nafld_df %>% mutate(outcome = "NAFLD")
) %>%
  rename(protein_raw = Protein) %>%
  mutate(
    gene_symbol = str_extract(protein_raw, "^[^;]+"),
    FDR = as.numeric(FDR),
    estimate = as.numeric(estimate) ,
    OR = exp(estimate)
  ) %>%
  filter(!is.na(FDR), FDR < 0.05)

# Categorize action needed based on beta estimate direction
incident_outcomes <- incident_outcomes %>%
  mutate(
    action_needed = case_when(
      estimate > 0 ~ "inhibit",   # protein up -> disease risk ↑ -> want to inhibit
      estimate < 0 ~ "activate",  # protein down -> disease risk ↓ -> want to activate
      TRUE ~ NA_character_
    )
  )

# Load DrugBank XML
xml_path <- "/Users/sivatejatang/HMS Dropbox/Sivateja Tangirala/arpa/data_packet/drugbank/full database.xml"
drugbank_xml <- read_xml(xml_path)
ns <- xml_ns_rename(xml_ns(drugbank_xml), d1 = "db")
drug_nodes <- xml_find_all(drugbank_xml, ".//db:drug", ns)

# Extract drug-target data from DrugBank
extract_drug_info <- function(drug_node) {
  groups <- xml_text(xml_find_all(drug_node, ".//db:groups/db:group", ns))
  if (!"approved" %in% groups) return(NULL)
  
  drug_id <- xml_text(xml_find_first(drug_node, ".//db:drugbank-id[@primary='true']", ns))
  drug_name <- xml_text(xml_find_first(drug_node, ".//db:name", ns))
  targets <- xml_find_all(drug_node, ".//db:targets/db:target", ns)
  
  map_dfr(targets, function(target_node) {
    known_action <- xml_text(xml_find_first(target_node, ".//db:known-action", ns))
    if (tolower(known_action) != "yes") return(NULL)
    
    actions <- xml_text(xml_find_all(target_node, ".//db:action", ns)) %>% tolower()
    if (length(actions) == 0) return(NULL)
    
    polypeptide <- xml_find_first(target_node, ".//db:polypeptide", ns)
    gene_symbol <- if (!is.na(polypeptide)) xml_text(xml_find_first(polypeptide, ".//db:gene-name", ns)) else NA
    protein_name <- if (!is.na(polypeptide)) xml_text(xml_find_first(polypeptide, ".//db:name", ns)) else NA
    
    tibble(
      drugbank_id = drug_id,
      drug_name = drug_name,
      gene_symbol = gene_symbol,
      protein_name = protein_name,
      actions = paste(actions, collapse = "; "),
      actions_list = list(actions)
    )
  })
}

# Run extraction
drug_df <- map_dfr(drug_nodes, extract_drug_info) %>%
  filter(!is.na(gene_symbol)) %>%
  distinct(drugbank_id, gene_symbol, .keep_all = TRUE)

# Define action type categories
inhibit_terms <- c("inhibitor", "antagonist", "blocker", "suppressor")
activate_terms <- c("agonist", "activator", "inducer", "stimulator")

# Join and filter directionally matched drugs
drug_hits <- inner_join(incident_outcomes, drug_df, by = "gene_symbol") %>%
  rowwise() %>%
  filter(
    (action_needed == "inhibit"  & any(actions_list %in% inhibit_terms)) |
      (action_needed == "activate" & any(actions_list %in% activate_terms))
  ) %>%
  ungroup()

# Save final hits
write_csv(drug_hits, "/Users/sivatejatang/Downloads/updated_directionally_matched_drug_hits.csv")



# Remove rows with missing drug or protein info
drug_hits <- drug_hits %>%
  filter(!is.na(drug_name), !is.na(protein_name), !is.na(Phenotype))

# 1. Unique drugs per outcome
unique_drugs_per_outcome <- drug_hits %>%
  group_by(outcome) %>%
  summarise(n_unique_drugs = n_distinct(drug_name)) %>%
  arrange(desc(n_unique_drugs))

# 2. Unique drug-protein pairs per outcome
unique_pairs_per_outcome <- drug_hits %>%
  group_by(outcome) %>%
  summarise(n_unique_drug_protein_pairs = n_distinct(paste(drug_name, protein_name, sep = "_"))) %>%
  arrange(desc(n_unique_drug_protein_pairs))

# 3. Overall counts
total_unique_drugs <- drug_hits %>% distinct(drug_name) %>% nrow()
total_unique_pairs <- drug_hits %>% distinct(drug_name, protein_name) %>% nrow()

# Print summary
cat("Total unique drugs across all outcomes:", total_unique_drugs, "\n")
cat("Total unique drug-protein pairs across all outcomes:", total_unique_pairs, "\n")

# View tables
print(unique_drugs_per_outcome)
print(unique_pairs_per_outcome)


