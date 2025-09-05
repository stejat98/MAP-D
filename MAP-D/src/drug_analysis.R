# =============================================================================
# MAP-D Drug Target Analysis Module
# Metabolic Atlas of Progression to Diabetes (MAP-D)
# =============================================================================

# Source required files
source("config.R")
source("src/utils.R")

#' Load and parse DrugBank XML database
#' 
#' @param drugbank_file Path to DrugBank XML file
#' @return Data frame with drug-target information
load_drugbank_data <- function(drugbank_file = DRUGBANK_XML_FILE) {
  
  log_analysis_step("Loading DrugBank Database", 
                   "Parsing XML file for approved drugs and targets")
  
  # Load required packages
  load_required_packages(c("xml2", "dplyr", "purrr", "stringr"))
  
  # Check if file exists
  if (!check_file_exists(drugbank_file, "DrugBank XML file")) {
    stop("DrugBank XML file not found. Please download from https://go.drugbank.com/releases/latest")
  }
  
  cat("Loading DrugBank XML file:", drugbank_file, "\n")
  cat("This may take several minutes...\n")
  
  # Load XML
  drugbank_xml <- read_xml(drugbank_file)
  ns <- xml_ns_rename(xml_ns(drugbank_xml), d1 = "db")
  
  # Find all drug nodes
  drug_nodes <- xml_find_all(drugbank_xml, ".//db:drug", ns)
  cat("Found", length(drug_nodes), "drugs in database\n")
  
  # Function to extract drug information
  extract_drug_info <- function(drug_node) {
    
    # Only include approved drugs
    groups <- xml_text(xml_find_all(drug_node, ".//db:groups/db:group", ns))
    if (!"approved" %in% groups) {
      return(NULL)
    }
    
    # Extract basic drug information
    drug_id <- xml_text(xml_find_first(drug_node, ".//db:drugbank-id[@primary='true']", ns))
    drug_name <- xml_text(xml_find_first(drug_node, ".//db:name", ns))
    
    # Find targets
    targets <- xml_find_all(drug_node, ".//db:targets/db:target", ns)
    
    if (length(targets) == 0) {
      return(NULL)
    }
    
    # Extract target information
    target_info <- map_dfr(targets, function(target_node) {
      
      # Only include targets with known action
      known_action <- xml_text(xml_find_first(target_node, ".//db:known-action", ns))
      if (tolower(known_action) != "yes") {
        return(NULL)
      }
      
      # Extract actions
      actions <- xml_text(xml_find_all(target_node, ".//db:action", ns))
      actions <- tolower(actions)
      
      if (length(actions) == 0) {
        return(NULL)
      }
      
      # Extract target protein information
      polypeptide <- xml_find_first(target_node, ".//db:polypeptide", ns)
      
      if (is.na(polypeptide)) {
        return(NULL)
      }
      
      gene_symbol <- xml_text(xml_find_first(polypeptide, ".//db:gene-name", ns))
      protein_name <- xml_text(xml_find_first(polypeptide, ".//db:name", ns))
      uniprot_id <- xml_text(xml_find_first(polypeptide, ".//db:external-identifier[db:resource='UniProtKB']/db:identifier", ns))
      
      return(data.frame(
        drugbank_id = drug_id,
        drug_name = drug_name,
        gene_symbol = gene_symbol,
        protein_name = protein_name,
        uniprot_id = uniprot_id,
        actions = paste(actions, collapse = "; "),
        actions_list = I(list(actions)),
        stringsAsFactors = FALSE
      ))
    })
    
    return(target_info)
  }
  
  # Extract information for all drugs (with progress)
  cat("Extracting drug-target information...\n")
  
  # Process in batches to show progress
  batch_size <- 1000
  n_batches <- ceiling(length(drug_nodes) / batch_size)
  
  all_drug_data <- list()
  
  for (i in seq_len(n_batches)) {
    start_idx <- (i - 1) * batch_size + 1
    end_idx <- min(i * batch_size, length(drug_nodes))
    
    batch_nodes <- drug_nodes[start_idx:end_idx]
    batch_data <- map_dfr(batch_nodes, extract_drug_info)
    
    if (nrow(batch_data) > 0) {
      all_drug_data[[i]] <- batch_data
    }
    
    if (i %% 5 == 0 || i == n_batches) {
      cat("  Processed", end_idx, "of", length(drug_nodes), "drugs\n")
    }
  }
  
  # Combine all data
  drug_data <- bind_rows(all_drug_data)
  
  # Remove duplicates and clean data
  drug_data <- drug_data %>%
    filter(!is.na(gene_symbol), gene_symbol != "") %>%
    distinct(drugbank_id, gene_symbol, .keep_all = TRUE)
  
  cat("Extracted", nrow(drug_data), "approved drug-target pairs\n")
  cat("Covering", length(unique(drug_data$drug_name)), "unique drugs\n")
  cat("Covering", length(unique(drug_data$gene_symbol)), "unique target genes\n")
  
  return(drug_data)
}

#' Classify drug actions as inhibitory or activating
#' 
#' @param drug_data Drug-target data frame from load_drugbank_data
#' @param inhibit_terms Vector of terms indicating inhibitory action
#' @param activate_terms Vector of terms indicating activating action
#' @return Drug data with action classification
classify_drug_actions <- function(drug_data, 
                                 inhibit_terms = INHIBIT_TERMS,
                                 activate_terms = ACTIVATE_TERMS) {
  
  log_analysis_step("Classifying Drug Actions", 
                   "Categorizing drugs as inhibitory or activating")
  
  # Function to classify actions for a single drug-target pair
  classify_single_action <- function(actions_list) {
    
    if (is.null(actions_list) || length(actions_list) == 0) {
      return("unknown")
    }
    
    actions <- unlist(actions_list)
    
    # Check for inhibitory terms
    is_inhibitory <- any(sapply(inhibit_terms, function(term) {
      any(grepl(term, actions, ignore.case = TRUE))
    }))
    
    # Check for activating terms
    is_activating <- any(sapply(activate_terms, function(term) {
      any(grepl(term, actions, ignore.case = TRUE))
    }))
    
    # Classify based on matches
    if (is_inhibitory && !is_activating) {
      return("inhibitory")
    } else if (is_activating && !is_inhibitory) {
      return("activating")
    } else if (is_inhibitory && is_activating) {
      return("mixed")
    } else {
      return("other")
    }
  }
  
  # Apply classification
  drug_data_classified <- drug_data %>%
    rowwise() %>%
    mutate(action_classification = classify_single_action(actions_list)) %>%
    ungroup()
  
  # Summary of classifications
  classification_summary <- table(drug_data_classified$action_classification)
  cat("Action classification summary:\n")
  print(classification_summary)
  
  return(drug_data_classified)
}

#' Match drugs to PWAS results based on directional consistency
#' 
#' @param pwas_results PWAS results data frame
#' @param drug_data Classified drug-target data
#' @param significance_threshold FDR threshold for significant associations
#' @return Data frame with matched drug-protein pairs
match_drugs_to_pwas <- function(pwas_results, drug_data, 
                               significance_threshold = FDR_THRESHOLD) {
  
  log_analysis_step("Matching Drugs to PWAS Results", 
                   "Finding directionally consistent drug-target pairs")
  
  # Filter for significant PWAS results
  significant_pwas <- pwas_results %>%
    filter(FDR < significance_threshold, !is.na(estimate))
  
  cat("Significant PWAS associations:", nrow(significant_pwas), "\n")
  
  if (nrow(significant_pwas) == 0) {
    warning("No significant PWAS results found")
    return(data.frame())
  }
  
  # Extract gene symbols from protein names
  significant_pwas <- significant_pwas %>%
    mutate(gene_symbol = str_extract(Protein, "^[^;]+"))
  
  # Determine required action based on effect direction
  significant_pwas <- significant_pwas %>%
    mutate(
      action_needed = case_when(
        estimate > 0 ~ "inhibit",   # protein up -> disease risk up -> need to inhibit
        estimate < 0 ~ "activate",  # protein down -> disease risk up -> need to activate
        TRUE ~ NA_character_
      )
    )
  
  # Join with drug data
  drug_matches <- inner_join(
    significant_pwas, 
    drug_data, 
    by = "gene_symbol"
  )
  
  cat("Initial drug-protein matches:", nrow(drug_matches), "\n")
  
  # Filter for directionally consistent matches
  consistent_matches <- drug_matches %>%
    filter(
      (action_needed == "inhibit" & action_classification == "inhibitory") |
      (action_needed == "activate" & action_classification == "activating")
    )
  
  cat("Directionally consistent matches:", nrow(consistent_matches), "\n")
  
  # Add additional information
  consistent_matches <- consistent_matches %>%
    mutate(
      OR = exp(estimate),  # Convert log odds to odds ratio for logistic models
      match_type = paste(action_needed, "with", action_classification, "drug")
    ) %>%
    arrange(FDR, p.value)
  
  return(consistent_matches)
}

#' Analyze drug repurposing opportunities
#' 
#' @param drug_matches Matched drug-protein pairs from match_drugs_to_pwas
#' @return List with repurposing analysis results
analyze_drug_repurposing <- function(drug_matches) {
  
  log_analysis_step("Drug Repurposing Analysis", 
                   "Identifying therapeutic opportunities")
  
  if (nrow(drug_matches) == 0) {
    warning("No drug matches provided")
    return(list())
  }
  
  # Summary by outcome
  outcome_summary <- drug_matches %>%
    group_by(Phenotype) %>%
    summarise(
      n_unique_drugs = n_distinct(drug_name),
      n_unique_proteins = n_distinct(gene_symbol),
      n_drug_protein_pairs = n(),
      top_drug = drug_name[which.min(FDR)][1],
      top_protein = gene_symbol[which.min(FDR)][1],
      best_fdr = min(FDR, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_unique_drugs))
  
  # Summary by drug
  drug_summary <- drug_matches %>%
    group_by(drug_name, drugbank_id) %>%
    summarise(
      n_targets = n_distinct(gene_symbol),
      n_phenotypes = n_distinct(Phenotype),
      phenotypes = paste(unique(Phenotype), collapse = "; "),
      targets = paste(unique(gene_symbol), collapse = "; "),
      best_fdr = min(FDR, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_targets), best_fdr)
  
  # Summary by protein target
  protein_summary <- drug_matches %>%
    group_by(gene_symbol, protein_name) %>%
    summarise(
      n_drugs = n_distinct(drug_name),
      n_phenotypes = n_distinct(Phenotype),
      drugs = paste(unique(drug_name), collapse = "; "),
      phenotypes = paste(unique(Phenotype), collapse = "; "),
      best_fdr = min(FDR, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_drugs), best_fdr)
  
  # Identify multi-target opportunities
  multi_target_drugs <- drug_summary %>%
    filter(n_targets > 1) %>%
    arrange(desc(n_targets))
  
  # Identify multi-phenotype opportunities  
  multi_phenotype_proteins <- protein_summary %>%
    filter(n_phenotypes > 1) %>%
    arrange(desc(n_phenotypes))
  
  cat("Repurposing analysis summary:\n")
  cat("  Total unique drugs:", length(unique(drug_matches$drug_name)), "\n")
  cat("  Total unique proteins:", length(unique(drug_matches$gene_symbol)), "\n")
  cat("  Multi-target drugs:", nrow(multi_target_drugs), "\n")
  cat("  Multi-phenotype proteins:", nrow(multi_phenotype_proteins), "\n")
  
  return(list(
    outcome_summary = outcome_summary,
    drug_summary = drug_summary,
    protein_summary = protein_summary,
    multi_target_drugs = multi_target_drugs,
    multi_phenotype_proteins = multi_phenotype_proteins
  ))
}

#' Run comprehensive drug target analysis
#' 
#' @param pwas_results PWAS results data frame
#' @param drugbank_file Path to DrugBank XML file
#' @param output_dir Directory to save results
#' @return List with all analysis results
run_drug_target_analysis <- function(pwas_results, 
                                   drugbank_file = DRUGBANK_XML_FILE,
                                   output_dir = RESULTS_DIR) {
  
  log_analysis_step("Comprehensive Drug Target Analysis", 
                   "Analyzing therapeutic targets from PWAS results")
  
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Step 1: Load DrugBank data
  drug_data <- load_drugbank_data(drugbank_file)
  
  # Step 2: Classify drug actions
  drug_data_classified <- classify_drug_actions(drug_data)
  
  # Step 3: Match drugs to PWAS results
  drug_matches <- match_drugs_to_pwas(pwas_results, drug_data_classified)
  
  # Step 4: Analyze repurposing opportunities
  repurposing_analysis <- analyze_drug_repurposing(drug_matches)
  
  # Save results
  cat("\nSaving drug analysis results...\n")
  
  # Save main results
  write_csv(drug_data_classified, file.path(output_dir, "drugbank_processed.csv"))
  write_csv(drug_matches, file.path(output_dir, "drug_target_matches.csv"))
  
  # Save summary analyses
  if (length(repurposing_analysis) > 0) {
    write_csv(repurposing_analysis$outcome_summary, 
             file.path(output_dir, "drug_repurposing_by_outcome.csv"))
    write_csv(repurposing_analysis$drug_summary, 
             file.path(output_dir, "drug_repurposing_by_drug.csv"))
    write_csv(repurposing_analysis$protein_summary, 
             file.path(output_dir, "drug_repurposing_by_protein.csv"))
  }
  
  # Save complete results as RDS
  complete_results <- list(
    drugbank_data = drug_data_classified,
    drug_matches = drug_matches,
    repurposing_analysis = repurposing_analysis
  )
  
  saveRDS(complete_results, file.path(output_dir, "drug_analysis_complete_results.rds"))
  
  cat("Drug target analysis completed successfully.\n")
  
  return(complete_results)
}

#' Create drug target network visualization
#' 
#' @param drug_matches Drug-protein matches data frame
#' @param top_n Number of top associations to include
#' @param output_file Output file path (optional)
#' @return Network plot object
create_drug_target_network <- function(drug_matches, top_n = 50, output_file = NULL) {
  
  load_required_packages(c("ggplot2", "ggraph", "igraph", "dplyr"))
  
  # Select top associations
  top_matches <- drug_matches %>%
    arrange(FDR) %>%
    head(top_n)
  
  if (nrow(top_matches) == 0) {
    warning("No drug matches to plot")
    return(NULL)
  }
  
  # Create network edges
  edges <- top_matches %>%
    select(from = drug_name, to = gene_symbol, 
           weight = -log10(FDR), phenotype = Phenotype) %>%
    mutate(edge_type = "drug_target")
  
  # Create nodes
  drug_nodes <- data.frame(
    name = unique(edges$from),
    type = "drug",
    stringsAsFactors = FALSE
  )
  
  protein_nodes <- data.frame(
    name = unique(edges$to),
    type = "protein",
    stringsAsFactors = FALSE
  )
  
  nodes <- bind_rows(drug_nodes, protein_nodes)
  
  # Create igraph object
  g <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
  
  # Create plot
  p <- ggraph(g, layout = "fr") +
    geom_edge_link(aes(alpha = weight), color = "gray60") +
    geom_node_point(aes(color = type, size = type)) +
    geom_node_text(aes(label = name), size = 3, repel = TRUE) +
    scale_color_manual(values = c("drug" = "#1f77b4", "protein" = "#ff7f0e")) +
    scale_size_manual(values = c("drug" = 3, "protein" = 4)) +
    theme_void() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(
      title = "Drug-Target Network",
      subtitle = paste("Top", top_n, "associations by FDR"),
      color = "Node Type",
      size = "Node Type"
    )
  
  # Save plot if output file specified
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = 12, height = 10, dpi = PLOT_DPI, units = "in")
    cat("Drug-target network plot saved to:", output_file, "\n")
  }
  
  return(p)
}

cat("Drug analysis module loaded successfully.\n")
