# =============================================================================
# MAP-D R Package Requirements
# Metabolic Atlas of Progression to Diabetes (MAP-D)
# =============================================================================

# Required R version
if (getRversion() < "4.0.0") {
  stop("R version 4.0.0 or higher is required")
}

# Core packages for data manipulation and analysis
core_packages <- c(
  "tidyverse",      # >= 1.3.0 - Data manipulation and visualization
  "dplyr",          # >= 1.0.0 - Data manipulation
  "tidyr",          # >= 1.1.0 - Data tidying
  "readr",          # >= 1.4.0 - Data reading
  "stringr",        # >= 1.4.0 - String manipulation
  "purrr"           # >= 0.3.0 - Functional programming
)

# Statistical analysis packages
stats_packages <- c(
  "glmnet",         # >= 4.1.0 - LASSO and elastic net
  "broom",          # >= 0.7.0 - Tidy statistical output
  "tableone",       # >= 0.13.0 - Table 1 creation
  "parallel"        # Base R - Parallel processing
)

# Data I/O packages
io_packages <- c(
  "fst",            # >= 0.9.4 - Fast data storage
  "readxl",         # >= 1.3.0 - Excel file reading
  "xml2"            # >= 1.3.0 - XML parsing for DrugBank
)

# Visualization packages
viz_packages <- c(
  "ggplot2",        # >= 3.3.0 - Grammar of graphics
  "ggrepel",        # >= 0.9.0 - Text repelling for plots
  "ggraph",         # >= 2.0.0 - Network visualization
  "igraph",         # >= 1.2.0 - Graph analysis
  "scales"          # >= 1.1.0 - Scale functions for visualization
)

# Bioconductor packages (optional, for enhanced analysis)
bioc_packages <- c(
  "limma",          # Linear models for microarray/proteomics data
  "qvalue"          # Q-value estimation for multiple testing
)

# Development and utility packages
dev_packages <- c(
  "devtools",       # Development tools
  "roxygen2",       # Documentation
  "testthat",       # Unit testing
  "knitr",          # Dynamic report generation
  "rmarkdown"       # R Markdown documents
)

# Combine all required packages
required_packages <- c(
  core_packages,
  stats_packages,
  io_packages,
  viz_packages
)

# Optional packages for enhanced functionality
optional_packages <- c(
  bioc_packages,
  dev_packages
)

#' Install required packages
#' 
#' @param packages Character vector of package names
#' @param install_optional Logical, whether to install optional packages
#' @return NULL (packages are installed)
install_mapd_packages <- function(packages = required_packages, 
                                 install_optional = FALSE) {
  
  cat("Installing MAP-D required packages...\n")
  
  # Check which packages are missing
  installed_packages <- rownames(installed.packages())
  missing_packages <- setdiff(packages, installed_packages)
  
  if (length(missing_packages) == 0) {
    cat("All required packages are already installed.\n")
  } else {
    cat("Installing", length(missing_packages), "missing packages:\n")
    cat(paste(missing_packages, collapse = ", "), "\n")
    
    # Install missing packages
    install.packages(missing_packages, dependencies = TRUE)
  }
  
  # Install optional packages if requested
  if (install_optional) {
    cat("\nInstalling optional packages...\n")
    
    # Install Bioconductor packages
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    
    missing_bioc <- setdiff(bioc_packages, installed_packages)
    if (length(missing_bioc) > 0) {
      BiocManager::install(missing_bioc)
    }
    
    # Install development packages
    missing_dev <- setdiff(dev_packages, installed_packages)
    if (length(missing_dev) > 0) {
      install.packages(missing_dev)
    }
  }
  
  cat("Package installation completed.\n")
}

#' Check package versions and compatibility
#' 
#' @param packages Character vector of package names to check
#' @return Data frame with package versions
check_package_versions <- function(packages = required_packages) {
  
  cat("Checking package versions...\n")
  
  version_info <- data.frame(
    package = character(),
    installed_version = character(),
    is_installed = logical(),
    stringsAsFactors = FALSE
  )
  
  for (pkg in packages) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      version <- as.character(packageVersion(pkg))
      version_info <- rbind(version_info, data.frame(
        package = pkg,
        installed_version = version,
        is_installed = TRUE,
        stringsAsFactors = FALSE
      ))
    } else {
      version_info <- rbind(version_info, data.frame(
        package = pkg,
        installed_version = "Not installed",
        is_installed = FALSE,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Print summary
  cat("\nPackage Status Summary:\n")
  cat("Installed:", sum(version_info$is_installed), "of", nrow(version_info), "packages\n")
  
  if (any(!version_info$is_installed)) {
    cat("\nMissing packages:\n")
    missing <- version_info[!version_info$is_installed, ]
    for (i in seq_len(nrow(missing))) {
      cat("  -", missing$package[i], "\n")
    }
  }
  
  return(version_info)
}

#' Generate session information for reproducibility
#' 
#' @param output_file File to save session info (optional)
#' @return Session info object
generate_session_info <- function(output_file = NULL) {
  
  session_info <- sessionInfo()
  
  if (!is.null(output_file)) {
    sink(output_file)
    cat("MAP-D Analysis Session Information\n")
    cat("Generated on:", as.character(Sys.time()), "\n")
    cat("=" , rep("=", 50), "=\n", sep = "")
    print(session_info)
    sink()
    cat("Session information saved to:", output_file, "\n")
  }
  
  return(session_info)
}

# Print package requirements if script is run directly
if (!interactive()) {
  cat("MAP-D Package Requirements\n")
  cat("=" , rep("=", 30), "=\n", sep = "")
  
  cat("\nRequired packages:\n")
  cat(paste(required_packages, collapse = ", "), "\n")
  
  cat("\nOptional packages:\n")
  cat(paste(optional_packages, collapse = ", "), "\n")
  
  cat("\nTo install all required packages, run:\n")
  cat("source('requirements.R')\n")
  cat("install_mapd_packages()\n")
  
  cat("\nTo check current package status, run:\n")
  cat("check_package_versions()\n")
}
