# Copyright (c) 2024 Omid Arhami o.arhami@gmail.com
# Licensed under Pre-Publication Academic License https://github.com/omid-arhami/topolow/blob/main/LICENSE

# Script to check and install required dependencies
# To be used in SLURM environment or other systems where dependencies might be missing

# List of required packages
required_packages <- c(
  "ggplot2", "dplyr", "data.table", "reshape2", "plotly", "Racmacs", 
  "parallel", "coda", "MASS", "vegan", "igraph", "lhs", "umap",
  "gridExtra", "scales"
)

# Check which packages are missing
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

# Install missing packages if needed
if (length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

# Verify all packages are now available
still_missing <- missing_packages[!sapply(missing_packages, requireNamespace, quietly = TRUE)]
if (length(still_missing) > 0) {
  stop("Failed to install the following packages: ", paste(still_missing, collapse = ", "))
} else {
  cat("All required packages are available.\n")
}
