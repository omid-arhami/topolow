# Copyright (c) 2024 Omid Arhami o.arhami@gmail.com
# Licensed under Pre-Publication Academic License https://github.com/omid-arhami/topolow/blob/main/LICENSE

# install_github.R
# Installation script for topolow package

# Install devtools if needed
if (!require("devtools")) install.packages("devtools")

# Install required dependencies
required_packages <- c(
  "ggplot2", "dplyr", "data.table", "reshape2", "plotly", "rgl", 
  "Racmacs", "parallel", "coda", "MASS", "vegan", "igraph", 
  "lhs", "umap", "gridExtra", "scales", "filelock"
)

for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
  }
}

# Check if all packages are installed
missing_pkgs <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Failed to install the following packages: ", paste(missing_pkgs, collapse = ", "))
}

# Install topolow from GitHub
devtools::install_github("omid-arhami/topolow")

# Test installation
library(topolow)
?topolow  # Should open package documentation

# Print success message
cat("topolow has been successfully installed with all required dependencies.\n")