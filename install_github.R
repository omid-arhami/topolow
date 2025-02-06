# install_github.R
# Installation script for topolow package
# Copyright (c) 2024 Omid Arhami o.arhami@gmail.com
# Licensed under MIT License

# Install devtools if needed
if (!require("devtools")) install.packages("devtools")

# Install required dependencies
required_packages <- c(
  "ggplot2", "dplyr", "data.table", "reshape2", "plotly", "rgl", 
  "Racmacs", "parallel", "coda", "MASS", "vegan", "igraph", 
  "lhs", "umap", "gridExtra", "scales", "colorspace"
)

for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
  }
}

# Install topolow from GitHub
devtools::install_github("omid-arhami/topolow")

# Test installation
library(topolow)
?topolow  # Should open package documentation
