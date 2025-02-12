# build_package.R
# Script for building, testing, and preparing the topolow package
# Copyright (c) 2024 Omid Arhami omid.arhami@uga.edu
# Licensed under MIT License

#' Setup and Installation
#' This script handles package development tasks including:
#' - Installing dependencies
#' - Building documentation
#' - Running tests
#' - Creating distributions
#' - Preparing for GitHub

# Install required development packages if needed
if (!require("devtools")) install.packages("devtools")
if (!require("roxygen2")) install.packages("roxygen2")
if (!require("testthat")) install.packages("testthat")
if (!require("covr")) install.packages("covr")
if (!require("knitr")) install.packages("knitr")
if (!require("rmarkdown")) install.packages("rmarkdown")

# Load devtools
library(devtools)

# Clean any previous builds
unlink("NAMESPACE")
devtools::clean_dll()
devtools::clean_vignettes()

# Process data files
source("data-raw/example_positions.R")
source("data-raw/assay_data.R")

# Update documentation
devtools::document()

# Run tests
devtools::test()

# Check test coverage
covr::package_coverage()

# Run comprehensive package checks
options(warn = 1)  # Show warnings immediately
Sys.setenv(R_COMPILE_AND_INSTALL_PACKAGES = "true")
devtools::check(document = TRUE, 
                manual = TRUE, 
                force_suggests = TRUE, 
                cran = TRUE)

# Build package distributions
devtools::build()                   # Source package (.tar.gz)
devtools::build(binary = TRUE)      # Windows binary (.zip)

# Install locally
devtools::install(build_vignettes = TRUE)

# GitHub preparation steps
usethis::use_github_action_check_standard()  # Add GitHub Actions CI workflow
usethis::use_github_action("test-coverage")  # Add coverage workflow