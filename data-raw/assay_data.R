# data-raw/assay_data.R
# Script to prepare HIV and H3N2 datasets

####################
# Load H3N2 data. REF: 10.1126/science.1097211
h3n2_raw <- read.csv("data-raw/Smith2004-data.csv")

# Basic validation
required_cols <- c("virusStrain", "serumStrain", "titer", 
                  "virusYear", "serumYear", "cluster", "color")
if(!all(required_cols %in% names(h3n2_raw))) {
  stop("Missing required columns in H3N2 data")
}

# Clean and compress
h3n2_data <- h3n2_raw
h3n2_data$cluster <- factor(h3n2_data$cluster)


####################
# Load DENV data. Ref: https://doi.org/10.7554/eLife.42496
denv_raw <- read.csv("data-raw/DENV_titers.csv")

# Basic validation
required_cols <- c("virus_strain", "serum_strain", "titer", 
                  "virusYear", "serumYear", "cluster", "color")
if(!all(required_cols %in% names(denv_raw))) {
  stop("Missing required columns in DENV data")
}

# Clean and compress
denv_data <- denv_raw
denv_data$cluster <- factor(denv_data$cluster)


####################
# Load HIV data
hiv_viruses <- read.csv("data-raw/HIV_viruses.csv")
hiv_titers <- read.csv("data-raw/hiv_assay_titers.csv")

# Validate HIV data
if(!all(c("Virus.name", "Country", "Subtype", "Year") %in% names(hiv_viruses))) {
  stop("Missing required columns in HIV viruses data")
}

if(!all(c("Antibody", "Virus", "IC50") %in% names(hiv_titers))) {
  stop("Missing required columns in HIV titers data")
}

# Save processed datasets
usethis::use_data(h3n2_data, overwrite = TRUE, compress = "xz")
usethis::use_data(denv_data, overwrite = TRUE, compress = "xz")
usethis::use_data(hiv_viruses, overwrite = TRUE, compress = "xz")
usethis::use_data(hiv_titers, overwrite = TRUE, compress = "xz")
