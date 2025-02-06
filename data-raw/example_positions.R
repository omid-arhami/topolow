# data-raw/example_positions.R
# Script to prepare example_positions dataset

# Load raw data
raw_positions <- read.csv("data-raw/example_positions.csv")

# Verify required columns
required_cols <- c(paste0("V", 1:5), "name", "antigen", "antiserum", 
                  "cluster", "color", "year")
if(!all(required_cols %in% names(raw_positions))) {
  stop("Missing required columns in raw data")
}

# Convert logical columns
raw_positions$antigen <- as.logical(raw_positions$antigen)
raw_positions$antiserum <- as.logical(raw_positions$antiserum)

# Convert cluster to factor
raw_positions$cluster <- factor(raw_positions$cluster)

# Verify coordinates are numeric
coord_cols <- paste0("V", 1:5)
raw_positions[coord_cols] <- lapply(raw_positions[coord_cols], as.numeric)

# Remove any rows with NAs
example_positions <- na.omit(raw_positions)

# Verify no remaining NAs
if(any(is.na(example_positions))) {
  stop("NA values remain in processed data")
}

# Compress to save space
example_positions <- as.data.frame(example_positions)

# Save to data/ directory
usethis::use_data(example_positions, overwrite = TRUE, compress = "xz")