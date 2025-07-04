# Copyright (c) 2024 Omid Arhami omid.arhami@uga.edu
# R/transformations.R

#' Data Format Transformations for Antigenic Cartography
#' 
#' @description
#' This file contains functions for transforming data between different formats 
#' used in antigenic cartography. Functions handle conversion between:
#' - Long and matrix formats
#' - Distance and titer measurements
#' - Handling of threshold measurements (< and >)
#' 
#' @importFrom data.table data.table setDT fread
#' @keywords internal
"_PACKAGE"


#' Convert Long Format Data to Distance Matrix
#'
#' @description
#' Converts a dataset from long format to a symmetric distance matrix. The function
#' handles antigenic cartography data where measurements may exist between antigens
#' and antisera points. Row and column names can be optionally sorted by a time 
#' variable.
#'
#' @param data Data frame in long format
#' @param chnames Character. Name of column holding the challenge point names.
#' @param chorder Character. Optional name of column for challenge point ordering.
#' @param rnames Character. Name of column holding reference point names.
#' @param rorder Character. Optional name of column for reference point ordering.
#' @param values_column Character. Name of column containing distance/difference values. It should be from the nature of "distance" (e.g., antigenic distance or IC50), not "similarity" (e.g., HI Titer.)
#' @param rc Logical. If TRUE, reference points are treated as a subset of challenge
#'        points. If FALSE, they are treated as distinct sets. Default is FALSE.
#' @param sort Logical. Whether to sort rows/columns by chorder/rorder. Default FALSE.
#' 
#' @details
#' The function expects data in long format with at least three columns:
#' - A column for challenge point names
#' - A column for reference point names  
#' - A column containing the distance/difference values
#' 
#' Optionally, ordering columns can be provided to sort the output matrix.
#' The 'rc' parameter determines how to handle shared names between references
#' and challenges.
#'
#' @return A symmetric `matrix` of distances with row and column names corresponding 
#'         to the unique points in the data. `NA` values represent unmeasured pairs.
#'
#' @examples
#' data <- data.frame(
#'   antigen = c("A", "B", "A"),
#'   serum = c("X", "X", "Y"), 
#'   distance = c(2.5, 1.8, 3.0),
#'   year = c(2000, 2001, 2000)
#' )
#' 
#' # Basic conversion
#' mat <- long_to_matrix(data, 
#'                      chnames = "antigen",
#'                      rnames = "serum",
#'                      values_column = "distance")
#'                      
#' # With sorting by year
#' mat_sorted <- long_to_matrix(data,
#'                             chnames = "antigen",
#'                             chorder = "year",
#'                             rnames = "serum", 
#'                             rorder = "year",
#'                             values_column = "distance",
#'                             sort = TRUE)
#' @importFrom data.table setDT
#' @importFrom stats na.omit
#' @export
long_to_matrix <- function(data, chnames, chorder = NULL, 
                          rnames, rorder = NULL, values_column, 
                          rc = FALSE, sort = FALSE) {
  
  # Validate inputs 
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  
  required_cols <- c(chnames, rnames, values_column)
  if (!all(required_cols %in% names(data))) {
    missing <- setdiff(required_cols, names(data))
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }
  
  # Validate order columns if specified
  if (!is.null(chorder) && !(chorder %in% names(data))) {
    stop("chorder column '", chorder, "' not found in data")
  }
  
  if (!is.null(rorder) && !(rorder %in% names(data))) {
    stop("rorder column '", rorder, "' not found in data") 
  }
  
  # Validate numeric/character columns
  if (!is.character(data[[chnames]]) && !is.factor(data[[chnames]])) {
    stop("Challenge names column must be character or factor")
  }
  
  if (!is.character(data[[rnames]]) && !is.factor(data[[rnames]])) {
    stop("Reference names column must be character or factor")
  }
  
  if (!is.numeric(data[[values_column]]) && 
      !all(grepl("^[0-9<>]", na.omit(data[[values_column]])))) {
    stop("Values column must be numeric or contain valid threshold indicators (< or >)")
  }
  
  # Convert to data.table for efficiency
  data.table::setDT(data)
  
  if (rc == FALSE) {
    # Mark antigens and antisera
    data[, (chnames) := paste0("V/", get(chnames))]
    data[, (rnames) := paste0("S/", get(rnames))]
  }
  
  # Get unique point names
  #all_points <- sort(unique(unlist(data[, .(get(chnames), get(rnames))])))
  all_points <- sort(unique(c(data[[chnames]], data[[rnames]])))

  # Create square matrix with NA values
  n <- length(all_points)
  distance_matrix <- matrix(NA, nrow = n, ncol = n)
  rownames(distance_matrix) <- all_points
  colnames(distance_matrix) <- all_points

  if (sort == TRUE) {
    # Get one rank per name
    ranks <- numeric(length(all_points))
    for (i in seq_along(all_points)) {
      name <- all_points[i]
      yr <- 0
      
      # Try to get rank from challenge order
      if (!is.null(chorder)) {
        name_rank <- unique(data[get(chnames) == name, get(chorder)])
        if (length(name_rank) > 0) yr <- min(name_rank)
      }
      
      # If not found, try reference order
      if (yr == 0 && !is.null(rorder)) {
        name_rank <- unique(data[get(rnames) == name, get(rorder)])
        if (length(name_rank) > 0) yr <- min(name_rank)
      }
      
      ranks[i] <- yr
    }
    
    ranks <- as.numeric(ranks)
    
    # Reorder matrix by ranks
    idx <- order(ranks)
    distance_matrix <- distance_matrix[idx, idx]
  }

  # Fill in the distances
  for (i in seq_len(nrow(data))) {
    r   <- data[i, get(chnames)]
    c   <- data[i, get(rnames)]
    val <- data[i, get(values_column)]

    # Set both matrix elements for symmetry
    distance_matrix[r, c] <- val
    distance_matrix[c, r] <- val
  }
  
  # Set diagonal to 0
  diag(distance_matrix) <- 0
  
  return(distance_matrix)
}



#' Convert Distance Matrix to Titer Panel Format
#'
#' @description
#' Converts a distance matrix to a titer panel format, handling threshold measurements
#' and logarithmic transformations common in antigenic cartography. The function 
#' identifies reference points (typically antisera) and challenge points (typically
#' antigens) based on row/column name prefixes.
#'
#' @param input_matrix Matrix of distances, with row/column names prefixed with 
#'        "V/" for antigens and "S/" for sera
#' @param base Numeric. Base for logarithmic transformation. Default exp(1). For HI Assay 2
#' @param tens Numeric. Scaling factor for final titers. Default 1. For HI Assay 10
#'
#' @details 
#' The function:
#' 1. Identifies antigen and serum entries from matrix row/column names
#' 2. Creates titer table from antigen-serum pairs
#' 3. Handles threshold indicators (< and >) in distance values
#' 4. Applies appropriate transformations to convert distances to titers
#'
#' Transformation steps:
#' 1. Extract numeric values from thresholded measurements
#' 2. Convert distances to titers via logarithmic transformation
#' 3. Apply scaling factor
#' 4. Reapply threshold indicators to transformed values
#'
#' @return A `matrix` of titers with:
#' - Rows corresponding to antigen strains (without "V/" prefix)
#' - Columns corresponding to antisera (without "S/" prefix)
#' - Values as character strings including threshold indicators where applicable
#' - `NA` values replaced with "*"
#'
#' @examples
#' # Create sample distance matrix
#' dist_mat <- matrix(c(0, 2, ">3", 2, 0, 4, "3", 4, 0), nrow=3)
#' rownames(dist_mat) <- c("V/strain1", "V/strain2", "S/serum1")
#' colnames(dist_mat) <- c("V/strain1", "V/strain2", "S/serum1")
#'
#' # Convert to titer panel
#' titer_panel <- dist_to_titer_table(dist_mat)
#'
#' @export
dist_to_titer_table <- function(input_matrix, base=exp(1), tens=1) {
  # Input validation
  if (!is.matrix(input_matrix)) {
    stop("input_matrix must be a matrix")
  }
  
  if (!all(apply(input_matrix, 2, function(x) is.numeric(x) || is.character(x)))) {
    stop("Matrix values must be numeric or character")
  }
  
  if (is.null(rownames(input_matrix)) || is.null(colnames(input_matrix))) {
    stop("Matrix must have row and column names")
  }

  # Validate prefixes
  if(!any(startsWith(rownames(input_matrix), "S/")) ||
     !any(startsWith(rownames(input_matrix), "V/"))) {
    stop("Matrix must have both antiserum (S/) and virus (V/) entries")
  }
  
  if (!is.numeric(base) || base <= 0) {
    stop("base must be a positive number")
  }
  
  if (!is.numeric(tens) || tens <= 0) {
    stop("tens must be a positive number")
  }

  # Creating the panel of data:
  # Find row names that start with "S/"
  serum_names_S <- rownames(input_matrix)[startsWith(rownames(input_matrix), "S/")]
  
  # Remove "S/" from the beginning of each element in serum_names
  serum_names <- sub("^S/", "", serum_names_S)
  
  # Find row names that start with "V/"
  antigen_names_V <- rownames(input_matrix)[startsWith(rownames(input_matrix), "V/")]
  
  # Remove "V/" from the beginning of each element 
  antigen_names <- sub("^V/", "", antigen_names_V)
  
  # # Create table for antigen-serum pairs
  # distance_table <- matrix(NA, 
  #                        nrow = length(antigen_names_V), 
  #                        ncol = length(serum_names_S),
  #                        dimnames = list(antigen_names_V, serum_names_S))
  
  # # Fill values from input matrix
  # for (antigen in antigen_names_V) {
  #   for (serum in serum_names_S) {
  #     if (antigen %in% rownames(input_matrix) && serum %in% colnames(input_matrix)) {
  #       distance_table[antigen, serum] <- input_matrix[antigen, serum]
  #     }
  #   }
  # }

  # Convert the distances to titers in a panel format
  distance_table = symmetric_to_nonsymmetric_matrix(input_matrix, serum_names_S)
  
  
  # Set row and column names
  rownames(distance_table) <- antigen_names
  colnames(distance_table) <- serum_names

  # Handle threshold values
  replace_greater_than <- function(x) {
    if (is.character(x) && grepl("^>", x)) {
      as.numeric(sub("^>", "", x))
    } else {
      x
    }
  }
  
  replace_smaller_than <- function(x) {
    if (is.character(x) && grepl("^<", x)) {
      as.numeric(sub("^<", "", x))
    } else {
      x
    }
  }
  
  # Process matrix values
  numeric_dist_table <- apply(distance_table, c(1, 2), replace_greater_than)
  numeric_dist_table <- apply(numeric_dist_table, c(1, 2), replace_smaller_than)
  numeric_dist_table <- apply(numeric_dist_table, 2, as.numeric)
  
  # Find column maxima
  column_maxes <- apply(numeric_dist_table, 2, max, na.rm = TRUE)

  # Create an empty matrix to hold the results, preserving dimensions
  titer_table <- matrix(NA, 
                        nrow = nrow(distance_table), 
                        ncol = ncol(distance_table),
                        dimnames = list(rownames(distance_table), colnames(distance_table)))

  # Loop through the matrix to calculate titers
  for (j in 1:ncol(distance_table)) {
    for (i in 1:nrow(distance_table)) {
      value <- distance_table[i, j]
      max_value <- column_maxes[j]

      if (is.na(value)) {
        titer_table[i, j] <- NA
      } else if (is.character(value) && startsWith(value, ">")) {
        num_part <- as.numeric(sub(">", "", value))
        titer_table[i, j] <- paste0("<", max_value - num_part)
      } else if (is.character(value) && startsWith(value, "<")) {
        num_part <- as.numeric(sub("<", "", value))
        titer_table[i, j] <- paste0(">", max_value - num_part)
      } else {
        titer_table[i, j] <- as.character(max_value - as.numeric(value))
      }
    }
  }

  # Final transformation to titer scale
  transform_value <- function(x, base=base, tens=tens) {
    if (is.na(x)) {
      return(NA)
    } else if (is.character(x) && grepl("^<", x)) {
      num_part <- as.numeric(sub("^<", "", x))
      return(paste0("<", (base^num_part) * tens))
    } else if (is.character(x) && grepl("^>", x)) {
      num_part <- as.numeric(sub("^>", "", x))
      return(paste0(">", (base^num_part) * tens))
    } else {
      return((base^as.numeric(x)) * tens)
    }
  }
  
  # Apply transformation
  titer_table <- apply(titer_table, c(1, 2), 
                       function(x) transform_value(x, base=base, tens=tens))

  # Replace NA with *
  titer_table[is.na(titer_table)] <- '*'
  
  return(titer_table)
}


#' Convert distance matrix to assay panel format
#'
#' @param dist_matrix Distance matrix
#' @param selected_names Names of reference points
#' @return A non-symmetric `matrix` in assay panel format, where rows are test antigens and columns are reference antigens.
#' @export
symmetric_to_nonsymmetric_matrix <- function(dist_matrix, selected_names) {
  if (!is.matrix(dist_matrix)) {
    stop("dist_matrix must be a matrix")
  }
  if (!is.character(selected_names)) {
    stop("selected_names must be a character vector") 
  }
  
  # Subset matrix keeping only virus rows and selected sera columns
  panel <- dist_matrix[!rownames(dist_matrix) %in% selected_names, selected_names, drop = FALSE]
  return(panel)
}


#' Convert coordinates to distance matrix
#'
#' Calculates pairwise Euclidean distances between points in coordinate space
#' 
#' @param positions Matrix of coordinates where rows are points and columns are dimensions;
#' It can be a matrix or a data frame.
#' @return A symmetric `matrix` of pairwise Euclidean distances between points.
#' @importFrom stats dist
#' @export
coordinates_to_matrix <- function(positions) {
  if (!is.matrix(positions) && !is.data.frame(positions)) {
    stop("positions must be a matrix or a data frame")
  }
  if (is.data.frame(positions)) {
    positions <- as.matrix(positions)
  }
  
  num_points <- nrow(positions)
  p_dist_mat <- as.matrix(stats::dist(positions))
  
  # Set row and column names if they exist
  if (!is.null(rownames(positions))) {
    rownames(p_dist_mat) <- rownames(positions)
    colnames(p_dist_mat) <- rownames(positions)
  }
  
  return(p_dist_mat)
}


#' Filter matrix to only virus vs antiserum distances
#'
#' @param dist_matrix Distance matrix
#' @param selected_names Names of selected reference points
#' @return A `matrix` of the same dimensions as the input, but with non-virus-vs-antiserum distances set to `NA`.
#' @export
only_virus_vs_as <- function(dist_matrix, selected_names) {
  if (!is.matrix(dist_matrix)) {
    stop("dist_matrix must be a matrix")
  }
  if (!is.character(selected_names)) {
    stop("selected_names must be a character vector")
  }
  
  # Make values NA for non-selected columns
  dist_matrix[, !colnames(dist_matrix) %in% selected_names] <- NA
  
  # Update symmetric values
  dist_matrix[rownames(dist_matrix) %in% selected_names, ] <- 
    t(dist_matrix[, colnames(dist_matrix) %in% selected_names])
  
  return(dist_matrix)
}


#' Log Transform Parameter Samples
#'
#' @description
#' Reads samples from a CSV file and log transforms specific parameters (N, k0, cooling_rate, c_repulsion)
#' if they exist in the data. If `output_file` is specified, the transformed data is saved. Otherwise,
#' the transformed data frame is returned.
#'
#' @param samples_file Character. Path to CSV file containing samples.
#' @param output_file Character. Optional path (including filename) for saving transformed data as a CSV. 
#' If NULL (the default), the function returns the transformed data frame without writing a file.
#' @return A `data.frame` with log-transformed parameters. If `output_file` is specified, the function also
#'         writes the data frame to the specified path and returns it invisibly.
#' @examples
#' # This example uses a sample file included with the package.
#' sample_file <- system.file("extdata", "sample_params.csv", package = "topolow")
#' 
#' # Ensure the file exists before running the example
#' if (nzchar(sample_file)) {
#'   # Transform the data from the sample file and return as a data frame
#'   transformed_data <- log_transform_parameters(sample_file, output_file = NULL)
#'   
#'   # Display the first few rows of the transformed data
#'   print(head(transformed_data))
#' }
#' 
#' @importFrom utils read.csv write.csv
#' @importFrom stats na.omit
#' @export
log_transform_parameters <- function(samples_file, output_file = NULL) {
  # Validate input file
  if (!file.exists(samples_file)) {
    stop("Input file not found: ", samples_file)
  }
  
  # Read samples 
  samples <- tryCatch({
    utils::read.csv(samples_file, stringsAsFactors = FALSE)
  }, error = function(e) {
    stop("Error reading file: ", e$message)
  })
  
  # Remove rows with NA values
  samples <- stats::na.omit(samples)
  
  # Check which parameters exist
  params_to_transform <- c("N", "k0", "cooling_rate", "c_repulsion")
  existing_params <- intersect(names(samples), params_to_transform)
  
  if (length(existing_params) == 0) {
    message("No parameters found to transform or all parameters are already transformed.")
    return(samples)
  }
  
  # Validate numeric columns
  for (param in existing_params) {
    if (!is.numeric(samples[[param]])) {
      samples[[param]] <- suppressWarnings(as.numeric(samples[[param]]))
      if (any(is.na(samples[[param]]))) {
        stop("Non-numeric values found in column: ", param)
      }
    }
    # Check for non-positive values
    if (any(samples[[param]] <= 0, na.rm = TRUE)) {
      stop("Non-positive values found in column: ", param,
           ". Log transform requires positive values.")
    }
  }
  
  # Create log-transformed columns
  for (param in existing_params) {
    log_param <- paste0("log_", param)
    samples[[log_param]] <- log(samples[[param]])
  }
  
  # Remove original columns
  samples <- samples[, !names(samples) %in% existing_params, drop = FALSE]
  
  # Reorder to fixed column order
  desired_order <- c(
    "log_N",
    "log_k0",
    "log_cooling_rate",
    "log_c_repulsion",
    "Holdout_MAE",
    "NLL"
  )
  missing_cols <- setdiff(desired_order, names(samples))
  if (length(missing_cols) > 0) {
    stop("Missing expected columns: ", paste(missing_cols, collapse = ", "))
  }
  samples <- samples[, desired_order, drop = FALSE]
  
  # Write output if a path is provided
  if (!is.null(output_file)) {
      if (!is.character(output_file) || length(output_file) != 1) {
          stop("'output_file' must be a single character string path.", call. = FALSE)
      }
  tryCatch({
    utils::write.csv(samples, file = output_file, row.names = FALSE)
  }, error = function(e) {
    stop("Error writing output file: ", e$message)
  })
      message("Log transformed parameters: ", paste(existing_params, collapse = ", "))
  message("Output saved to: ", output_file)
      return(invisible(samples))
  }
  
  # Return the transformed data frame if not writing to file
  return(samples)
}
