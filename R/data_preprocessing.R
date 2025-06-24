# Copyright (c) 2024 Omid Arhami omid.arhami@uga.edu
# License: free of charge access granted to any academic researcher to use this software for non-commercial, academic research purposes **only**.  Nobody may modify, distribute, sublicense, or publicly share the Software or any derivative works, until the paper is published by the original authors.  The Software is provided "as is" without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose and noninfringement.  In no event shall the authors or copyright holders be liable for any claim, damages or other liability, whether in an action of contract, tort or otherwise, arising from, out of or in connection with the Software or the use or other dealings in the Software.

# R/data_preprocessing.R

#' Antigenic Data Preprocessing Functions
#' 
#' @description
#' Functions for standardizing and preprocessing antigenic assay data from various 
#' sources into consistent formats. Handles titer and IC50 measurements, threshold
#' values, and produces both long and matrix formats suitable for mapping.
#'
#' @keywords internal
"_PACKAGE"

#' Process Raw Antigenic Assay Data
#'
#' @description
#' Processes raw antigenic assay data from CSV files into standardized long and matrix
#' formats. Handles both titer data (which needs conversion to distances) and direct
#' distance measurements like IC50. Preserves threshold indicators (<, >) and handles
#' repeated measurements by averaging.
#'
#' @param file_path Character. Path to CSV file containing raw data.
#' @param antigen_col Character. Name of column containing virus/antigen identifiers.
#' @param serum_col Character. Name of column containing serum/antibody identifiers.
#' @param value_col Character. Name of column containing measurements (titers or distances).
#' @param is_titer Logical. Whether values are titers (TRUE) or distances like IC50 (FALSE).
#' @param metadata_cols Character vector. Names of additional columns to preserve.
#' @param id_prefix Logical. Whether to prefix IDs with V/ and S/ (default: TRUE).
#' @param base Numeric. Base for logarithm transformation (default: 2 for titers, e for IC50).
#' @param scale_factor Numeric. Scale factor for titers (default: 10).
#'
#' @return List containing:
#'   \item{long}{Data frame in long format with standardized columns}
#'   \item{matrix}{Distance matrix}
#'
#' @details
#' The function handles these key steps:
#' 1. Reads and validates input data
#' 2. Transforms values to log scale
#' 3. Converts titers to distances if needed
#' 4. Averages repeated measurements
#' 5. Creates standardized long format
#' 6. Creates distance matrix
#' 7. Preserves metadata and threshold indicators
#' 8. Preserves virusYear and serumYear columns if present
#' 
#' Input requirements and constraints:
#' * CSV file must contain required columns
#' * Column names must match specified parameters in the function input
#' * Values can include threshold indicators (< or >)
#' * Metadata columns must exist if specified
#' * Allowed Year-related column names are "virusYear" and "serumYear"
#'
#' @examples
#' \dontrun{
#' # Process titer data (e.g., HI assay)
#' results <- process_antigenic_data(
#'   "smith2004.csv",
#'   antigen_col = "virusStrain",
#'   serum_col = "serumStrain", 
#'   value_col = "titer",
#'   is_titer = TRUE,
#'   metadata_cols = c("cluster", "color")
#' )
#'
#' # Process IC50 data
#' results <- process_antigenic_data(
#'   "hiv_assays.csv",
#'   antigen_col = "Virus",
#'   serum_col = "Antibody",
#'   value_col = "IC50",
#'   is_titer = FALSE
#' )
#' }
#' @importFrom utils read.csv
#' @importFrom dplyr %>% group_by mutate ungroup summarise select distinct left_join
#' @importFrom rlang sym
#' @importFrom stats na.omit
#' @export
process_antigenic_data <- function(file_path, antigen_col, serum_col, 
                                   value_col,
                                   is_titer = TRUE, 
                                   metadata_cols = NULL,
                                   id_prefix = FALSE,
                                   base = NULL, 
                                   scale_factor = 10) {
  # Input validation
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  
  # Read data
  data <- utils::read.csv(file_path)
  
  # Validate required columns
  req_cols <- c(antigen_col, serum_col, value_col)
  if (!all(req_cols %in% names(data))) {
    missing <- setdiff(req_cols, names(data))
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }
  
  # Clean invalid values
  data <- data[!is.na(data[[value_col]]), ]  # Remove NA values
  data <- data[data[[value_col]] != "", ]    # Remove empty strings
  
  # Keep only rows where value starts with a number or < or >
  data <- data[grepl("^[0-9<>]", data[[value_col]]), ]
  
  if (nrow(data) == 0) {
    stop("No valid measurements remaining after cleaning")
  }
  
  # Check for year columns
  year_cols <- intersect(c("virusYear", "serumYear"), names(data))
  
  # Add year columns to metadata if they exist
  if (length(year_cols) > 0) {
    metadata_cols <- unique(c(metadata_cols, year_cols))
  }
  
  # Validate metadata columns if specified
  if (!is.null(metadata_cols)) {
    missing_meta <- setdiff(metadata_cols, names(data))
    if (length(missing_meta) > 0) {
      stop("Missing metadata columns: ", paste(missing_meta, collapse = ", "))
    }
  }
  
  
  if (!is.logical(is_titer)) {
    stop("is_titer must be logical")
  }
  
  if (!is.null(base) && (!is.numeric(base) || base <= 0)) {
    stop("base must be NULL or a positive number")
  }
  
  if (!is.numeric(scale_factor) || scale_factor <= 0) {
    stop("scale_factor must be a positive number")
  }

  # Set default base if not provided
  if (is.null(base)) {
    base <- if(is_titer) 2 else exp(1)
  }

  # Helper function definitions from your code
  remove_sign <- function(x) {
    as.numeric(gsub("[<>]", "", x))
  }

  reapply_sign <- function(values, avg) {
    if (any(grepl("[<>]", values))) {
      sign <- ifelse(any(grepl("<", values)), "<", ">")
      return(paste0(sign, avg))
    } else {
      return(as.character(avg))
    }
  }

  # Process values
  if (is_titer) {
    # Convert titers to log scale
    data$log_value <- sapply(data[[value_col]], function(x) {
      if (grepl("^<", x)) {
        paste0("<", log(as.numeric(sub("<", "", x)) / scale_factor, base = base))
      } else if(grepl("^>", x)){
        paste0(">", log(as.numeric(sub(">", "", x)) / scale_factor, base = base))
      } else if (is.numeric(as.numeric(x))) {
        log(as.numeric(x) / scale_factor, base = base)
      } else {
        NA
      }
    })
    
    # Calculate distances using Smith's method
    data$processed_value <- sapply(data$log_value, function(x) {
      if (grepl("^<", x)) {
        as.numeric(sub("<", "", x))
      } else if(grepl("^>", x)){
        as.numeric(sub(">", "", x))
      } else {
        as.numeric(x)
      }
    })
    
    # Calculate distances per reference
    data <- data %>%
      dplyr::group_by(!!sym(serum_col)) %>%
      dplyr::mutate(
        max_value = max(processed_value, na.rm = TRUE),
        distance = max_value - processed_value
      ) %>%
      dplyr::ungroup()
    
    # Adjust distances for threshold values
    data$distance <- sapply(1:nrow(data), function(i) {
      if (grepl("^<", data$log_value[i])) {
        paste0(">", data$distance[i])
      } else if (grepl("^>", data$log_value[i])) {
        paste0("<", data$distance[i])
      } else {
        data$distance[i]
      }
    })
    
  } else {
    # For IC50, calculate log directly
    data$distance <- sapply(data[[value_col]], function(x) {
      if (grepl("^<", x)) {
        x_num <- as.numeric(sub("<", "", x))
        paste0("<", log(1 + x_num, base = base))
      } else if (grepl("^>", x)) {
        x_num <- as.numeric(sub(">", "", x))
        paste0(">", log(1 + x_num, base = base))
      } else if (is.numeric(as.numeric(x))) {
        log(1 + as.numeric(x), base = base)
      } else {
        NA
      }
    })
  }

  # Combine repeated measurements
  long_data <- data %>%
    dplyr::group_by(!!sym(antigen_col), !!sym(serum_col)) %>%
    dplyr::summarise(
      raw_value = reapply_sign(!!sym(value_col), 
                               if(is_titer) {
                                 scale_factor * base^(mean(log(remove_sign(!!sym(value_col))/scale_factor, 
                                                               base = base), na.rm = TRUE))
                               } else {
                                 exp(mean(log(remove_sign(!!sym(value_col))), na.rm = TRUE))
                               }),
      distance = reapply_sign(distance, 
                              mean(remove_sign(distance), na.rm = TRUE)),
      .groups = 'drop'
    )
  
  # Add back metadata columns if specified
  if (!is.null(metadata_cols)) {
    # Get combined distinct metadata for antigens
    metadata <- data %>%
      dplyr::select(!!sym(antigen_col), !!sym(serum_col), dplyr::any_of(metadata_cols)) %>%
      dplyr::distinct()
    
    # Join metadata based on both antigen and serum columns
    long_data <- long_data %>%
      dplyr::left_join(metadata, by = c(antigen_col, serum_col))
  }
  
  # Remove the non complete rows
  long_data <- na.omit(long_data)
  
  # sort by year to conform with our assumptions used at various places
  if("virusYear" %in% names(long_data)) {
    long_data <- long_data[order(long_data$virusYear), ]
  }

  # # Add prefixes if requested
  # if (id_prefix) {
  #   long_data[[antigen_col]] <- paste0("V/", long_data[[antigen_col]])
  #   long_data[[serum_col]] <- paste0("S/", long_data[[serum_col]])
  # }
  
  # Create matrix using long_to_matrix function
  # Determine ordering columns
  virus_year_col <- "virusYear"
  serum_year_col <- "serumYear"
  
  distance_matrix <- long_to_matrix(
    long_data,
    chnames = antigen_col,
    chorder = if("virusYear" %in% names(long_data)) "virusYear" else NULL,
    rnames = serum_col,
    rorder = if("serumYear" %in% names(long_data)) "serumYear" else NULL,
    values_column = "distance",
    rc = FALSE,
    sort = ("virusYear" %in% names(long_data)) || ("serumYear" %in% names(long_data))
  )

  # Return both formats
  return(list(
    long = long_data,
    matrix = distance_matrix
  ))
}



#' Process Raw Antigenic Assay Data without transformations
#'
#' @description
#' Processes raw antigenic assay data from CSV files into standardized long and matrix
#' formats. Handles both titer data (which needs conversion to distances) and direct
#' distance measurements like IC50. Preserves threshold indicators (<, >) and handles
#' repeated measurements by averaging.
#'
#' @param file_path Character. Path to CSV file containing raw data.
#' @param antigen_col Character. Name of column containing virus/antigen identifiers.
#' @param serum_col Character. Name of column containing serum/antibody identifiers.
#' @param value_col Character. Name of column containing measurements (titers or distances).
#' @param is_titer Logical. Whether values are titers (TRUE) or distances like IC50 (FALSE).
#' @param metadata_cols Character vector. Names of additional columns to preserve.
#' @param id_prefix Logical. Whether to prefix IDs with V/ and S/ (default: TRUE).
#' @param base Numeric. Base for logarithm transformation (default: 2 for titers, e for IC50).
#' @param scale_factor Numeric. Scale factor for titers (default: 10).
#'
#' @return List containing:
#'   \item{long}{Data frame in long format with standardized columns}
#'   \item{matrix}{Distance matrix}
#'
#' @details
#' The function handles these key steps:
#' 1. Reads and validates input data
#' 2. Transforms values to log scale
#' 3. Converts titers to distances if needed
#' 4. Averages repeated measurements
#' 5. Creates standardized long format
#' 6. Creates distance matrix
#' 7. Preserves metadata and threshold indicators
#' 8. Preserves virusYear and serumYear columns if presen
#' 
#' Input requirements and constraints:
#' * CSV file must contain required columns
#' * Column names must match specified parameters in the function input
#' * Values can include threshold indicators (< or >)
#' * Metadata columns must exist if specified
#' * Allowed Year-related column names are "virusYear" and "serumYear"
#'
#' @examples
#' \dontrun{
#' # Process titer data (e.g., HI assay)
#' results <- process_antigenic_data(
#'   "smith2004.csv",
#'   antigen_col = "virusStrain",
#'   serum_col = "serumStrain", 
#'   value_col = "titer",
#'   is_titer = TRUE,
#'   metadata_cols = c("cluster", "color")
#' )
#'
#' # Process IC50 data
#' results <- process_antigenic_data(
#'   "hiv_assays.csv",
#'   antigen_col = "Virus",
#'   serum_col = "Antibody",
#'   value_col = "IC50",
#'   is_titer = FALSE
#' )
#' }
#' @importFrom utils read.csv
#' @importFrom dplyr %>% group_by mutate ungroup summarise select distinct left_join
#' @importFrom rlang sym
#' @importFrom stats na.omit
#' @export
process_antigenic_data_notransform <- function(file_path, antigen_col, serum_col, 
                                   value_col,
                                   is_titer = TRUE, 
                                   metadata_cols = NULL,
                                   id_prefix = FALSE,
                                   base = NULL, 
                                   scale_factor = 10) {
  # Input validation
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  
  # Read data
  data <- utils::read.csv(file_path)
  
  # Validate required columns
  req_cols <- c(antigen_col, serum_col, value_col)
  if (!all(req_cols %in% names(data))) {
    missing <- setdiff(req_cols, names(data))
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }
  
  # Clean invalid values
  data <- data[!is.na(data[[value_col]]), ]  # Remove NA values
  data <- data[data[[value_col]] != "", ]    # Remove empty strings
  
  # Keep only rows where value starts with a number or < or >
  data <- data[grepl("^[0-9<>]", data[[value_col]]), ]
  
  if (nrow(data) == 0) {
    stop("No valid measurements remaining after cleaning")
  }
  
  # Check for year columns
  year_cols <- intersect(c("virusYear", "serumYear"), names(data))
  
  # Add year columns to metadata if they exist
  if (length(year_cols) > 0) {
    metadata_cols <- unique(c(metadata_cols, year_cols))
  }
  
  # Validate metadata columns if specified
  if (!is.null(metadata_cols)) {
    missing_meta <- setdiff(metadata_cols, names(data))
    if (length(missing_meta) > 0) {
      stop("Missing metadata columns: ", paste(missing_meta, collapse = ", "))
    }
  }
  
  
  if (!is.logical(is_titer)) {
    stop("is_titer must be logical")
  }
  
  if (!is.null(base) && (!is.numeric(base) || base <= 0)) {
    stop("base must be NULL or a positive number")
  }
  
  if (!is.numeric(scale_factor) || scale_factor <= 0) {
    stop("scale_factor must be a positive number")
  }
  
  # Set default base if not provided
  if (is.null(base)) {
    base <- if(is_titer) 2 else exp(1)
  }
  
  # Helper function definitions from your code
  remove_sign <- function(x) {
    as.numeric(gsub("[<>]", "", x))
  }
  
  reapply_sign <- function(values, avg) {
    if (any(grepl("[<>]", values))) {
      sign <- ifelse(any(grepl("<", values)), "<", ">")
      return(paste0(sign, avg))
    } else {
      return(as.character(avg))
    }
  }
  
  # Process values
  if (is_titer) {
    # Convert titers to log scale
    data$log_value <- sapply(data[[value_col]], function(x) {
      if (grepl("^<", x)) {
        paste0("<", log(as.numeric(sub("<", "", x)) / scale_factor, base = base))
      } else if (grepl("^>", x)) {
        paste0(">", log(as.numeric(sub(">", "", x)) / scale_factor, base = base))
      } else if (is.numeric(as.numeric(x))) {
        log(as.numeric(x) / scale_factor, base = base)
      } else {
        NA
      }
    })
    
    data$processed_value <- sapply(data$log_value, function(x) {
      if (grepl("^<", x)) {
        as.numeric(sub("<", "", x))
      } else if (grepl("^>", x)) {
        as.numeric(sub(">", "", x))
      } else {
        as.numeric(x)
      }
    })
    
    # Calculate distances per reference
    data <- data %>%
      dplyr::group_by(!!sym(serum_col)) %>%
      dplyr::mutate(
        max_value = max(processed_value, na.rm = TRUE),
        distance = max_value - processed_value
      ) %>%
      dplyr::ungroup()
    
    # Adjust distances for threshold values
    data$distance <- sapply(1:nrow(data), function(i) {
      if (grepl("^<", data$log_value[i])) {
        paste0(">", data$distance[i])
      } else if (grepl("^>", data$log_value[i])) {
        paste0("<", data$distance[i])
      } else {
        data$distance[i]
      }
    })
    
  } else {
    # For IC50, calculate log directly
    data$distance <- sapply(data[[value_col]], function(x) {
      if (grepl("^<", x)) {
        x_num <- as.numeric(sub("<", "", x))
        paste0("<", x_num)
      } else if (grepl("^>", x)) {
        x_num <- as.numeric(sub(">", "", x))
        paste0(">", x_num)
      } else if (is.numeric(as.numeric(x))) {
        as.numeric(x)
      } else {
        NA
      }
    })
  }
  
  # Combine repeated measurements
  long_data <- data %>%
    dplyr::group_by(!!sym(antigen_col), !!sym(serum_col)) %>%
    dplyr::summarise(
      raw_value = reapply_sign(!!sym(value_col), 
                               if(is_titer) {
                                 scale_factor * base^(mean(log(remove_sign(!!sym(value_col))/scale_factor, 
                                                               base = base), na.rm = TRUE))
                               } else {
                                 mean(remove_sign(!!sym(value_col)), na.rm = TRUE)
                               }),
      distance = reapply_sign(distance, 
                              mean(remove_sign(distance), na.rm = TRUE)),
      .groups = 'drop'
    )
  
  # Add back metadata columns if specified
  if (!is.null(metadata_cols)) {
    # Get combined distinct metadata for antigens
    metadata <- data %>%
      dplyr::select(!!sym(antigen_col), !!sym(serum_col), dplyr::any_of(metadata_cols)) %>%
      dplyr::distinct()
    
    # Join metadata based on both antigen and serum columns
    long_data <- long_data %>%
      dplyr::left_join(metadata, by = c(antigen_col, serum_col))
  }
  
  # Remove the non complete rows
  long_data <- na.omit(long_data)
  
  # sort by year to conform with our assumptions used at various places
  if("virusYear" %in% names(long_data)) {
    long_data <- long_data[order(long_data$virusYear), ]
  }
  
  # # Add prefixes if requested
  # if (id_prefix) {
  #   long_data[[antigen_col]] <- paste0("V/", long_data[[antigen_col]])
  #   long_data[[serum_col]] <- paste0("S/", long_data[[serum_col]])
  # }
  
  # Create matrix using long_to_matrix function
  # Determine ordering columns
  virus_year_col <- "virusYear"
  serum_year_col <- "serumYear"
  
  distance_matrix <- long_to_matrix(
    long_data,
    chnames = antigen_col,
    chorder = if("virusYear" %in% names(long_data)) "virusYear" else NULL,
    rnames = serum_col,
    rorder = if("serumYear" %in% names(long_data)) "serumYear" else NULL,
    values_column = "distance",
    rc = FALSE,
    sort = ("virusYear" %in% names(long_data)) || ("serumYear" %in% names(long_data))
  )
  
  # Return both formats
  return(list(
    long = long_data,
    matrix = distance_matrix
  ))
}



#' Validate Antigenic Dataset
#'
#' @description
#' Validates the structure and content of an antigenic assay dataset.
#' Checks for required columns, data types, and valid measurements.
#'
#' @param data Data frame to validate
#' @param antigen_col Character. Name of antigen column
#' @param serum_col Character. Name of serum column
#' @param value_col Character. Name of value column
#' @param metadata_cols Character vector. Optional metadata column names
#' @return Invisibly returns TRUE if validation passes, otherwise stops with error
#' @keywords internal
validate_antigenic_data <- function(data, antigen_col, serum_col, value_col, 
                                  metadata_cols = NULL) {
  # Check column existence
  req_cols <- c(antigen_col, serum_col, value_col)
  if (!all(req_cols %in% names(data))) {
    missing <- setdiff(req_cols, names(data))
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }
  
  # Check metadata columns if specified
  if (!is.null(metadata_cols)) {
    missing_meta <- setdiff(metadata_cols, names(data))
    if (length(missing_meta) > 0) {
      stop("Missing metadata columns: ", paste(missing_meta, collapse = ", "))
    }
  }
  
  # Check for empty/invalid values
  if (any(data[[antigen_col]] == "" | is.na(data[[antigen_col]]))) {
    stop("Empty or NA values found in antigen column")
  }
  if (any(data[[serum_col]] == "" | is.na(data[[serum_col]]))) {
    stop("Empty or NA values found in serum column")
  }
  
  # Validate measurements
  values <- data[[value_col]]
  clean_values <- gsub("^[<>]", "", values)
  numeric_values <- suppressWarnings(as.numeric(clean_values))
  
  if (all(is.na(numeric_values))) {
    stop("No valid numeric measurements found")
  }
  
  if (any(numeric_values < 0, na.rm = TRUE)) {
    stop("Negative values found in measurements")
  }
  
  # Validate threshold indicators
  invalid_thresh <- values[!is.na(values) & !grepl("^[0-9<>]", values)]
  if (length(invalid_thresh) > 0) {
    stop("Invalid threshold indicators found: ", 
         paste(unique(invalid_thresh)[1:5], collapse = ", "), 
         if(length(invalid_thresh) > 5) "...")
  }
  
  invisible(TRUE)
}


#' Prune Distance Data for Network Quality
#'
#' @description
#' Iteratively removes viruses and antibodies with insufficient connections to create a 
#' well-connected network subset. The pruning continues until all remaining points have
#' at least the specified minimum number of connections.
#'
#' @param data Data frame in long format containing:
#'        - Column for viruses/antigens
#'        - Column for antibodies/antisera
#'        - Distance measurements (can contain NAs)
#'        - Optional metadata columns
#' @param virus_col Character name of virus/antigen column
#' @param antibody_col Character name of antibody/antiserum column 
#' @param min_connections Integer minimum required connections per point
#' @param iterations Integer maximum pruning iterations (default 100)
#' @return List containing:
#'   \item{pruned_data}{Data frame of pruned measurements}
#'   \item{stats}{List of pruning statistics including:
#'     \itemize{
#'       \item original_points: Number of points before pruning
#'       \item remaining_points: Number of points after pruning
#'       \item iterations: Number of pruning iterations performed
#'       \item min_connections: Minimum connections in final set
#'     }
#'   }
#' @examples
#' \dontrun{
#' # Basic pruning keeping points with at least 10 connections
#' pruned <- prune_distance_network(hiv_data, 
#'                                 virus_col = "Virus",
#'                                 antibody_col = "Antibody", 
#'                                 min_connections = 10)
#'                                 
#' # Check pruning statistics
#' print(pruned$stats)
#' }
#' @importFrom dplyr %>% filter
#' @importFrom rlang sym
#' @importFrom igraph graph_from_data_frame is_connected components
#' @export
prune_distance_network <- function(data, virus_col, antibody_col,
                                   min_connections, iterations = 100) {
  
  # Input validation
  if (!is.data.frame(data)) {
    stop("data must be a data frame")
  }
  
  if (!all(c(virus_col, antibody_col) %in% names(data))) {
    stop("Specified virus and antibody columns must exist in data")
  }
  
  if (!is.numeric(min_connections) || min_connections < 1) {
    stop("min_connections must be a positive integer")
  }
  
  if (!is.numeric(iterations) || iterations < 1) {
    stop("iterations must be a positive integer")
  }
  
  # Keep track of original data size
  n_original <- length(unique(c(data[[virus_col]], data[[antibody_col]])))
  
  # Initialize pruned data
  pruned_data <- data
  iteration <- 0
  continue_pruning <- TRUE
  
  while(continue_pruning && iteration < iterations) {
    iteration <- iteration + 1
    
    # Count connections per point
    virus_counts <- table(pruned_data[[virus_col]])
    antibody_counts <- table(pruned_data[[antibody_col]])
    
    # Identify points to remove
    low_virus <- names(virus_counts)[virus_counts < min_connections]
    low_antibody <- names(antibody_counts)[antibody_counts < min_connections]
    
    if(length(low_virus) == 0 && length(low_antibody) == 0) {
      continue_pruning <- FALSE
    } else {
      # Remove measurements involving points with too few connections
      pruned_data <- pruned_data %>%
        filter(!(!!sym(virus_col) %in% low_virus) & 
                 !(!!sym(antibody_col) %in% low_antibody))
      
      # Check if we've removed all data
      if(nrow(pruned_data) == 0) {
        stop("No data remains after pruning. Try reducing min_connections.")
      }
    }
  }
  
  # Calculate final statistics
  n_remaining <- length(unique(c(pruned_data[[virus_col]], 
                                 pruned_data[[antibody_col]])))
  
  final_stats <- list(
    original_points = n_original,
    remaining_points = n_remaining,
    iterations = iteration,
    min_connections = min_connections
  )
  
  # Verify network connectivity
  check_connectivity <- function(data, virus_col, antibody_col) {
    # Create adjacency list representation
    edges <- data.frame(
      from = data[[virus_col]],
      to = data[[antibody_col]]
    )
    
    # Create graph
    g <- igraph::graph_from_data_frame(edges, directed = FALSE)
    
    # Check if graph is connected
    is_connected <- igraph::is_connected(g)
    
    if(!is_connected) {
      # Get number of components
      components <- igraph::components(g)
      warning(sprintf(
        "Network has %d disconnected components. Consider increasing min_connections.",
        components$no
      ))
    }
    
    return(is_connected)
  }
  
  # Add connectivity check to stats
  final_stats$is_connected <- check_connectivity(pruned_data, virus_col, antibody_col)
  
  return(list(
    pruned_data = pruned_data,
    stats = final_stats
  ))
}


#' Detect Outliers Using Median Absolute Deviation
#'
#' @description
#' Detects outliers in numeric data using the Median Absolute Deviation (MAD) method.
#' This robust method is less sensitive to extreme values than standard deviation
#' and works well for non-normally distributed data.
#'
#' @details
#' The function:
#' 1. Calculates median and MAD of the data
#' 2. Uses scaled MAD (constant = 1.4826) for normal distribution consistency
#' 3. Identifies points > k MADs from median as outliers 
#' 4. Returns both outlier mask and summary statistics
#'
#' MAD scaling constant 1.4826 is calculated as 1/Phi^(-1)(3/4), where Phi is the 
#' standard normal CDF. This makes MAD consistent with standard deviation for
#' normal distributions.
#'
#' @param data Numeric vector of values to analyze
#' @param k Numeric threshold for outlier detection (default: 3). Points more than
#'        k MADs from median are considered outliers.
#' @param take_log Logical. Whether to log transform data before (and only for) outlier detection (default: FALSE)
#' @return A list containing:
#'   \item{outlier_mask}{Logical vector indicating outliers}
#'   \item{stats}{List containing:
#'     \itemize{
#'       \item median: Median of data
#'       \item mad: Median absolute deviation
#'       \item n_outliers: Number of outliers detected
#'     }
#'   }
#' @examples
#' \dontrun{
#' # Detect outliers in parameter values
#' params <- c(0.01, 0.012, 0.011, 0.1, 0.009, 0.011, 0.15)
#' outliers <- detect_outliers_mad(params)
#' print(outliers$stats$n_outliers) # Number of outliers
#' clean_params <- params[!outliers$outlier_mask] # Remove outliers
#' }
#' @importFrom stats median mad
#' @keywords internal
detect_outliers_mad <- function(data, k = 3, take_log=FALSE) {
  # Extract numeric values and handle thresholds
  process_value <- function(x) {
    if (is.character(x)) {
      if (grepl("^<", x)) {
        as.numeric(sub("<", "", x))
      } else if (grepl("^>", x)) {
        as.numeric(sub(">", "", x)) 
      } else {
        as.numeric(x)
      }
    } else {
      as.numeric(x)
    }
  }
  
  if (!is.numeric(data)) {
    # Convert to numeric
    data <- sapply(data, process_value)
  }
  if (!is.numeric(k) || k <= 0) {
    stop("k must be a positive number")
  }
  
  if(take_log) {
    data <- log(data, base=2)
  }
  
  # Calculate robust statistics
  med <- median(data, na.rm = TRUE)
  mad_val <- stats::mad(data, constant = 1.4826, na.rm = TRUE)
  
  # Identify outliers
  is_outlier <- abs(data - med) > k * mad_val
  
  # Return results
  list(
    outlier_mask = is_outlier,
    stats = list(
      median = med,
      mad = mad_val,
      n_outliers = sum(is_outlier)
    )
  )
}


#' Clean Data by Removing MAD-based Outliers
#'
#' @description
#' Removes outliers from numeric data using the Median Absolute Deviation method.
#' Outliers are replaced with NA values. This function is particularly useful
#' for cleaning parameter tables where each column may contain outliers.
#'
#' @param x Numeric vector to clean
#' @param k Numeric threshold for outlier detection (default: 3)
#' @param take_log Logical. Whether to log transform data before outlier detection (default: FALSE)
#' @return Numeric vector with outliers replaced by NA
#' @examples
#' # Clean parameter values
#' params <- c(0.01, 0.012, 0.011, 0.1, 0.009, 0.011, 0.15)
#' clean_params <- clean_data(params)
#'
#' # Clean multiple parameter columns
#' param_table <- data.frame(
#'   k0 = runif(100),
#'   cooling_rate = runif(100),
#'   c_repulsion = runif(100)
#' )
#' clean_table <- as.data.frame(lapply(param_table, clean_data))
#'
#' @seealso \code{\link{detect_outliers_mad}} for the underlying outlier detection
#' @export
clean_data <- function(x, k = 3, take_log = FALSE) {
  if (!is.numeric(k) || k <= 0) {
    stop("k must be a positive number")
  }
  
  outlier_results <- detect_outliers_mad(x, k = k, take_log)
  x[outlier_results$outlier_mask] <- NA
  return(x)
}