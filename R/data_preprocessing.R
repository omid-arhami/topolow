# Copyright (c) 2024 Omid Arhami omid.arhami@uga.edu
# R/data_preprocessing.R

#' topolow Data Preprocessing Functions
#'


#' Convert List Format Data to Dissimilarity Matrix
#'
#' @description
#' Converts data from long/list format (one measurement per row) to a symmetric
#' dissimilarity matrix. The function handles both similarity and dissimilarity
#' data, with optional conversion from similarity to dissimilarity.
#'
#' @param data Data frame in long format with columns for objects, references, and values.
#' @param object_col Character. Name of the column containing object identifiers.
#' @param reference_col Character. Name of the column containing reference identifiers.
#' @param value_col Character. Name of the column containing measurement values.
#' @param is_similarity Logical. Whether values are similarities (TRUE) or dissimilarities (FALSE).
#'   If TRUE, similarities will be converted to dissimilarities by subtracting from the
#'   maximum value per reference. Default: FALSE.
#'
#' @details
#' The function expects data in long format with at least three columns:
#' - A column for object names
#' - A column for reference names
#' - A column containing the (dis)similarity values
#'
#' When `is_similarity = TRUE`, the function converts similarities to dissimilarities
#' by subtracting each similarity value from the maximum similarity value within
#' each reference group. Threshold indicators (< or >) are handled appropriately
#' and inverted during similarity-to-dissimilarity conversion.
#'
#' @return A symmetric matrix of dissimilarities with row and column names corresponding
#'         to the union of unique objects and references in the data. NA values represent
#'         unmeasured pairs, and the diagonal is set to 0.
#'
#' @examples
#' # Example with dissimilarity data
#' data_dissim <- data.frame(
#'   object = c("A", "B", "A", "C"),
#'   reference = c("X", "X", "Y", "Y"),
#'   dissimilarity = c(2.5, 1.8, 3.0, 4.2)
#' )
#'
#' mat_dissim <- list_to_matrix(
#'   data = data_dissim,
#'   object_col = "object",
#'   reference_col = "reference",
#'   value_col = "dissimilarity",
#'   is_similarity = FALSE
#' )
#'
#' # Example with similarity data (will be converted to dissimilarity)
#' data_sim <- data.frame(
#'   object = c("A", "B", "A", "C"),
#'   reference = c("X", "X", "Y", "Y"),
#'   similarity = c(7.5, 8.2, 7.0, 5.8)
#' )
#'
#' mat_from_sim <- list_to_matrix(
#'   data = data_sim,
#'   object_col = "object",
#'   reference_col = "reference",
#'   value_col = "similarity",
#'   is_similarity = TRUE
#' )
#'
#' @importFrom data.table setDT
#' @importFrom dplyr %>% group_by mutate ungroup summarise
#' @importFrom rlang sym .data
#' @importFrom stats na.omit
#' @export
list_to_matrix <- function(data, object_col, reference_col, value_col,
                           is_similarity = FALSE) {

  # ===== INPUT VALIDATION =====
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }

  required_cols <- c(object_col, reference_col, value_col)
  if (!all(required_cols %in% names(data))) {
    missing <- setdiff(required_cols, names(data))
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }

  # Validate numeric/character columns
  if (!is.character(data[[object_col]]) && !is.factor(data[[object_col]])) {
    stop("Object names column must be character or factor")
  }

  if (!is.character(data[[reference_col]]) && !is.factor(data[[reference_col]])) {
    stop("Reference names column must be character or factor")
  }

  if (!is.numeric(data[[value_col]]) &&
      !all(grepl("^[0-9<>]", stats::na.omit(data[[value_col]])))) {
    stop("Values column must be numeric or contain valid threshold indicators (< or >)")
  }

  if (!is.logical(is_similarity)) {
    stop("is_similarity must be logical")
  }

  # ===== DATA CLEANING =====
  # Clean invalid values
  data <- data[!is.na(data[[value_col]]), ]  # Remove NA values
  data <- data[data[[value_col]] != "", ]    # Remove empty strings

  # Keep only rows where value starts with a number or < or >
  data <- data[grepl("^[0-9<>]", data[[value_col]]), ]

  if (nrow(data) == 0) {
    stop("No valid measurements remaining after cleaning")
  }

  # ===== HELPER FUNCTIONS =====
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

  # ===== SIMILARITY TO DISSIMILARITY CONVERSION =====
  if (is_similarity) {
    # Convert similarity to dissimilarity 
    # Step 1: Convert character values to numeric, preserving threshold signs
    data$signed_value <- sapply(data[[value_col]], function(x) {
      if (grepl("^<", x)) {
        paste0("<", as.numeric(sub("<", "", x)))
      } else if(grepl("^>", x)){
        paste0(">", as.numeric(sub(">", "", x)))
      } else if (is.numeric(as.numeric(x))) {
        as.numeric(x)
      } else {
        NA
      }
    })

    # Step 2: Get raw numeric value for calculations
    data$processed_value <- sapply(data$signed_value, remove_sign)

    # Step 3: Calculate dissimilarity per reference by subtracting from max similarity
    data <- data %>%
      dplyr::group_by(!!rlang::sym(reference_col)) %>%
      dplyr::mutate(
        max_value = max(.data$processed_value, na.rm = TRUE),
        dissimilarity = .data$max_value - .data$processed_value
      ) %>%
      dplyr::ungroup()

    # Step 4: Re-apply threshold signs, inverting them for dissimilarity
    data$final_value <- sapply(1:nrow(data), function(i) {
      if (grepl("^<", data$signed_value[i])) {
        paste0(">", data$dissimilarity[i])
      } else if (grepl("^>", data$signed_value[i])) {
        paste0("<", data$dissimilarity[i])
      } else {
        data$dissimilarity[i]
      }
    })

  } else {
    # If data is already dissimilarity, just format it correctly
    data$final_value <- sapply(data[[value_col]], function(x) {
      if (grepl("^<", x) || grepl("^>", x)) {
        paste0(substr(x, 1, 1), remove_sign(x))
      } else if (is.numeric(as.numeric(x))) {
        as.numeric(x)
      } else {
        NA
      }
    })
  }

  # ===== COMBINE REPEATED MEASUREMENTS =====
  long_data <- data %>%
    dplyr::group_by(!!rlang::sym(object_col), !!rlang::sym(reference_col)) %>%
    dplyr::summarise(
      dissimilarity = reapply_sign(final_value,
                                   mean(remove_sign(final_value), na.rm = TRUE)),
      .groups = 'drop'
    )

  # Remove rows with NAs that might have been introduced
  long_data <- stats::na.omit(long_data)

  # ===== MATRIX CREATION (from titers_list_to_matrix logic) =====

  # Convert to data.table for efficiency
  data.table::setDT(long_data)

  # Get unique point names
  all_points <- unique(c(long_data[[object_col]], long_data[[reference_col]]))

  # Create square matrix with NA values
  n <- length(all_points)
  dissimilarity_matrix <- matrix(NA, nrow = n, ncol = n)
  rownames(dissimilarity_matrix) <- all_points
  colnames(dissimilarity_matrix) <- all_points

  # ===== FILL IN THE DISSIMILARITIES =====
  for (i in seq_len(nrow(long_data))) {
    r   <- long_data[i, get(object_col)]
    c   <- long_data[i, get(reference_col)]
    val <- long_data[i, get("dissimilarity")]

    # Set both matrix elements for symmetry
    dissimilarity_matrix[r, c] <- val
    dissimilarity_matrix[c, r] <- val
  }

  # Set diagonal to 0
  diag(dissimilarity_matrix) <- 0

  return(dissimilarity_matrix)
}


#' Convert Table Format Data to Dissimilarity Matrix
#'
#' @description
#' Converts data from table/matrix format (objects as rows, references as columns)
#' to a symmetric dissimilarity matrix. The function creates a matrix where both
#' rows and columns contain the union of all object and reference names.
#'
#' @param data Matrix or data frame where rownames represent objects, columnnames represent
#'   references, and cells contain (dis)similarity values.
#' @param is_similarity Logical. Whether values are similarities (TRUE) or dissimilarities (FALSE).
#'   If TRUE, similarities will be converted to dissimilarities by subtracting from the
#'   maximum value per column (reference). Default: FALSE.
#'
#' @details
#' The function takes a table where:
#' - Rows represent objects
#' - Columns represent references
#' - Values represent (dis)similarities
#'
#' It creates a symmetric matrix where both rows and columns contain the union of
#' all object names (row names) and reference names (column names). The original
#' measurements are preserved, and the matrix is made symmetric by filling both
#' (i,j) and (j,i) positions with the same value.
#'
#' When `is_similarity = TRUE`, similarities are converted to dissimilarities by
#' subtracting each value from the maximum value in its respective column (reference).
#' Threshold indicators (< or >) are handled and inverted during conversion.
#'
#' @return A symmetric matrix of dissimilarities with row and column names
#'         corresponding to the union of all object and reference names.
#'         NA values represent unmeasured pairs, and the diagonal is set to 0.
#'
#' @examples
#' # Example with dissimilarity data in table format
#' dissim_table <- matrix(c(1.2, 2.1, 3.4, 1.8, 2.9, 4.1),
#'                       nrow = 2, ncol = 3,
#'                       dimnames = list(c("Obj1", "Obj2"),
#'                                     c("Ref1", "Ref2", "Ref3")))
#'
#' mat_dissim <- table_to_matrix(dissim_table, is_similarity = FALSE)
#'
#' # Example with similarity data (will be converted to dissimilarity)
#' sim_table <- matrix(c(8.8, 7.9, 6.6, 8.2, 7.1, 5.9),
#'                    nrow = 2, ncol = 3,
#'                    dimnames = list(c("Obj1", "Obj2"),
#'                                  c("Ref1", "Ref2", "Ref3")))
#'
#' mat_from_sim <- table_to_matrix(sim_table, is_similarity = TRUE)
#'
#' @importFrom stats na.omit
#' @export
table_to_matrix <- function(data, is_similarity = FALSE) {

  # ===== INPUT VALIDATION =====
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("'data' must be a matrix or data frame")
  }

  if (!is.logical(is_similarity)) {
    stop("is_similarity must be logical")
  }

  # Convert to matrix if data frame
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }

  # Validate that data contains valid values
  if (nrow(data) == 0 || ncol(data) == 0) {
    stop("Input data must have at least one row and one column")
  }

  # ===== GET ROW AND COLUMN NAMES =====
  object_names <- rownames(data)
  reference_names <- colnames(data)

  if (is.null(object_names)) {
    object_names <- paste0("Obj", 1:nrow(data))
    rownames(data) <- object_names
  }

  if (is.null(reference_names)) {
    reference_names <- paste0("Ref", 1:ncol(data))
    colnames(data) <- reference_names
  }

  # ===== HELPER FUNCTIONS =====
  remove_sign <- function(x) {
    as.numeric(gsub("[<>]", "", x))
  }

  extract_numeric <- function(x) {
    if (is.character(x)) {
      if (grepl("^<", x) || grepl("^>", x)) {
        as.numeric(gsub("[<>]", "", x))
      } else {
        as.numeric(x)
      }
    } else {
      as.numeric(x)
    }
  }

  # ===== SIMILARITY TO DISSIMILARITY CONVERSION =====
  if (is_similarity) {
    processed_data <- data

    for (j in 1:ncol(data)) {
      col_vals <- data[, j]

      # Extract numeric values for finding maximum
      numeric_vals <- sapply(col_vals, extract_numeric)

      # Find maximum for this reference (column), ignoring NA values
      max_val <- max(numeric_vals, na.rm = TRUE)

      # Convert similarities to dissimilarities for this column
      for (i in 1:nrow(data)) {
        val <- data[i, j]
        if (is.na(val)) {
          processed_data[i, j] <- NA
        } else if (is.character(val)) {
          if (grepl("^<", val)) {
            # < becomes > after conversion
            numeric_part <- as.numeric(gsub("<", "", val))
            processed_data[i, j] <- paste0(">", max_val - numeric_part)
          } else if (grepl("^>", val)) {
            # > becomes < after conversion
            numeric_part <- as.numeric(gsub(">", "", val))
            processed_data[i, j] <- paste0("<", max_val - numeric_part)
          } else {
            # Regular numeric value
            processed_data[i, j] <- max_val - as.numeric(val)
          }
        } else {
          # Numeric value
          processed_data[i, j] <- max_val - as.numeric(val)
        }
      }
    }

    data <- processed_data
  }

  # ===== CREATE SYMMETRIC MATRIX =====
  # Create union of all names
  all_names <- unique(c(object_names, reference_names))
  n <- length(all_names)

  # Create symmetric matrix filled with NA
  symmetric_matrix <- matrix(NA, nrow = n, ncol = n)
  rownames(symmetric_matrix) <- all_names
  colnames(symmetric_matrix) <- all_names

  # ===== FILL IN VALUES FROM ORIGINAL TABLE =====
  for (i in 1:nrow(data)) {
    for (j in 1:ncol(data)) {
      obj_name <- object_names[i]
      ref_name <- reference_names[j]
      value <- data[i, j]

      # Set both symmetric positions
      symmetric_matrix[obj_name, ref_name] <- value
      symmetric_matrix[ref_name, obj_name] <- value
    }
  }

  # Set diagonal to 0
  diag(symmetric_matrix) <- 0

  return(symmetric_matrix)
}


#' Process Raw Antigenic Assay Data
#'
#' @description
#' Processes raw antigenic assay data from data frames into standardized long and matrix
#' formats. Handles both similarity data (like titers, which need conversion to distances) 
#' and direct dissimilarity measurements like IC50. Preserves threshold indicators (<, >) 
#' and handles repeated measurements by averaging.
#'
#' @param data Data frame containing raw data.
#' @param antigen_col Character. Name of column containing virus/antigen identifiers.
#' @param serum_col Character. Name of column containing serum/antibody identifiers.
#' @param value_col Character. Name of column containing measurements (titers or distances).
#' @param is_similarity Logical. Whether values are measures of similarity such as titers 
#'        or binding affinities (TRUE) or dissimilarities like IC50 (FALSE). Default: FALSE.
#' @param metadata_cols Character vector. Names of additional columns to preserve.
#' @param base Numeric. Base for logarithm transformation (default: 2 for similarities, e for dissimilarities).
#' @param scale_factor Numeric. Scale factor for similarities. This is the base value that all other 
#'        dilutions are multiples of. E.g., 10 for HI assay where titers are 10, 20, 40,... Default: 1.
#'
#' @return A list containing two elements:
#'   \item{long}{A `data.frame` in long format with standardized columns, including the original identifiers, processed values, and calculated distances. Any specified metadata is also included.}
#'   \item{matrix}{A numeric `matrix` representing the processed symmetric distance matrix, with antigens and sera on columns and rows.}
#'
#' @details
#' The function handles these key steps:
#' 1. Validates input data and required columns
#' 2. Transforms values to log scale
#' 3. Converts similarities to distances using Smith's method if needed
#' 4. Averages repeated measurements
#' 5. Creates standardized long format
#' 6. Creates symmetric distance matrix
#' 7. Preserves metadata and threshold indicators
#' 8. Preserves virusYear and serumYear columns if present
#' 
#' Input requirements and constraints:
#' * Data frame must contain required columns
#' * Column names must match specified parameters
#' * Values can include threshold indicators (< or >)
#' * Metadata columns must exist if specified
#' * Allowed Year-related column names are "virusYear" and "serumYear"
#'
#' @examples
#' # Example 1: Processing HI titer data (similarities)
#' antigen_data <- data.frame(
#'   virus = c("A/H1N1/2009", "A/H1N1/2010", "A/H1N1/2011", "A/H1N1/2009", "A/H1N1/2010"),
#'   serum = c("anti-2009", "anti-2009", "anti-2009", "anti-2010", "anti-2010"),
#'   titer = c(1280, 640, "<40", 2560, 1280),  # Some below detection limit
#'   cluster = c("A", "A", "B", "A", "A"),
#'   color = c("red", "red", "blue", "red", "red")
#' )
#' 
#' # Process HI titer data (similarities -> distances)
#' results <- process_antigenic_data(
#'   data = antigen_data,
#'   antigen_col = "virus",
#'   serum_col = "serum", 
#'   value_col = "titer",
#'   is_similarity = TRUE,  # Titers are similarities
#'   metadata_cols = c("cluster", "color"),
#'   scale_factor = 10  # Base dilution factor
#' )
#' 
#' # View the long format data
#' print(results$long)
#' # View the distance matrix
#' print(results$matrix)
#' 
#' # Example 2: Processing IC50 data (already dissimilarities)
#' ic50_data <- data.frame(
#'   virus = c("HIV-1", "HIV-2", "HIV-3"),
#'   antibody = c("mAb1", "mAb1", "mAb2"),
#'   ic50 = c(0.05, ">10", 0.2)
#' )
#' 
#' results_ic50 <- process_antigenic_data(
#'   data = ic50_data,
#'   antigen_col = "virus",
#'   serum_col = "antibody",
#'   value_col = "ic50",
#'   is_similarity = FALSE  # IC50 values are dissimilarities
#' )
#'
#' @importFrom dplyr %>% group_by mutate ungroup summarise select distinct left_join any_of
#' @importFrom rlang sym
#' @importFrom stats na.omit
#' @export
process_antigenic_data <- function(data, antigen_col, serum_col, 
                                   value_col,
                                   is_similarity = FALSE, 
                                   metadata_cols = NULL,
                                   base = NULL, 
                                   scale_factor = 1) {
  
  # Validate inputs 
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }

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
  
  if (!is.logical(is_similarity)) {
    stop("is_similarity must be logical")
  }
  
  if (!is.null(base) && (!is.numeric(base) || base <= 0)) {
    stop("base must be NULL or a positive number")
  }
  
  if (!is.numeric(scale_factor) || scale_factor <= 0) {
    stop("scale_factor must be a positive number")
  }

  # Set default base if not provided
  if (is.null(base)) {
    base <- if(is_similarity) 2 else exp(1)
  }

  # Helper function definitions
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
  if (is_similarity) {
    # Convert similarities to log scale
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
    # For dissimilarity, calculate log directly
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
  # Create a column name for the original values to use in reapply_sign
  data$original_values <- data[[value_col]]
  
  long_data <- data %>%
    dplyr::group_by(!!sym(antigen_col), !!sym(serum_col)) %>%
    dplyr::summarise(
      raw_value = reapply_sign(.data$original_values, 
                               if(is_similarity) {
                                 scale_factor * base^(mean(log(remove_sign(.data$original_values)/scale_factor, 
                                                               base = base), na.rm = TRUE))
                               } else {
                                 exp(mean(log(1 + remove_sign(.data$original_values)), na.rm = TRUE)) - 1
                               }),
      distance = reapply_sign(.data$distance, 
                              mean(remove_sign(.data$distance), na.rm = TRUE)),
      .groups = 'drop'
    )
  
  # Add back metadata columns if specified
  if (!is.null(metadata_cols)) {
    # Get combined distinct metadata for antigens and sera
    metadata <- data %>%
      dplyr::select(!!sym(antigen_col), !!sym(serum_col), dplyr::any_of(metadata_cols)) %>%
      dplyr::distinct()
    
    # Join metadata based on both antigen and serum columns
    long_data <- long_data %>%
      dplyr::left_join(metadata, by = c(antigen_col, serum_col))
  }
  
  # Remove the non complete rows
  long_data <- na.omit(long_data)
  
  # Sort by year to conform with our assumptions used at various places
  if("virusYear" %in% names(long_data)) {
    long_data <- long_data[order(long_data$virusYear), ]
  }
  
  # Create matrix using titers_list_to_matrix function
  distance_matrix <- titers_list_to_matrix(
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
#' mat <- titers_list_to_matrix(data, 
#'                      chnames = "antigen",
#'                      rnames = "serum",
#'                      values_column = "distance")
#'                      
#' # With sorting by year
#' mat_sorted <- titers_list_to_matrix(data,
#'                             chnames = "antigen",
#'                             chorder = "year",
#'                             rnames = "serum", 
#'                             rorder = "year",
#'                             values_column = "distance",
#'                             sort = TRUE)
#' @importFrom data.table setDT
#' @importFrom stats na.omit
#' @export
titers_list_to_matrix <- function(data, chnames, chorder = NULL, 
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


#' Clean Data by Removing MAD-based Outliers
#'
#' @description
#' Removes outliers from numeric data using the Median Absolute Deviation (MAD) method.
#' Outliers are replaced with NA values.
#'
#' @param x Numeric vector to clean.
#' @param k Numeric threshold for outlier detection (default: 3).
#' @param take_log Logical. Deprecated parameter. Log transformation should be done before calling this function.
#' @return A numeric vector of the same length as `x`, where detected outliers have been replaced with `NA`.
#' @examples
#' # Clean parameter values
#' params <- c(0.01, 0.012, 0.011, 0.1, 0.009, 0.011, 0.15)
#' clean_params <- clean_data(params)
#'
#' @seealso `detect_outliers_mad` for the underlying outlier detection.
#' @export
clean_data <- function(x, k = 3, take_log = FALSE) {
  if (take_log) {
    lifecycle::deprecate_warn(
      "2.1.0",
      "clean_data(take_log)",
      details = "Log transformation should be done before calling this function."
    )
  }
  if (!is.numeric(k) || k <= 0) {
    stop("k must be a positive number")
  }

  outlier_results <- detect_outliers_mad(x, k = k)
  x[outlier_results$outlier_mask] <- NA
  return(x)
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


#' Convert Coordinates to a Distance Matrix
#'
#' @description
#' Calculates pairwise Euclidean distances between points in a coordinate space.
#' @param positions Matrix or Data Frame of coordinates where rows are points and
#'        columns are dimensions.
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

  p_dist_mat <- as.matrix(stats::dist(positions))

  # Set row and column names if they exist
  if (!is.null(rownames(positions))) {
    rownames(p_dist_mat) <- rownames(positions)
    colnames(p_dist_mat) <- rownames(positions)
  }

  return(p_dist_mat)
}



# Newed
#' Detect Outliers Using Median Absolute Deviation
#'
#' @description
#' Detects outliers in numeric data using the Median Absolute Deviation (MAD) method.
#' This robust method is less sensitive to extreme values than standard deviation
#' and works well for non-normally distributed data.
#'
#' @details
#' The function calculates the median and MAD of the data and identifies points
#' that are more than `k` MADs from the median as outliers.
#'
#' @param data Numeric vector of values to analyze
#' @param k Numeric threshold for outlier detection (default: 3).
#' @return A list containing:
#'   \item{outlier_mask}{Logical vector indicating outliers}
#'   \item{stats}{List containing:
#'     \itemize{
#'       \item median: Median of data
#'       \item mad: Median absolute deviation
#'       \item n_outliers: Number of outliers detected
#'     }
#'   }#' @importFrom stats median mad
#' @keywords internal
detect_outliers_mad <- function(data, k = 3) {
  # Helper to extract numeric values from character strings with thresholds
  process_value <- function(x) {
    if (is.character(x)) {
      if (grepl("^<", x)) {
        suppressWarnings(as.numeric(sub("<", "", x)))
      } else if (grepl("^>", x)) {
        suppressWarnings(as.numeric(sub(">", "", x)))
      } else {
        suppressWarnings(as.numeric(x))
      }
    } else {
      as.numeric(x)
    }
  }

  if (!is.numeric(data)) {
    # Convert to numeric if not already
    data <- sapply(data, process_value)
  }
  if (!is.numeric(k) || k <= 0) {
    stop("k must be a positive number")
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
      n_outliers = sum(is_outlier, na.rm = TRUE)
    )
  )
}


#' Log Transform Parameter Samples
#'
#' @description
#' Reads parameter samples from a CSV file and applies a log transformation to
#' specified parameter columns (e.g., N, k0, cooling_rate, c_repulsion).
#'
#' **Note**: As of version 2.0.0, this function is primarily for backward compatibility
#' with existing parameter files. The `initial_parameter_optimization()` function now
#' returns log-transformed parameters directly, eliminating the need for this separate
#' transformation step in the normal workflow.
#'
#' @details
#' This function is maintained for users who have existing parameter files from
#' older versions of the package or who need to work with parameter files that
#' contain original-scale parameters. In the current workflow:
#' 
#' - `initial_parameter_optimization()` --> returns log-transformed parameters directly
#' - `run_adaptive_sampling()` --> works with log-transformed parameters
#' - `euclidean_embedding()` --> works with original-scale parameters
#' 
#' If you are working with the current workflow (using `Euclidify()` or calling
#' `initial_parameter_optimization()` directly), you typically do not need to call
#' this function.
#'
#' @param samples_file Character. Path to the CSV file containing the parameter samples.
#' @param output_file Character. Optional path to save the transformed data as a new CSV file.
#' @return A `data.frame` with the log-transformed parameters. If `output_file` is
#'         specified, the data frame is also written to a file and returned invisibly.
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
#' @note 
#' **Backward Compatibility Note**: This function is maintained for compatibility
#' with existing workflows and parameter files. For new workflows, consider using
#' `initial_parameter_optimization()` which returns log-transformed parameters directly.
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

  # Check which parameters exist to be transformed
  params_to_transform <- c("N", "k0", "cooling_rate", "c_repulsion")
  existing_params <- intersect(names(samples), params_to_transform)

  if (length(existing_params) == 0) {
    message("No parameters found to transform or all parameters are already transformed.")
    return(samples)
  }

  # Validate that parameter columns are numeric and positive
  for (param in existing_params) {
    if (!is.numeric(samples[[param]])) {
      samples[[param]] <- suppressWarnings(as.numeric(samples[[param]]))
      if (any(is.na(samples[[param]]))) {
        stop("Non-numeric values found in column: ", param)
      }
    }
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

  # Attempt to reorder to a standard format if columns exist
  desired_order <- c(
    "log_N",
    "log_k0",
    "log_cooling_rate",
    "log_c_repulsion",
    "Holdout_MAE",
    "NLL"
  )
  
  # Only reorder if all desired columns are present
  if(all(desired_order %in% names(samples))) {
      samples <- samples[, desired_order, drop = FALSE]
  }


  # Write output if a path is provided
  if (!is.null(output_file)) {
      if (!is.character(output_file) || length(output_file) != 1) {
          stop("'output_file' must be a single character string path.", call. = FALSE)
      }
      tryCatch({
        utils::write.csv(samples, file = output_file, row.names = FALSE)
        message("Log transformed parameters: ", paste(existing_params, collapse = ", "))
        message("Output saved to: ", output_file)
      }, error = function(e) {
        stop("Error writing output file: ", e$message)
      })
      return(invisible(samples))
  }

  # Return the transformed data frame if not writing to file
  return(samples)
}




#' Prune Sparse Dissimilarity Matrix to Well-Connected Subset
#'
#' @description
#' Iteratively removes poorly connected points from a sparse dissimilarity matrix
#' to create a well-connected, denser subset suitable for optimization and subsampling.
#' This is particularly useful for very sparse datasets (>95% missing) where random
#' subsampling would create disconnected graphs.
#'
#' @param dissimilarity_matrix Square symmetric dissimilarity matrix. Can contain
#'   NA values and threshold indicators (< or >).
#' @param min_connections Integer. Minimum number of observed dissimilarities 
#'   required per point. Points with fewer connections are iteratively removed.
#'   Default: 4. Higher values create denser but smaller subsets.
#' @param target_completeness Numeric. Target network completeness (0-1) to achieve.
#'   If specified, overrides min_connections and adaptively increases the threshold
#'   until target is reached. Default: NULL.
#' @param max_iterations Integer. Maximum pruning iterations to prevent infinite loops.
#'   Default: 100.
#' @param min_points Integer. Minimum number of points to retain. Stops pruning if
#'   fewer points remain. Default: 10.
#' @param ensure_connected Logical. If TRUE, verifies final subset is connected and
#'   warns if not. Requires igraph package. Default: TRUE.
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return A list containing:
#'   \item{pruned_matrix}{Matrix. The pruned dissimilarity matrix.}
#'   \item{kept_indices}{Integer vector. Row/column indices of kept points.}
#'   \item{kept_names}{Character vector. Names of kept points (if original had names).}
#'   \item{removed_indices}{Integer vector. Indices of removed points.}
#'   \item{stats}{List of pruning statistics:
#'     \itemize{
#'       \item \code{original_size}: Original matrix dimensions
#'       \item \code{final_size}: Final matrix dimensions  
#'       \item \code{points_removed}: Number of points removed
#'       \item \code{percent_removed}: Percentage of points removed
#'       \item \code{original_completeness}: Original network completeness
#'       \item \code{final_completeness}: Final network completeness
#'       \item \code{original_measurements}: Original number of observed values
#'       \item \code{final_measurements}: Final number of observed values
#'       \item \code{iterations}: Number of pruning iterations
#'       \item \code{min_connections_used}: Final min_connections threshold
#'       \item \code{is_connected}: Whether final subset is connected (if checked)
#'       \item \code{n_components}: Number of connected components (if checked)
#'     }
#'   }
#'
#' @details
#' The function works by:
#' 1. Counting observed measurements per point (row/column)
#' 2. Identifying points below the threshold
#' 3. Removing those points and their measurements
#' 4. Repeating until all remaining points meet the threshold
#' 5. Optionally verifying connectivity
#'
#' **Choosing min_connections:**
#' - For very sparse data, start with min_connections = 3-10
#' - For target completeness, use target_completeness instead
#'
#' **Adaptive thresholding with target_completeness:**
#' If target_completeness is specified, the function:
#' 1. Starts with min_connections = 4
#' 2. Prunes the matrix
#' 3. Checks if completeness target is met
#' 4. If not, increases min_connections by 1 and repeats
#' 5. Stops when target is reached or min_points is hit
#'
#' @examples
#' # Create a sparse matrix
#' set.seed(123)
#' n <- 1000
#' mat <- matrix(NA, n, n)
#' diag(mat) <- 0
#' # Add only 1% observations
#' n_obs <- floor(n * n * 0.15)
#' indices <- sample(which(upper.tri(mat)), n_obs)
#' mat[indices] <- runif(n_obs, 0, 10)
#' mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
#'
#' # Prune with minimum connections
#' result <- prune_sparse_matrix(mat, min_connections = 2)
#' print(result$stats)
#'
#' # Prune to target completeness
#' result2 <- prune_sparse_matrix(mat, target_completeness = 0.05)
#' print(result2$stats)
#'
#' @seealso
#' \code{\link{check_matrix_connectivity}} for checking connectivity,
#' \code{\link{analyze_network_structure}} for network analysis
#'
#' @export
prune_sparse_matrix <- function(dissimilarity_matrix,
                                min_connections = 4,
                                target_completeness = NULL,
                                max_iterations = 100,
                                min_points = 10,
                                ensure_connected = TRUE,
                                verbose = TRUE) {
  
  # ==========================================================================
  # INPUT VALIDATION
  # ==========================================================================
  if (!is.matrix(dissimilarity_matrix)) {
    stop("dissimilarity_matrix must be a matrix")
  }
  
  if (nrow(dissimilarity_matrix) != ncol(dissimilarity_matrix)) {
    stop("dissimilarity_matrix must be square")
  }
  
  original_size <- nrow(dissimilarity_matrix)
  
  if (original_size < min_points) {
    stop(sprintf("Matrix size (%d) is smaller than min_points (%d)",
                 original_size, min_points))
  }
  
  if (!is.numeric(min_connections) || min_connections < 1) {
    stop("min_connections must be a positive integer")
  }
  
  if (!is.null(target_completeness)) {
    if (!is.numeric(target_completeness) || 
        target_completeness <= 0 || target_completeness > 1) {
      stop("target_completeness must be between 0 and 1")
    }
  }
  
  # Check for igraph if connectivity check is requested
  if (ensure_connected && !requireNamespace("igraph", quietly = TRUE)) {
    warning("igraph package not available. Skipping connectivity check.")
    ensure_connected <- FALSE
  }
  
  # ==========================================================================
  # CALCULATE ORIGINAL STATISTICS
  # ==========================================================================
  original_measurements <- sum(!is.na(dissimilarity_matrix) & 
                                 row(dissimilarity_matrix) != col(dissimilarity_matrix)) / 2
  original_completeness <- (2 * original_measurements) / (original_size * (original_size - 1))
  
  if (verbose) {
    cat("\n", rep("=", 80), "\n", sep = "")
    cat("SPARSE MATRIX PRUNING\n")
    cat(rep("=", 80), "\n", sep = "")
    cat(sprintf("\nOriginal matrix: %d x %d\n", original_size, original_size))
    cat(sprintf("Original measurements: %d\n", original_measurements))
    cat(sprintf("Original completeness: %.3f (%.1f%%)\n", 
                original_completeness, original_completeness * 100))
  }
  
  # ==========================================================================
  # ADAPTIVE THRESHOLDING (if target_completeness specified)
  # ==========================================================================
  if (!is.null(target_completeness)) {
    if (verbose) {
      cat(sprintf("\nUsing adaptive thresholding to reach %.1f%% completeness...\n",
                  target_completeness * 100))
    }
    
    current_threshold <- 4
    threshold_increment <- 2
    max_threshold_attempts <- 10
    
    for (attempt in 1:max_threshold_attempts) {
      if (verbose) {
        cat(sprintf("\nAttempt %d: min_connections = %d\n", 
                    attempt, current_threshold))
      }
      
      # Try pruning with current threshold
      temp_result <- prune_sparse_matrix(
        dissimilarity_matrix = dissimilarity_matrix,
        min_connections = current_threshold,
        target_completeness = NULL,  # Disable recursion
        max_iterations = max_iterations,
        min_points = min_points,
        ensure_connected = FALSE,  # Check only at end
        verbose = FALSE
      )
      
      # Check if target reached
      if (temp_result$stats$final_completeness >= target_completeness) {
        if (verbose) {
          cat(sprintf("  Target reached: %.3f >= %.3f\n",
                      temp_result$stats$final_completeness,
                      target_completeness))
        }
        
        # Do final connectivity check if requested
        if (ensure_connected) {
          connectivity <- check_matrix_connectivity(
            temp_result$pruned_matrix,
            min_completeness = target_completeness
          )
          temp_result$stats$is_connected <- connectivity$is_connected
          temp_result$stats$n_components <- connectivity$n_components
        }
        
        return(temp_result)
      }
      
      # Check if we've hit minimum size
      if (temp_result$stats$final_size <= min_points) {
        warning(sprintf(
          paste0("Reached min_points (%d) before target completeness. ",
                 "Final completeness: %.3f, Target: %.3f"),
          min_points, temp_result$stats$final_completeness, target_completeness
        ))
        return(temp_result)
      }
      
      if (verbose) {
        cat(sprintf("  Completeness: %.3f < %.3f, increasing threshold...\n",
                    temp_result$stats$final_completeness, target_completeness))
      }
      
      current_threshold <- current_threshold + threshold_increment
    }
    
    stop("Could not reach target completeness after maximum attempts")
  }
  
  # ==========================================================================
  # STANDARD PRUNING (fixed min_connections)
  # ==========================================================================
  if (verbose) {
    cat(sprintf("\nPruning with min_connections = %d\n", min_connections))
    cat(rep("-", 80), "\n", sep = "")
  }
  
  # Initialize
  current_matrix <- dissimilarity_matrix
  kept_indices <- seq_len(original_size)
  iteration <- 0
  continue_pruning <- TRUE
  
  while (continue_pruning && iteration < max_iterations) {
    iteration <- iteration + 1
    current_size <- nrow(current_matrix)
    
    # Count connections per point (excluding diagonal)
    connection_counts <- apply(current_matrix, 1, function(row) {
      sum(!is.na(row) & row != row[1])  # Exclude diagonal (self-connection)
    })
    
    # Identify poorly connected points
    poorly_connected <- which(connection_counts < min_connections)
    n_to_remove <- length(poorly_connected)
    
    if (n_to_remove == 0) {
      # All points meet threshold
      continue_pruning <- FALSE
      if (verbose) {
        cat(sprintf("Iteration %d: All %d points have >= %d connections\n",
                    iteration, current_size, min_connections))
      }
    } else if ((current_size - n_to_remove) < min_points) {
      # Would drop below minimum
      warning(sprintf(
        paste0("Stopping at iteration %d: removing %d points would leave only %d ",
               "points (min_points = %d)"),
        iteration, n_to_remove, current_size - n_to_remove, min_points
      ))
      continue_pruning <- FALSE
    } else {
      # Remove poorly connected points
      keep_mask <- connection_counts >= min_connections
      current_matrix <- current_matrix[keep_mask, keep_mask, drop = FALSE]
      kept_indices <- kept_indices[keep_mask]
      
      if (verbose && iteration %% 10 == 0) {
        cat(sprintf("Iteration %d: Removed %d points, %d remaining\n",
                    iteration, n_to_remove, nrow(current_matrix)))
      }
    }
    
    # Safety check
    if (nrow(current_matrix) == 0) {
      stop("All points removed during pruning. Try lower min_connections.")
    }
  }
  
  if (iteration >= max_iterations) {
    warning(sprintf("Reached max_iterations (%d) before convergence", max_iterations))
  }
  
  # ==========================================================================
  # CALCULATE FINAL STATISTICS
  # ==========================================================================
  final_size <- nrow(current_matrix)
  final_measurements <- sum(!is.na(current_matrix) & 
                             row(current_matrix) != col(current_matrix)) / 2
  final_completeness <- (2 * final_measurements) / (final_size * (final_size - 1))
  
  removed_indices <- setdiff(seq_len(original_size), kept_indices)
  
  stats <- list(
    original_size = original_size,
    final_size = final_size,
    points_removed = length(removed_indices),
    percent_removed = (length(removed_indices) / original_size) * 100,
    original_completeness = original_completeness,
    final_completeness = final_completeness,
    original_measurements = original_measurements,
    final_measurements = final_measurements,
    completeness_increase = final_completeness / original_completeness,
    iterations = iteration,
    min_connections_used = min_connections
  )
  
  # ==========================================================================
  # CONNECTIVITY CHECK
  # ==========================================================================
  if (ensure_connected) {
    if (verbose) cat("\nChecking connectivity...\n")
    
    connectivity <- tryCatch({
      check_matrix_connectivity(current_matrix, min_completeness = 0.01)
    }, error = function(e) {
      warning("Connectivity check failed: ", e$message)
      list(is_connected = NA, n_components = NA)
    })
    
    stats$is_connected <- connectivity$is_connected
    stats$n_components <- connectivity$n_components
    
    if (!is.na(connectivity$is_connected) && !connectivity$is_connected) {
      warning(sprintf(
        "Final matrix has %d disconnected components. Consider increasing min_connections.",
        connectivity$n_components
      ))
    }
  }
  
  # ==========================================================================
  # PREPARE RETURN OBJECT
  # ==========================================================================
  kept_names <- NULL
  if (!is.null(rownames(dissimilarity_matrix))) {
    kept_names <- rownames(dissimilarity_matrix)[kept_indices]
  }
  
  if (verbose) {
    cat("\n", rep("=", 80), "\n", sep = "")
    cat("PRUNING COMPLETE\n")
    cat(rep("=", 80), "\n", sep = "")
    cat(sprintf("\nFinal matrix: %d x %d (%d points removed, %.1f%%)\n",
                final_size, final_size, stats$points_removed, stats$percent_removed))
    cat(sprintf("Final measurements: %d\n", final_measurements))
    cat(sprintf("Final completeness: %.3f (%.1f%%) - %.1fx increase\n",
                final_completeness, final_completeness * 100, stats$completeness_increase))
    if (!is.na(stats$is_connected)) {
      cat(sprintf("Connected: %s", stats$is_connected))
      if (!stats$is_connected) {
        cat(sprintf(" (%d components)", stats$n_components))
      }
      cat("\n")
    }
    cat("\n")
  }
  
  return(list(
    pruned_matrix = current_matrix,
    kept_indices = kept_indices,
    kept_names = kept_names,
    removed_indices = removed_indices,
    stats = stats
  ))
}