# Copyright (c) 2024 Omid Arhami omid.arhami@uga.edu
# R/core.R

#' TopoLow Core Functions

#' Vectorized Processing of Dissimilarity Matrix for Convergence Error Calculations
#'
#' @description
#' Efficiently processes elements of the dissimilarity matrix for calculating convergence error
#' using pre-processed numeric representations of thresholds. This optimized version
#' eliminates expensive string operations during optimization.
#'
#' @details
#' This function handles threshold logic for convergence error calculation by using
#' pre-processed numeric matrices:
#'
#' - For "greater than" thresholds (threshold_mask = 1): Returns the numeric value if the
#'   calculated distance is less than the threshold, otherwise returns NA
#'
#' - For "less than" thresholds (threshold_mask = -1): Returns the numeric value if the
#'   calculated distance is greater than the threshold, otherwise returns NA
#'
#' - For regular values (threshold_mask = 0): Returns the numeric value
#'
#' This function operates on entire matrices at once using vectorized operations,
#' which is significantly faster than processing each element individually.
#'
#' @param distances_numeric Numeric matrix. The numeric dissimilarity values (without threshold
#'        indicators)
#' @param threshold_mask Integer matrix. Codes representing threshold types:
#'        1 for "greater than" (>), -1 for "less than" (<), or 0 for exact values
#' @param p_dist_mat Numeric matrix. The calculated distance matrix to compare against
#'
#' @return Numeric matrix with processed distance values. Elements where threshold
#'         conditions are not satisfied will contain NA.
#'
#' @keywords internal
vectorized_process_distance_matrix <- function(distances_numeric, threshold_mask, p_dist_mat) {
  # Create result matrix (defaults to NA)
  result <- matrix(NA, nrow = nrow(distances_numeric), ncol = ncol(distances_numeric))

  # Find indices for each threshold type
  gt_indices <- which(threshold_mask == 1)  # Greater than (>)
  lt_indices <- which(threshold_mask == -1) # Less than (<)
  normal_indices <- which(threshold_mask == 0) # Normal values

  # Process greater than thresholds
  if (length(gt_indices) > 0) {
    valid <- !is.na(distances_numeric[gt_indices]) &
      !is.na(p_dist_mat[gt_indices]) &
      p_dist_mat[gt_indices] < distances_numeric[gt_indices]

    if (any(valid)) {
      result[gt_indices[valid]] <- distances_numeric[gt_indices[valid]]
    }
  }

  # Process less than thresholds
  if (length(lt_indices) > 0) {
    valid <- !is.na(distances_numeric[lt_indices]) &
      !is.na(p_dist_mat[lt_indices]) &
      p_dist_mat[lt_indices] > distances_numeric[lt_indices]

    if (any(valid)) {
      result[lt_indices[valid]] <- distances_numeric[lt_indices[valid]]
    }
  }

  # Process normal values
  if (length(normal_indices) > 0) {
    result[normal_indices] <- distances_numeric[normal_indices]
    result[is.infinite(result)] <- NA
  }

  return(result)
}


#' Main topolow algorithm implementation
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' topolow (topological stochastic pairwise reconstruction for Euclidean embedding) optimizes
#' point positions in an N-dimensional space to match a target dissimilarity matrix.
#' The algorithm uses a physics-inspired approach with spring and repulsive forces
#' to find optimal point configurations while handling missing and thresholded measurements.
#'
#' @details
#' The algorithm iteratively updates point positions using:
#' * Spring forces between points with measured dissimilarities.
#' * Repulsive forces between points without measurements.
#' * Conditional forces for thresholded measurements (< or >).
#' * An adaptive spring constant that decays over iterations.
#' * Convergence monitoring based on relative error change.
#' * Automatic matrix reordering to optimize convergence.
#' Consider if downstream analyses depend on specific point ordering: The order of points in 
#' the output is adjusted to put high-dissimilarity points in the opposing ends.
#'
#' This function replaces the deprecated [create_topolow_map()]. The core algorithm
#' is identical, but includes performance improvements and enhanced validation.
#'
#' @param dissimilarity_matrix Matrix. A square, symmetric dissimilarity matrix. Can contain
#'        NA values for missing measurements and character strings with < or > prefixes for
#'        thresholded measurements.
#' @param ndim Integer. Number of dimensions for the embedding space.
#' @param mapping_max_iter Integer. Maximum number of map optimization iterations.
#' @param k0 Numeric. Initial spring constant controlling spring forces.
#' @param cooling_rate Numeric. Rate of spring constant decay per iteration (0 < cooling_rate < 1).
#' @param c_repulsion Numeric. Repulsion constant controlling repulsive forces.
#' @param relative_epsilon Numeric. Convergence threshold for relative change in error.
#'        Default is 1e-4.
#' @param convergence_counter Integer. Number of iterations below threshold before declaring
#'        convergence. Default is 5.
#' @param initial_positions Matrix or NULL. Optional starting coordinates. If NULL,
#'        random initialization is used. Matrix should have nrow = nrow(dissimilarity_matrix)
#'        and ncol = ndim.
#' @param write_positions_to_csv Logical. Whether to save point positions to a CSV file.
#'        Default is FALSE.
#' @param output_dir Character. Directory to save the CSV file. Required if
#'        `write_positions_to_csv` is TRUE.
#' @param verbose Logical. Whether to print progress messages. Default is FALSE.
#'
#' @return A `list` object of class `topolow`. This list contains the results of the
#'   optimization and includes the following components:
#' \itemize{
#'   \item `positions`: A `matrix` of the optimized point coordinates in the N-dimensional space.
#'   \item `est_distances`: A `matrix` of the Euclidean distances between points in the final optimized configuration.
#'   \item `mae`: The final Mean Absolute Error between the target dissimilarities and the estimated distances.
#'   \item `iter`: The total number of iterations performed before the algorithm terminated.
#'   \item `parameters`: A `list` containing the input parameters used for the optimization run.
#'   \item `convergence`: A `list` containing the final convergence status, including a logical `achieved` flag and the final `error` value.
#' }
#'
#' @examples
#' # Create a simple dissimilarity matrix
#' dist_mat <- matrix(c(0, 2, 3, 2, 0, 4, 3, 4, 0), nrow=3)
#'
#' # Run topolow in 2D
#' result <- euclidean_embedding(
#'   dissimilarity_matrix = dist_mat,
#'   ndim = 2,
#'   mapping_max_iter = 100,
#'   k0 = 1.0,
#'   cooling_rate = 0.001,
#'   c_repulsion = 0.01,
#'   verbose = FALSE
#' )
#'
#' # View results
#' head(result$positions)
#' print(result$mae)
#'
#' # Example with thresholded measurements
#' thresh_mat <- matrix(c(0, ">2", 3, ">2", 0, "<5", 3, "<5", 0), nrow=3)
#' result_thresh <- euclidean_embedding(
#'   dissimilarity_matrix = thresh_mat,
#'   ndim = 2,
#'   mapping_max_iter = 50,
#'   k0 = 0.5,
#'   cooling_rate = 0.01,
#'   c_repulsion = 0.001
#' )
#'
#' @seealso [create_topolow_map()] for the deprecated predecessor function.
#'
#' @importFrom stats runif dist hclust
#' @importFrom utils write.csv
#' @export
euclidean_embedding <- function(dissimilarity_matrix,
                               ndim,
                               mapping_max_iter,
                               k0,
                               cooling_rate,
                               c_repulsion,
                               relative_epsilon = 1e-4,
                               convergence_counter = 5,
                               initial_positions = NULL,
                               write_positions_to_csv = FALSE,
                               output_dir,
                               verbose = FALSE) {

  # Input validation with informative messages
  if (!is.matrix(dissimilarity_matrix)) {
    stop("dissimilarity_matrix must be a matrix")
  }
  if (nrow(dissimilarity_matrix) != ncol(dissimilarity_matrix)) {
    stop("dissimilarity_matrix must be square")
  }
  # Check if there are any finite non-zero dissimilarities
  finite_dissim <- suppressWarnings(as.numeric(dissimilarity_matrix))
  finite_dissim[is.infinite(finite_dissim)] <- NA
  if (sum(!is.na(finite_dissim) & finite_dissim != 0) == 0) {
    warning("No finite non-zero dissimilarities found. Results may be unreliable.")
  }
  if (!is.numeric(ndim) || ndim < 1 || ndim != round(ndim)) {
    stop("ndim must be a positive integer")
  }
  if (!is.numeric(mapping_max_iter) || mapping_max_iter < 1 || mapping_max_iter != round(mapping_max_iter)) {
    stop("mapping_max_iter must be a positive integer")
  }
  if (!is.numeric(k0) || k0 <= 0) {
    stop("k0 must be a positive number")
  }
  if (k0 > 30) {
    warning("High k0 value (> 30) may lead to instability")
  }
  if (!is.numeric(cooling_rate) || cooling_rate <= 0 || cooling_rate >= 1) {
    stop("cooling_rate must be between 0 and 1")
  }
  if (!is.numeric(c_repulsion) || c_repulsion <= 0) {
    stop("c_repulsion must be a positive number")
  }
  if (!is.numeric(relative_epsilon) || relative_epsilon <= 0) {
    stop("relative_epsilon must be a positive number")
  }
  if (!is.numeric(convergence_counter) ||
      convergence_counter < 1 ||
      convergence_counter != round(convergence_counter)) {
    stop("convergence_counter must be a positive integer")
  }

  # Validate initial positions if provided
  if (!is.null(initial_positions)) {
    if (!is.matrix(initial_positions)) {
      stop("initial_positions must be a matrix")
    }
    if (nrow(initial_positions) != nrow(dissimilarity_matrix)) {
      stop("initial_positions must have same number of rows as dissimilarity_matrix")
    }
    if (ncol(initial_positions) != ndim) {
      stop("initial_positions must have ndim columns")
    }
  }

  # Initialize variables
  n <- nrow(dissimilarity_matrix)

  if (n < 2) {
    stop("dissimilarity_matrix must have at least 2 rows/columns")
  }


  ###=== Reorder rows and columns so largest values are furthest from diagonal
  if (n > 1) {
    # Extract numeric values for clustering, handling threshold indicators
    numeric_matrix <- matrix(NA, n, n)
    for(i in 1:n) {
      for(j in 1:n) {
        val <- dissimilarity_matrix[i,j]
        if(!is.na(val)) {
          if(is.character(val)) {
            # Remove < or > prefix and extract numeric value
            numeric_matrix[i,j] <- as.numeric(gsub("^[<>]", "", val))
          } else {
            numeric_matrix[i,j] <- as.numeric(val)
          }
        }
      }
    }

    # Create spectral ordering to concentrate largest values in corners
    tryCatch({
      # Calculate average dissimilarity for each point using only non-NA values
      avg_dissim <- numeric(n)
      for(i in 1:n) {
        # Get all values for this point (excluding diagonal)
        row_vals <- numeric_matrix[i, -i]
        col_vals <- numeric_matrix[-i, i]
        all_vals <- c(row_vals, col_vals)

        # Calculate average using only non-NA values
        non_na_vals <- all_vals[!is.na(all_vals)]
        if(length(non_na_vals) > 0) {
          avg_dissim[i] <- mean(non_na_vals)
        } else {
          avg_dissim[i] <- 0  # If no measurements, put in middle
        }
      }

      # Check if we have meaningful averages (not all zeros)
      if(sum(avg_dissim > 0) > 1) {
        # Order points by average dissimilarity (descending)
        # This puts high-dissimilarity points at extremes, low at center
        new_order <- order(avg_dissim)

        # Apply ordering to the original matrix (with NAs and thresholds intact)
        dissimilarity_matrix <- dissimilarity_matrix[new_order, new_order]
        if (verbose) cat("Matrix reordered for spectral pattern (largest values in corners)\n")
      } else {
        if (verbose) cat("Insufficient data for meaningful spectral ordering\n")
      }
    }, error = function(e) {
      # If ordering fails, skip reordering
      if (verbose) cat("Could not reorder matrix:", e$message, "\n")
    })
  }

  # If initial positions provided, ensure they match the reordered matrix
  if (!is.null(initial_positions)) {
    if (!is.null(rownames(initial_positions)) &&
        !is.null(rownames(dissimilarity_matrix))) {
      if (!identical(rownames(initial_positions), rownames(dissimilarity_matrix))) {
        # Reorder initial_positions to match
        initial_positions <- initial_positions[rownames(dissimilarity_matrix), ]
      }
    }
  }
  #====

  # Pre-compute NA status once
  is_na_matrix <- is.na(dissimilarity_matrix)
  node_degrees <- rowSums(!is_na_matrix)
  node_degrees_1 <- node_degrees + 1

  distances_numeric <- matrix(Inf, n, n)
  threshold_mask <- matrix(0, n, n)  # 0 = number, 1 = >, -1 = <

  for(i in 1:n) {
    for(j in 1:n) {
      if(!is.na(dissimilarity_matrix[i,j])) {
        if(is.character(dissimilarity_matrix[i,j])) {
          if(startsWith(dissimilarity_matrix[i,j], ">")) {
            distances_numeric[i,j] <- as.numeric(sub(">", "", dissimilarity_matrix[i,j]))
            threshold_mask[i,j] <- 1
          } else if(startsWith(dissimilarity_matrix[i,j], "<")) {
            distances_numeric[i,j] <- as.numeric(sub("<", "", dissimilarity_matrix[i,j]))
            threshold_mask[i,j] <- -1
          } else {
            distances_numeric[i,j] <- as.numeric(dissimilarity_matrix[i,j])
          }
        } else {
          distances_numeric[i,j] <- dissimilarity_matrix[i,j]
        }
      }
    }
  }

  # Initialize positions
  positions <- if(is.null(initial_positions)) {
    dissimilarity_matrix_numeric <- suppressWarnings(as.numeric(as.character(dissimilarity_matrix)))
    dissimilarity_matrix_numeric[is.na(dissimilarity_matrix_numeric)] <- NA
    init_step <- max(dissimilarity_matrix_numeric, na.rm=TRUE) / n
    # random initialization
    random_steps <- matrix(stats::runif((n-1) * ndim, 0, 2*init_step), nrow=n-1, ncol=ndim)
    # First row stays zero, subsequent rows are cumulative sums
    pos <- rbind(matrix(0, 1, ndim), apply(random_steps, 2, cumsum))
    pos
  } else {
    initial_positions
  }

  rownames(positions) <- rownames(dissimilarity_matrix)

  # Initialize convergence tracking
  convergence_error0 <- 1e12
  prev_convergence_error <- convergence_error0
  k <- k0

  # OPTIMIZATION 3: Create point pairs once
  point_pairs <- matrix(0, n*(n-1)/2, 2)
  pair_idx <- 1
  for(i in 1:(n-1)) {
    for(j in (i+1):n) {
      point_pairs[pair_idx, ] <- c(i, j)
      pair_idx <- pair_idx + 1
    }
  }

  # OPTIMIZATION 4: Pre-compute which pairs have measured dissimilarities
  has_measurement <- matrix(FALSE, nrow(point_pairs), 1)
  for(idx in 1:nrow(point_pairs)) {
    i <- point_pairs[idx, 1]
    j <- point_pairs[idx, 2]
    has_measurement[idx] <- distances_numeric[i, j] != Inf
  }

  ################## Main optimization loop
  for(iter in 1:mapping_max_iter) {
    k_2 <- k / 2
    stop <- FALSE

    # Randomize point pair ordering for this iteration
    random_order <- sample(nrow(point_pairs))
    shuffled_pairs <- point_pairs[random_order,]
    shuffled_has_measurement <- has_measurement[random_order]

    # Calculate forces between pairs in random order
    for(pair_idx in 1:nrow(shuffled_pairs)) {
      i <- shuffled_pairs[pair_idx, 1]
      j <- shuffled_pairs[pair_idx, 2]

      delta <- positions[j,] - positions[i,]
      distance <- sqrt(sum(delta^2))
      distance_01 <- distance + 0.01

      # Use pre-computed measurement status
      if(shuffled_has_measurement[pair_idx]) {
        # Inline logic for processing ideal distance based on threshold mask
        if(threshold_mask[i, j] == 1) {  # ">" case
          ideal_distance_processed <- if(distance < distances_numeric[i, j]) distances_numeric[i, j] else NULL
        } else if(threshold_mask[i, j] == -1) {  # "<" case
          ideal_distance_processed <- if(distance > distances_numeric[i, j]) distances_numeric[i, j] else NULL
        } else {  # normal value
          ideal_distance_processed <- distances_numeric[i, j]
        }

        if (!is.null(ideal_distance_processed)) {
          # Spring force
          adjustment_factor <- 2*k*(ideal_distance_processed-distance)*delta/distance_01
          positions[i,] <- positions[i,] - adjustment_factor/(4*node_degrees_1[i]+k)
          positions[j,] <- positions[j,] + adjustment_factor/(4*node_degrees_1[j]+k)
        } else {
          # Repulsive force for >threshold measurements
          force <- c_repulsion/(2*distance_01^2)*(delta/distance_01)
          positions[i,] <- positions[i,] - force/node_degrees_1[i]
          positions[j,] <- positions[j,] + force/node_degrees_1[j]
        }
      } else {
        # Repulsive force for missing measurements
        force <- c_repulsion/(2*distance_01^2)*(delta/distance_01)
        positions[i,] <- positions[i,] - force/node_degrees_1[i]
        positions[j,] <- positions[j,] + force/node_degrees_1[j]
      }
    }

    # Check for numerical stability less frequently
    if(iter %% 5 == 0) {
      if(any(!is.finite(positions))) {
        warning("Numerical instability detected at iteration ", iter)
        stop <- TRUE
        break
      }
    }

    if(stop) break

    # Calculate current distances and convergence error
    p_dist_mat <- as.matrix(stats::dist(positions))
    rownames(p_dist_mat) <- rownames(positions)
    colnames(p_dist_mat) <- rownames(positions)

    dissimilarity_matrix_convergence_error <- vectorized_process_distance_matrix(distances_numeric, threshold_mask, p_dist_mat)

    # Vectorized convergence error calculation
    valid_mask <- !is.na(dissimilarity_matrix_convergence_error)
    convergence_error <- mean(abs(dissimilarity_matrix_convergence_error[valid_mask] - p_dist_mat[valid_mask]))

    # Calculate evaluation metrics similarly
    dissimilarity_matrix_raw <- suppressWarnings(as.numeric(dissimilarity_matrix))
    valid_indices_raw <- !is.na(dissimilarity_matrix_raw)
    mae <- mean(abs(dissimilarity_matrix_raw[valid_indices_raw] - p_dist_mat[valid_indices_raw]))

    if (verbose) {
      cat(sprintf(
        "Iteration=%d, MAE=%.4f, convergence_error=%.4f\n",
        iter, mae, convergence_error
      ))
    }

    # Check convergence
    if (iter > 1) {
      if (is.na(prev_convergence_error) ||
          is.na(convergence_error) ||
          convergence_error > convergence_error0 ||
          is.na((prev_convergence_error - convergence_error)/prev_convergence_error)) {

        if (verbose) {
          cat(paste(
            "! Skipping convergence check for this iteration.",
            "Please check model parameters k, cooling_rate, and c_repulsion.\n"
          ))
          cat(sprintf(
            "ndim=%d, k0=%.4f, cooling_rate=%.4f, c_repulsion=%.4f, MAE=%.4f, convergence_error=%.4f\n",
            ndim, k0, cooling_rate, c_repulsion, mae, convergence_error
          ))
        }
      } else {
        if ((prev_convergence_error - convergence_error)/prev_convergence_error <
            relative_epsilon) {
          convergence_counter <- convergence_counter - 1
          if(convergence_counter <= 0) {
            if(verbose) cat("Convergence achieved\n")
            break
          }
        }
        prev_convergence_error <- convergence_error
      }
    }

    # Update spring constant
    k <- k * (1 - cooling_rate)
  }

  # Save positions if requested
  if(write_positions_to_csv) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }
    if (missing(output_dir)) {
      stop("An 'output_dir' must be provided when 'write_positions_to_csv' is TRUE.", call. = FALSE)
    }
    csv_filename <- sprintf(
      "Positions_dim_%d_k0_%.4f_cooling_%.4f_c_repulsion_%.4f.csv",
      ndim, k0, cooling_rate, c_repulsion
    )
    # Write to the SPECIFIED directory
    full_path <- file.path(output_dir, csv_filename)
    utils::write.csv(positions, file = full_path, row.names = TRUE)
  }

  # Create result object
  result <- structure(
    list(
      positions = positions,
      est_distances = p_dist_mat,
      mae = mae,
      iter = iter,
      parameters = list(
        ndim = ndim,
        k0 = k0,
        cooling_rate = cooling_rate,
        c_repulsion = c_repulsion
      ),
      convergence = list(
        achieved = convergence_counter <= 0,
        error = convergence_error
      )
    ),
    class = "topolow"
  )

  return(result)
}



#' Main TopoLow algorithm implementation (DEPRECATED)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `create_topolow_map()` was deprecated in version 2.0.0 and will be removed in 
#' a future version. Please use [euclidean_embedding()] instead, which provides
#' the same functionality with improved performance and additional features.
#'
#' @details
#' This function has been superseded by [euclidean_embedding()], which offers:
#' * Enhanced matrix reordering for better optimization
#' * Improved parameter validation with informative warnings
#' * Consistent naming convention (dissimilarity vs distance)
#' * Better documentation and examples
#'
#' The core algorithm remains identical, ensuring your results will be equivalent.
#' The main changes are:
#' * Parameter name: `distance_matrix` --> `dissimilarity_matrix`
#' * Function name: `create_topolow_map()` --> `euclidean_embedding()`
#'
#' @param distance_matrix Matrix. Square, symmetric distance matrix. Can contain NA values
#'        for missing measurements and character strings with < or > prefixes for thresholded 
#'        measurements.
#' @param ndim Integer. Number of dimensions for the embedding space.
#' @param mapping_max_iter Integer. Maximum number of map optimization iterations.
#' @param k0 Numeric. Initial spring constant controlling spring forces.
#' @param cooling_rate Numeric. Rate of spring constant decay per iteration (0 < cooling_rate < 1).
#' @param c_repulsion Numeric. Repulsion constant controlling repulsive forces.
#' @param relative_epsilon Numeric. Convergence threshold for relative change in error.
#'        Default is 1e-4.
#' @param convergence_counter Integer. Number of iterations below threshold before declaring
#'        convergence. Default is 5.
#' @param initial_positions Matrix or NULL. Optional starting coordinates. If NULL,
#'        random initialization is used. Matrix should have nrow = nrow(distance_matrix)
#'        and ncol = ndim.
#' @param write_positions_to_csv Logical. Whether to save point positions to CSV file.
#'        Default is FALSE.
#' @param output_dir Character. Directory to save CSV file. Required if 
#'        `write_positions_to_csv` is TRUE.
#' @param verbose Logical. Whether to print progress messages. Default is FALSE.
#'
#' @return A `list` object of class `topolow`. This list contains the results of the
#'   optimization and includes the following components:
#' \itemize{
#'   \item `positions`: A `matrix` of the optimized point coordinates in the n-dimensional space.
#'   \item `est_distances`: A `matrix` of the Euclidean distances between points in the final optimized configuration.
#'   \item `mae`: The final Mean Absolute Error between the target distances and the estimated distances.
#'   \item `iter`: The total number of iterations performed before the algorithm terminated.
#'   \item `parameters`: A `list` containing the input parameters used for the optimization run.
#'   \item `convergence`: A `list` containing the final convergence status, including a logical `achieved` flag and the final `error` value.
#' }
#'
#' @examples
#' \donttest{
#' # Simple example (deprecated - use euclidean_embedding() instead)
#' dist_mat <- matrix(c(0, 2, 3, 2, 0, 4, 3, 4, 0), nrow=3)
#' 
#' # This will generate a deprecation warning
#' result <- create_topolow_map(
#'   dist_mat, 
#'   ndim = 2, 
#'   mapping_max_iter = 100,
#'   k0 = 1.0, 
#'   cooling_rate = 0.001, 
#'   c_repulsion = 0.01, 
#'   verbose = FALSE
#' )
#' 
#' # Recommended approach with new function:
#' result_new <- euclidean_embedding(
#'   dissimilarity_matrix = dist_mat,
#'   ndim = 2,
#'   mapping_max_iter = 100,
#'   k0 = 1.0,
#'   cooling_rate = 0.001,
#'   c_repulsion = 0.01,
#'   verbose = FALSE
#' )
#' }
#'
#' @seealso [euclidean_embedding()] for the replacement function.
#'
#' @keywords internal
#' @export
create_topolow_map <- function(distance_matrix,
                              ndim,
                              mapping_max_iter=1000,
                              k0,
                              cooling_rate,
                              c_repulsion,
                              relative_epsilon = 1e-4,
                              convergence_counter = 3,
                              initial_positions = NULL,
                              write_positions_to_csv = FALSE,
                              output_dir,
                              verbose = FALSE) {
  
  # Issue deprecation warning
  lifecycle::deprecate_warn(
    when = "2.0.0",
    what = "create_topolow_map()",
    with = "euclidean_embedding()",
    details = c(
      "i" = "The new function provides the same functionality with improvements:",
      "*" = "Parameter name: 'distance_matrix' --> 'dissimilarity_matrix'",
      "*" = "Enhanced matrix reordering for better optimization"
    )
  )
  
  # Handle missing output_dir parameter for backward compatibility
  if (write_positions_to_csv && missing(output_dir)) {
    output_dir <- getwd()
    warning("output_dir not specified, using current working directory")
  }
  
  # Call the new function with parameter mapping
  result <- euclidean_embedding(
    dissimilarity_matrix = distance_matrix,
    ndim = ndim,
    mapping_max_iter = mapping_max_iter,
    k0 = k0,
    cooling_rate = cooling_rate,
    c_repulsion = c_repulsion,
    relative_epsilon = relative_epsilon,
    convergence_counter = convergence_counter,
    initial_positions = initial_positions,
    write_positions_to_csv = write_positions_to_csv,
    output_dir = if(write_positions_to_csv) output_dir else NULL,
    verbose = verbose
  )
  return(result)
}


#' Print method for topolow objects
#'
#' Provides a concise display of key optimization results from `euclidean_embedding`.
#'
#' @param x A `topolow` object returned by `euclidean_embedding()`.
#' @param ... Additional arguments passed to print (not used).
#' @return The original `topolow` object (invisibly). This function is called for its
#' side effect of printing a summary to the console.
#' @examples
#' # Create a simple dissimilarity matrix and run the optimization
#' dist_mat <- matrix(c(0, 2, 3, 2, 0, 4, 3, 4, 0), nrow=3)
#' result <- euclidean_embedding(dist_mat, ndim=2, mapping_max_iter=50,
#'                            k0=1.0, cooling_rate=0.001, c_repulsion=0.1,
#'                            verbose = FALSE)
#' # Print the result object
#' print(result)
#' @export
print.topolow <- function(x, ...) {
  cat("topolow optimization result:\n")
  cat(sprintf("Dimensions: %d\n", x$parameters$ndim))
  cat(sprintf("Iterations: %d\n", x$iter))
  cat(sprintf("MAE: %.4f\n", x$mae))
  cat(sprintf("Convergence achieved: %s\n", x$convergence$achieved))
  cat(sprintf("Final convergence error: %.4f\n", x$convergence$error))
  invisible(x)
}


#' Summary method for topolow objects
#'
#' Provides a more detailed summary of the optimization results from `euclidean_embedding`,
#' including parameters, convergence, and performance metrics.
#'
#' @param object A `topolow` object returned by `euclidean_embedding()`.
#' @param ... Additional arguments passed to summary (not used).
#' @return No return value. This function is called for its side effect of
#' printing a detailed summary to the console.
#' @examples
#' # Create a simple dissimilarity matrix and run the optimization
#' dist_mat <- matrix(c(0, 2, 3, 2, 0, 4, 3, 4, 0), nrow=3)
#' result <- euclidean_embedding(dist_mat, ndim=2, mapping_max_iter=50,
#'                            k0=1.0, cooling_rate=0.001, c_repulsion=0.1,
#'                            verbose = FALSE)
#' # Summarize the result object
#' summary(result)
#' @export
summary.topolow <- function(object, ...) {
  print(object)
  cat("\nParameters:\n")
  cat(sprintf("k0: %.4f\n", object$parameters$k0))
  cat(sprintf("cooling_rate: %.4f\n", object$parameters$cooling_rate))
  cat(sprintf("c_repulsion: %.4f\n", object$parameters$c_repulsion))
}



#' Automatic Euclidean Embedding with Parameter Optimization
#'
#' @description
#' A user-friendly wrapper function that automatically optimizes parameters and 
#' performs Euclidean embedding on a dissimilarity matrix. This function handles
#' the entire workflow from parameter optimization to final embedding.
#'
#' @param dissimilarity_matrix Square symmetric dissimilarity matrix. Can contain
#'        NA values for missing measurements and threshold indicators (< or >).
#' @param output_dir Character. Directory for saving optimization files and results.
#'        Required - no default.
#' @param ndim_range Integer vector of length 2. Range for number of dimensions 
#'        (minimum, maximum). Default: c(2, 10)
#' @param k0_range Numeric vector of length 2. Range for initial spring constant
#'        (minimum, maximum). Default: c(0.1, 15)
#' @param cooling_rate_range Numeric vector of length 2. Range for cooling rate
#'        (minimum, maximum). Default: c(0.001, 0.07)
#' @param c_repulsion_range Numeric vector of length 2. Range for repulsion constant
#'        (minimum, maximum). Default: c(0.001, 0.4)
#' @param n_initial_samples Integer. Number of samples for initial parameter 
#'        optimization. Default: 100
#' @param n_adaptive_samples Integer. Number of samples for adaptive refinement.
#'        Default: 250  
#' @param max_cores Integer. Maximum number of cores to use. Default: NULL (auto-detect)
#' @param folds Integer. Number of cross-validation folds. Default: 20
#' @param mapping_max_iter Integer. Maximum iterations for final embedding. 
#'        Half this value is used for parameter search. Default: 1000
#' @param clean_intermediate Logical. Whether to remove intermediate files. Default: TRUE
#' @param verbose Character. Verbosity level: "off" (no output), "standard" (progress updates),
#'        or "full" (detailed output including from internal functions). Default: "standard"
#' @param fallback_to_defaults Logical. Whether to use default parameters if 
#'        optimization fails. Default: TRUE
#' @param save_results Logical. Whether to save the final positions as CSV. Default: FALSE
#'
#' @return A list containing:
#'   \item{positions}{Matrix of optimized coordinates}
#'   \item{est_distances}{Matrix of estimated distances}
#'   \item{mae}{Mean absolute error}
#'   \item{optimal_params}{List of optimal parameters found, including cross-validation MAE during optimization}
#'   \item{optimization_summary}{Summary of the optimization process}
#'   \item{data_characteristics}{Summary of input data characteristics}
#'   \item{runtime}{Total runtime in seconds}
#'
#' @examples
#' # Example 1: Basic usage with small matrix
#' test_data <- data.frame(
#' object = rep(paste0("Obj", 1:4), each = 4),
#' reference = rep(paste0("Ref", 1:4), 4),
#' score = sample(c(1, 2, 4, 8, 16, 32, 64, "<1", ">12"), 16, replace = TRUE)
#' )
#' dist_mat <- list_to_matrix(
#'   data = test_data,  # Pass the data frame, not file path
#'   object_col = "object",
#'   reference_col = "reference",
#'   value_col = "score",
#'   is_similarity = TRUE
#' )
#' \dontrun{
#' # Note: output_dir is required for actual use
#' result <- Euclidify(
#'   dissimilarity_matrix = dist_mat,
#'   output_dir = tempdir()  # Use temp directory for example
#' )
#' coordinates <- result$positions
#' }
#' 
#' # Example 2: Using custom parameter ranges
#' \dontrun{
#' result <- Euclidify(
#'   dissimilarity_matrix = dist_mat,
#'   output_dir = tempdir(),
#'   n_initial_samples = 10,
#'   n_adaptive_samples = 7,
#'   verbose = "off"
#' )
#' }
#' 
#' # Example 3: Handling missing data
#' dist_mat_missing <- dist_mat
#' dist_mat_missing[1, 3] <- dist_mat_missing[3, 1] <- NA
#' \dontrun{
#' result <- Euclidify(
#'   dissimilarity_matrix = dist_mat_missing,
#'   output_dir = tempdir(),
#'   n_initial_samples = 10,
#'   n_adaptive_samples = 7,
#'   verbose = "off"
#' )
#' }
#' 
#' # Example 4: Using threshold indicators
#' dist_mat_threshold <- dist_mat
#' dist_mat_threshold[1, 2] <- ">2"
#' dist_mat_threshold[2, 1] <- ">2"
#' \dontrun{
#' result <- Euclidify(
#'   dissimilarity_matrix = dist_mat_threshold,
#'   output_dir = tempdir(),
#'   n_initial_samples = 10,
#'   n_adaptive_samples = 7,
#'   verbose = "off"
#' )
#' }
#' 
#' # Example 5: Parallel processing with custom cores
#' \dontrun{
#' result <- Euclidify(
#'   dissimilarity_matrix = dist_mat,
#'   output_dir = tempdir(),
#'   max_cores = 4,
#'   n_adaptive_samples = 100,
#'   save_results = TRUE  # Save positions to CSV
#' )
#' }
#'
#' @export
Euclidify <- function(dissimilarity_matrix,
                      output_dir,
                      ndim_range = c(2, 10),
                      k0_range = c(0.1, 20),
                      cooling_rate_range = c(0.0001, 0.1),
                      c_repulsion_range = c(0.0001, 1),
                      n_initial_samples = 50,
                      n_adaptive_samples = 150,
                      max_cores = NULL,
                      folds = 20,
                      mapping_max_iter = 500,
                      clean_intermediate = TRUE,
                      verbose = "standard",
                      fallback_to_defaults = FALSE,
                      save_results = FALSE) {
  
  # Start timing
  start_time <- Sys.time()
  
  # Input validation
  if (missing(output_dir)) {
    stop("output_dir must be specified.")
  }
  
  if (!is.matrix(dissimilarity_matrix)) {
    stop("dissimilarity_matrix must be a matrix")
  }
  if (nrow(dissimilarity_matrix) != ncol(dissimilarity_matrix)) {
    stop("dissimilarity_matrix must be square")
  }
  if (length(ndim_range) != 2 || ndim_range[1] > ndim_range[2]) {
    stop("ndim_range must be a vector of length 2 with min <= max")
  }
  if (length(k0_range) != 2 || k0_range[1] > k0_range[2]) {
    stop("k0_range must be a vector of length 2 with min <= max")
  }
  if (length(cooling_rate_range) != 2 || cooling_rate_range[1] < 0 || cooling_rate_range[2] > 1 ||
      cooling_rate_range[1] > cooling_rate_range[2]) {
    stop("cooling_rate_range must be a vector of length 2 with 0 <= min <= max <= 1")
  }
  if (length(c_repulsion_range) != 2 || c_repulsion_range[1] < 0 || c_repulsion_range[2] < c_repulsion_range[1]) {
    stop("c_repulsion_range must be a vector of length 2 with min >= 0 and min <= max")
  }
  if (!verbose %in% c("off", "standard", "full")) {
    stop("verbose must be 'off', 'standard', or 'full'")
  }
  
  # Set verbose flags for internal functions
  verbose_internal <- (verbose == "full")
  verbose_main <- (verbose %in% c("standard", "full"))
  
  # Calculate mapping_max_iter for search (half of final)
  mapping_max_iter_search <- round(mapping_max_iter / 2)
  
  # Determine cores
  if (is.null(max_cores)) {
    available_cores <- future::availableCores()
    max_cores <- max(1, available_cores - 1)
  }
  
  if (verbose_main) {
    cat("=== TOPOLOW AUTOMATIC EUCLIDEAN EMBEDDING ===\n")
    cat("Dataset size:", nrow(dissimilarity_matrix), "x", ncol(dissimilarity_matrix), "\n")
    cat("Missing values:", sum(is.na(dissimilarity_matrix)), 
        "(", round(sum(is.na(dissimilarity_matrix)) / (nrow(dissimilarity_matrix)^2) * 100, 1), "%)\n")
    cat("Using", max_cores, "cores\n")
    cat("Output directory:", output_dir, "\n\n")
  }
  
  # Create output directory structure
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  optimization_dir <- file.path(output_dir, "optimization")
  if (dir.exists(optimization_dir) && clean_intermediate) {
    unlink(optimization_dir, recursive = TRUE)
  }
  if (!dir.exists(optimization_dir)) {
    dir.create(optimization_dir, recursive = TRUE)
  }
  
  # Initialize results storage
  optimization_summary <- list()
  optimal_params <- NULL
  
  # Step 1: Assess data characteristics
  if (verbose_main) cat("Step 1: Assessing data characteristics...\n")
  
  # Calculate non-Euclidean character
  n <- nrow(dissimilarity_matrix)
  # Convert to numeric matrix properly (fixing the bug)
  dissim_numeric <- matrix(suppressWarnings(as.numeric(as.vector(dissimilarity_matrix))), 
                           nrow = n, ncol = n)
  dissim_numeric[is.na(dissim_numeric)] <- median(dissim_numeric, na.rm = TRUE)
  
  D_squared <- dissim_numeric^2
  J <- diag(n) - (1/n) * matrix(1, n, n)
  B <- -0.5 * J %*% D_squared %*% J
  
  eigenvals <- eigen(B, symmetric = TRUE)$values
  positive_eigenvals <- eigenvals[eigenvals > 1e-12]
  negative_eigenvals <- eigenvals[eigenvals < -1e-12]
  
  if (length(positive_eigenvals) > 0) {
    deviation_score <- sum(abs(negative_eigenvals)) / sum(positive_eigenvals)
  } else {
    deviation_score <- 1.0
  }
  
  # Determine optimal dimensions based on eigenvalues
  cumulative_variance <- cumsum(positive_eigenvals) / sum(positive_eigenvals)
  dims_90_percent <- which(cumulative_variance >= 0.90)[1]
  if (is.na(dims_90_percent)) dims_90_percent <- length(positive_eigenvals)
  dims_90_percent <- max(1, dims_90_percent)  # Ensure at least 1 dimension
  
  if (verbose_main) {
    cat("  Non-Euclidean deviation score:", round(deviation_score, 4), "\n")
    cat("    (Higher values indicate greater deviation from Euclidean geometry.\n")
    cat("     Values > 0.1 suggest significant non-Euclidean structure.)\n")
    cat("  Dimensions for 90% variance:", dims_90_percent, "\n")
  }
  
  # Adjust ndim_range based on data characteristics (improved logic)
  original_max <- ndim_range[2]
  suggested_max <- dims_90_percent + 2
  
  if (suggested_max < ndim_range[1]) {
    if (verbose_main) {
      cat("  WARNING: Suggested max dimensions (", suggested_max, 
          ") is less than your minimum (", ndim_range[1], ").\n", sep = "")
      cat("  Consider lowering ndim_range minimum for better results.\n")
    }
  } else if (suggested_max < ndim_range[2]) {
    ndim_range[2] <- suggested_max
    if (verbose_main) {
      cat("  Adjusted max dimensions from", original_max, "to", ndim_range[2], 
          "based on data characteristics.\n")
    }
  }
  
  # Step 2: Initial parameter optimization (returns log-transformed parameters)
  if (verbose_main) cat("\nStep 2: Initial parameter optimization...\n")
  
  tryCatch({
    initial_results <- initial_parameter_optimization(
      dissimilarity_matrix = dissimilarity_matrix,
      mapping_max_iter = mapping_max_iter_search,
      relative_epsilon = 1e-4,
      convergence_counter = 3,
      scenario_name = "auto_optimization",
      N_min = ndim_range[1],
      N_max = ndim_range[2],
      k0_min = k0_range[1],
      k0_max = k0_range[2],
      cooling_rate_min = cooling_rate_range[1],
      cooling_rate_max = cooling_rate_range[2],
      c_repulsion_min = c_repulsion_range[1],
      c_repulsion_max = c_repulsion_range[2],
      num_samples = n_initial_samples,
      max_cores = max_cores,
      folds = folds,
      verbose = verbose_internal,
      write_files = TRUE,
      output_dir = optimization_dir
    )
    
    optimization_summary$initial_results <- initial_results
    
    if (nrow(initial_results) > 0) {
      if (verbose_main) {
        cat("  Initial optimization completed successfully\n")
        cat("  Best initial MAE:", round(min(initial_results$Holdout_MAE, na.rm = TRUE), 4), "\n")
        cat("  Parameters are in log scale for adaptive sampling\n")
      }
      
      # Step 3: Adaptive sampling refinement (parameters are already log-transformed)
      samples_file <- file.path(optimization_dir, "model_parameters", 
                                "auto_optimization_model_parameters.csv")
      
      if (file.exists(samples_file)) {
        if (verbose_main) cat("\nStep 3: Adaptive parameter refinement...\n")

        tryCatch({
          run_adaptive_sampling(
            initial_samples_file = samples_file,
            scenario_name = "auto_optimization",
            dissimilarity_matrix = dissimilarity_matrix,
            max_cores = max_cores,
            num_samples = n_adaptive_samples,
            mapping_max_iter = mapping_max_iter,
            relative_epsilon = 1e-5,
            folds = folds,
            output_dir = optimization_dir,
            verbose = verbose_internal
          )
        }, error = function(e) {
          if (verbose_main) cat("  WARNING: Adaptive sampling failed:", e$message, "\n")
          if (verbose_main) cat("  Continuing with initial optimization results\n")
          optimization_summary$adaptive_error <- e$message
        })

        # Step 4: Extract optimal parameters (improved handling)
        final_params_file <- file.path(optimization_dir, "model_parameters",
                                        "auto_optimization_model_parameters.csv")

        if (file.exists(final_params_file)) {
          final_params <- tryCatch({
            results <- read.csv(final_params_file, stringsAsFactors = FALSE)
            
            # More robust cleaning
            results <- results %>%
              filter(!is.na(.data$NLL) & !is.na(.data$Holdout_MAE) & 
                      is.finite(.data$NLL) & is.finite(.data$Holdout_MAE)) %>%
              na.omit()
            
            if (nrow(results) > 5) {  # Only apply outlier cleaning if we have enough data
              results <- as.data.frame(lapply(results, clean_data, k = 3))
              results <- na.omit(results)
            }
            
            results
          }, error = function(e) {
            if (verbose_main) cat("  WARNING: Error reading final results:", e$message, "\n")
            data.frame()  # Return empty data frame
          })
          
          if (nrow(final_params) > 0) {
            # Find best parameters (they are already in log scale)
            best_idx <- which.min(final_params$Holdout_MAE)
            
            optimal_params <- list(
              ndim = round(exp(final_params$log_N[best_idx])),
              k0 = exp(final_params$log_k0[best_idx]),
              cooling_rate = exp(final_params$log_cooling_rate[best_idx]),
              c_repulsion = exp(final_params$log_c_repulsion[best_idx]),
              CV_MAE = final_params$Holdout_MAE[best_idx],
              nll = final_params$NLL[best_idx]
            )
            
            # Validate optimal parameters are reasonable
            if (optimal_params$ndim < 1 || optimal_params$ndim > 50 ||
                optimal_params$k0 <= 0 || optimal_params$k0 > 100 ||
                optimal_params$cooling_rate <= 0 || optimal_params$cooling_rate >= 1 ||
                optimal_params$c_repulsion <= 0 || optimal_params$c_repulsion > 10) {
              
              if (verbose_main) cat("  WARNING: Optimal parameters are unreasonable, using fallback\n")
              optimal_params <- NULL
            } else {
              if (verbose_main) {
                cat("  Adaptive refinement completed successfully\n")
                cat("  Optimal parameters found:\n")
                cat("    - Dimensions:", optimal_params$ndim, "\n")
                cat("    - k0:", round(optimal_params$k0, 4), "\n")
                cat("    - cooling_rate:", round(optimal_params$cooling_rate, 6), "\n")
                cat("    - c_repulsion:", round(optimal_params$c_repulsion, 6), "\n")
                cat("    - Cross-validation MAE:", round(optimal_params$CV_MAE, 4), "\n")
                cat("  We recommend saving these parameters to set better ranges and \n")
                cat("  reduced iterations in future runs. It saves time.\n")

                # Issue warning if parameters are too close (0.05 of the corresponding range) to their range limits 
                if (abs(optimal_params$ndim - ndim_range[1]) < 0.05 * (ndim_range[2] - ndim_range[1]) ||
                    abs(optimal_params$k0 - k0_range[1]) < 0.05 * (k0_range[2] - k0_range[1]) ||
                    abs(optimal_params$cooling_rate - cooling_rate_range[1]) < 0.05 * (cooling_rate_range[2] - cooling_rate_range[1]) ||
                    abs(optimal_params$c_repulsion - c_repulsion_range[1]) < 0.05 * (c_repulsion_range[2] - c_repulsion_range[1]) ||
                    abs(optimal_params$ndim - ndim_range[2]) < 0.05 * (ndim_range[2] - ndim_range[1]) ||
                    abs(optimal_params$k0 - k0_range[2]) < 0.05 * (k0_range[2] - k0_range[1]) ||
                    abs(optimal_params$cooling_rate - cooling_rate_range[2]) < 0.05 * (cooling_rate_range[2] - cooling_rate_range[1]) ||
                    abs(optimal_params$c_repulsion - c_repulsion_range[2]) < 0.05 * (c_repulsion_range[2] - c_repulsion_range[1])) {
                  cat("  WARNING: Some optimal parameters are very close to their range limits.\n")
                }
                cat("  This may indicate the need for wider parameter ranges in future runs.\n")
              }
              
              optimization_summary$final_params <- final_params
            }
          }
        }
      }
    }
    
  }, error = function(e) {
    if (verbose_main) cat("  WARNING: Parameter optimization failed:", e$message, "\n")
    optimization_summary$error <- e$message
  })
  
  # Fallback to initial results if adaptive failed
  if (is.null(optimal_params) && !is.null(optimization_summary$initial_results)) {
    initial_results <- optimization_summary$initial_results
    if (nrow(initial_results) > 0) {
      # Initial results are now in log scale, so we need to extract them properly
      valid_results <- initial_results %>%
        filter(!is.na(.data$Holdout_MAE) & !is.na(.data$log_N) & !is.na(.data$log_k0) & 
              !is.na(.data$log_cooling_rate) & !is.na(.data$log_c_repulsion) &
              is.finite(.data$Holdout_MAE))
      
      if (nrow(valid_results) > 0) {
        best_idx <- which.min(valid_results[["Holdout_MAE"]])
        optimal_params <- list(
          ndim = round(exp(valid_results$log_N[best_idx])),
          k0 = exp(valid_results$log_k0[best_idx]),
          cooling_rate = exp(valid_results$log_cooling_rate[best_idx]),
          c_repulsion = exp(valid_results$log_c_repulsion[best_idx]),
          CV_MAE = valid_results$Holdout_MAE[best_idx]
        )
        if (verbose_main) cat("  Using initial optimization results\n")
      }
    }
  }

  # Warn if no optimal parameters found and stop the operation
  if (is.null(optimal_params) && !fallback_to_defaults) {
    if (verbose_main) {
      cat("  WARNING: No optimal parameters found from initial or adaptive optimization.\n")
      cat("  Consider adjusting parameter ranges or if you just want to get a less\n")
      cat("  accurate map with non-optimized parameters, set fallback_to_defaults = TRUE.\n")
    }
    stop(" STOP: Parameter optimization failed and fallback is disabled")
  }
  # Enhanced fallback to default parameters with better data-driven selection
  if (is.null(optimal_params) && fallback_to_defaults) {
    # Respect user's ndim_range while considering data characteristics
    fallback_ndim <- max(ndim_range[1], min(ndim_range[2], max(2, dims_90_percent)))
    
    # Scale default parameters based on matrix size and characteristics
    matrix_size <- nrow(dissimilarity_matrix)
    sparsity <- sum(is.na(dissimilarity_matrix)) / (matrix_size^2)
    
    # Adjust parameters based on data characteristics
    base_k0 <- if (sparsity > 0.5) 8.0 else 6.0  # Higher k0 for sparse data
    base_cooling <- if (sparsity > 0.5) 0.005 else 0.01  # Slower cooling for sparse data
    base_repulsion <- if (deviation_score > 0.2) 0.01 else 0.02  # Higher repulsion for non-Euclidean data
    
    optimal_params <- list(
      ndim = fallback_ndim,
      k0 = base_k0,
      cooling_rate = base_cooling,
      c_repulsion = base_repulsion
    )
    
    if (verbose_main) {
      cat("\n  Using data-driven default parameters:\n")
      cat("    - Dimensions:", optimal_params$ndim, "(based on eigenvalue analysis)\n")
      cat("    - k0:", optimal_params$k0, "(adjusted for sparsity:", round(sparsity*100, 1), "%)\n")
      cat("    - cooling_rate:", optimal_params$cooling_rate, "\n")
      cat("    - c_repulsion:", optimal_params$c_repulsion, "(adjusted for non-Euclidean score:", round(deviation_score, 3), ")\n")
    }
    optimization_summary$used_defaults <- TRUE
    optimization_summary$default_reasoning <- list(
      sparsity = sparsity,
      deviation_score = deviation_score,
      dims_90_percent = dims_90_percent
    )
  }
  
  
  # Step 4: Final embedding with optimal parameters
  if (verbose_main) cat("\nStep 4: Running final embedding with optimal parameters...\n")
  
  embedding_result <- tryCatch({
    euclidean_embedding(
      dissimilarity_matrix = dissimilarity_matrix,
      ndim = optimal_params$ndim,
      mapping_max_iter = mapping_max_iter * 2,  # Extra iterations for final
      k0 = optimal_params$k0,
      cooling_rate = optimal_params$cooling_rate,
      c_repulsion = optimal_params$c_repulsion,
      relative_epsilon = 1e-6,
      convergence_counter = 5,
      verbose = verbose_internal,
      write_positions_to_csv = save_results,
      output_dir = output_dir
    )
  }, error = function(e) {
    if (verbose_main) cat("  ERROR: Final embedding failed:", e$message, "\n")
    NULL
  })
  
  if (is.null(embedding_result)) {
    stop("Final embedding failed")
  }
  
  if (verbose_main) {
    cat("  Embedding completed successfully\n")
    cat("  Final MAE:", round(embedding_result$mae, 4), "\n")
    cat("  Iterations:", embedding_result$iter, "\n")
    cat("  Convergence achieved:", embedding_result$convergence$achieved, "\n")
  }
  
  # Clean up intermediate files
  if (clean_intermediate) {
    unlink(optimization_dir, recursive = TRUE)
    if (verbose_main) cat("\nCleaned up intermediate files\n")
  }
  
  # Prepare final results (without redundant embedding_result)
  results <- list(
    positions = embedding_result$positions,
    est_distances = embedding_result$est_distances,
    mae = embedding_result$mae,
    optimal_params = optimal_params,
    optimization_summary = optimization_summary,
    data_characteristics = list(
      n_objects = nrow(dissimilarity_matrix),
      missing_percentage = sum(is.na(dissimilarity_matrix)) / (nrow(dissimilarity_matrix)^2) * 100,
      non_euclidean_score = deviation_score,
      eigenvalue_dims_90 = dims_90_percent
    ),
    runtime = difftime(Sys.time(), start_time, units = "secs")
  )
  
  # Save positions as CSV if requested
  if (save_results) {
    positions_file <- file.path(output_dir, "euclidify_positions.csv")
    write.csv(embedding_result$positions, positions_file, row.names = TRUE)
    if (verbose_main) cat("\nPositions saved to:", positions_file, "\n")
  }
  
  if (verbose_main) {
    cat("\n=== EUCLIDIFY COMPLETED ===\n")
    cat("Total runtime:", round(as.numeric(results$runtime), 1), "seconds\n")
    cat("Final embedding dimensions:", optimal_params$ndim, "\n")
    cat("Final MAE:", round(embedding_result$mae, 4), "\n")
  }
  
  return(results)
}
