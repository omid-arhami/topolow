# Copyright (c) 2024 Omid Arhami omid.arhami@uga.edu
# R/core.R

#' @useDynLib topolow, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

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
#' This version uses an optimized C++ backend with:
#' * Compressed edge list (COO format) for memory efficiency
#' * Negative sampling to approximate unmeasured pair repulsion
#' * Immediate (Gauss-Seidel) position updates for better convergence
#' * Stochastic edge shuffling for escaping local optima
#'
#' @details
#' The algorithm iteratively updates point positions using:
#' * Spring forces between points with measured dissimilarities.
#' * Repulsive forces between points without measurements (via negative sampling).
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
#' @param convergence_counter Integer. Number of consecutive iterations below threshold before 
#'        declaring convergence. Default is 5.
#' @param initial_positions Matrix or NULL. Optional starting coordinates. If NULL,
#'        random initialization is used. Matrix should have nrow = nrow(dissimilarity_matrix)
#'        and ncol = ndim.
#' @param write_positions_to_csv Logical. Whether to save point positions to a CSV file.
#'        Default is FALSE.
#' @param output_dir Character. Directory to save the CSV file. Required if
#'        `write_positions_to_csv` is TRUE.
#' @param verbose Logical. Whether to print progress messages. Default is FALSE.
#' @param n_negative_samples Integer. Number of negative samples per edge endpoint. 
#'        Higher values better approximate the original O(N^2) algorithm but increase 
#'        computation time. Default is 5.
#' @param convergence_check_freq Integer. How often to check for convergence (every N iterations).
#'        Lower values give more precise stopping but add overhead. Default is 10.
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
                                mapping_max_iter = 1000,
                                k0,
                                cooling_rate,
                                c_repulsion,
                                relative_epsilon = 1e-4,
                                convergence_counter = 5,
                                initial_positions = NULL,
                                write_positions_to_csv = FALSE,
                                output_dir,
                                verbose = FALSE,
                                n_negative_samples = 5,
                                convergence_check_freq = 10) {

  # ===========================================================================
  # INPUT VALIDATION
  # ===========================================================================
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
  if (!is.numeric(n_negative_samples) || n_negative_samples < 0) {
    stop("n_negative_samples must be a non-negative integer")
  }
  if (!is.numeric(convergence_check_freq) || convergence_check_freq < 1) {
    stop("convergence_check_freq must be a positive integer")
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

  # ===========================================================================
  # MATRIX REORDERING : Reorder rows and columns so largest values are furthest from diagonal
  # ===========================================================================
  if (n > 1) {
    # Extract numeric values for clustering, handling threshold indicators
    # VECTORIZED: Replace nested for loops with matrix operations
    numeric_matrix <- matrix(NA_real_, n, n)
    
    # Get non-NA mask
    non_na_mask <- !is.na(dissimilarity_matrix)
    
    # Check if any values are character type
    if (is.character(dissimilarity_matrix)) {
      # For character matrices: remove < or > prefixes and convert to numeric
      # gsub operates element-wise on the entire matrix
      numeric_matrix[non_na_mask] <- as.numeric(gsub("^[<>]", "", dissimilarity_matrix[non_na_mask]))
    } else {
      # For numeric matrices: direct conversion
      numeric_matrix[non_na_mask] <- as.numeric(dissimilarity_matrix[non_na_mask])
    }

    # Create spectral ordering to concentrate largest values in corners
    tryCatch({
      # Calculate average dissimilarity for each point using only non-NA values
      # VECTORIZED: Replace for loop with matrix operations
      # Set diagonal to NA to exclude self-comparisons
      diag(numeric_matrix) <- NA
      
      # Compute row and column means, treating NA appropriately
      row_means <- rowMeans(numeric_matrix, na.rm = TRUE)
      col_means <- colMeans(numeric_matrix, na.rm = TRUE)
      
      # Average of row and column means (each point appears in both)
      avg_dissim <- (row_means + col_means) / 2
      
      # Handle points with no measurements (NaN from all-NA rows/cols)
      avg_dissim[is.nan(avg_dissim)] <- 0

      # Check if we have meaningful averages (not all zeros)
      if (sum(avg_dissim > 0) > 1) {
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

  # ===========================================================================
  # PREPROCESSING FOR C++ (using COO format instead of CSC to avoid zero-dropping)
  # ===========================================================================
  
  # Pre-compute NA status once
  is_na_matrix <- is.na(dissimilarity_matrix)
  node_degrees <- rowSums(!is_na_matrix)
  
  # Parse dissimilarity matrix to numeric values and threshold codes
  # VECTORIZED: Replace nested for loops with logical indexing
  distances_numeric <- matrix(Inf, n, n)
  threshold_mask <- matrix(0L, n, n)  # 0 = number, 1 = >, -1 = <
  
  # Get non-NA mask
  non_na_mask <- !is_na_matrix
  
  if (is.character(dissimilarity_matrix)) {
    # For character matrices, identify threshold types using vectorized string operations
    gt_mask <- non_na_mask & startsWith(dissimilarity_matrix, ">")
    lt_mask <- non_na_mask & startsWith(dissimilarity_matrix, "<")
    normal_mask <- non_na_mask & !gt_mask & !lt_mask
    
    # Set threshold codes
    threshold_mask[gt_mask] <- 1L
    threshold_mask[lt_mask] <- -1L
    
    # Extract numeric values (remove < or > prefix where applicable)
    if (any(gt_mask)) {
      distances_numeric[gt_mask] <- as.numeric(sub("^>", "", dissimilarity_matrix[gt_mask]))
    }
    if (any(lt_mask)) {
      distances_numeric[lt_mask] <- as.numeric(sub("^<", "", dissimilarity_matrix[lt_mask]))
    }
    if (any(normal_mask)) {
      distances_numeric[normal_mask] <- as.numeric(dissimilarity_matrix[normal_mask])
    }
  } else {
    # For numeric matrices, direct assignment
    distances_numeric[non_na_mask] <- dissimilarity_matrix[non_na_mask]
  }
  
  # ===========================================================================
  # BUILD COO EDGE LIST (UPPER TRIANGLE ONLY, PRE-FILTERED IN R)
  # ===========================================================================
  # This avoids the sparse matrix zero-dropping bug and reduces memory transfer
  # VECTORIZED: Replace nested for loops with which() and logical indexing
  
  # Create upper triangle mask
  upper_tri_mask <- upper.tri(distances_numeric)
  
  # Find valid edges: upper triangle AND not Inf
  valid_edges <- upper_tri_mask & (distances_numeric != Inf)
  
  # Get indices of valid edges
  edge_indices <- which(valid_edges, arr.ind = TRUE)
  
  if (nrow(edge_indices) == 0) {
    stop("No valid off-diagonal measurements found in dissimilarity matrix")
  }
  
  # Convert to 0-based indexing for C++
  edge_i <- as.integer(edge_indices[, 1] - 1L)
  edge_j <- as.integer(edge_indices[, 2] - 1L)
  
  # Extract corresponding distances and thresholds using linear indexing
  linear_idx <- edge_indices[, 1] + (edge_indices[, 2] - 1L) * n
  edge_dist <- distances_numeric[linear_idx]
  edge_thresh <- threshold_mask[linear_idx]
  
  # Create flattened has_measurement vector for O(1) lookup in C++
  # Stored row-major: index = i * n + j
  has_measurement <- !is_na_matrix
  has_measurement_flat <- as.integer(as.vector(t(has_measurement)))

  # ===========================================================================
  # INITIALIZATION OF POINT POSITIONS
  # ===========================================================================
  if (is.null(initial_positions)) {
    dissimilarity_matrix_numeric <- suppressWarnings(as.numeric(as.character(dissimilarity_matrix)))
    dissimilarity_matrix_numeric[is.na(dissimilarity_matrix_numeric)] <- NA
    init_step <- max(dissimilarity_matrix_numeric, na.rm = TRUE) / n
    # Random initialization
    random_steps <- matrix(stats::runif((n - 1) * ndim, 0, 2 * init_step), nrow = n - 1, ncol = ndim)
    # First row stays zero, subsequent rows are cumulative sums
    initial_positions <- rbind(matrix(0, 1, ndim), apply(random_steps, 2, cumsum))
  }

  rownames(initial_positions) <- rownames(dissimilarity_matrix)

  # ===========================================================================
  # CALL C++ OPTIMIZATION
  # ===========================================================================
  if (verbose) cat("Starting C++ optimization...\n")

  start_time <- Sys.time()
  
  cpp_result <- optimize_layout_cpp(
    initial_positions = initial_positions,
    edge_i = edge_i,
    edge_j = edge_j,
    edge_dist = edge_dist,
    edge_thresh = edge_thresh,
    has_measurement_flat = has_measurement_flat,
    degrees = as.integer(node_degrees),
    n_points = as.integer(n),
    n_iter = as.integer(mapping_max_iter),
    k0 = k0,
    cooling_rate = cooling_rate,
    c_repulsion = c_repulsion,
    n_negative_samples = as.integer(n_negative_samples),
    relative_epsilon = relative_epsilon,
    convergence_window = as.integer(convergence_counter),
    convergence_check_freq = as.integer(convergence_check_freq),
    verbose = verbose
  )
  
  end_time <- Sys.time()
  
  if (verbose) {
    cat(sprintf("Optimization finished in %.2f seconds.\n", 
                as.numeric(difftime(end_time, start_time, units = "secs"))))
  }

  # ===========================================================================
  # POST-PROCESSING
  # ===========================================================================
  positions <- cpp_result$positions
  rownames(positions) <- rownames(dissimilarity_matrix)

  # Calculate current distances and final metrics
  p_dist_mat <- as.matrix(stats::dist(positions))
  rownames(p_dist_mat) <- rownames(positions)
  colnames(p_dist_mat) <- rownames(positions)

  # Calculate evaluation metrics
  dissimilarity_matrix_raw <- suppressWarnings(as.numeric(dissimilarity_matrix))
  valid_indices_raw <- !is.na(dissimilarity_matrix_raw)
  mae <- mean(abs(dissimilarity_matrix_raw[valid_indices_raw] - p_dist_mat[valid_indices_raw]))

  # ===========================================================================
  # SAVE POSITIONS IF REQUESTED
  # ===========================================================================
  if (write_positions_to_csv) {
    if (missing(output_dir)) {
      stop("An 'output_dir' must be provided when 'write_positions_to_csv' is TRUE.", call. = FALSE)
    }
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }
    csv_filename <- sprintf(
      "Positions_dim_%d_k0_%.4f_cooling_%.4f_c_repulsion_%.4f.csv",
      ndim, k0, cooling_rate, c_repulsion
    )
    full_path <- file.path(output_dir, csv_filename)
    utils::write.csv(positions, file = full_path, row.names = TRUE)
    if (verbose) cat("Positions saved to:", full_path, "\n")
  }

  # ===========================================================================
  # CREATE RESULT OBJECT (original + new fields)
  # ===========================================================================
  result <- structure(
    list(
      positions = positions,
      est_distances = p_dist_mat,
      mae = mae,
      iter = cpp_result$iterations,
      parameters = list(
        ndim = ndim,
        k0 = k0,
        cooling_rate = cooling_rate,
        c_repulsion = c_repulsion,
        n_negative_samples = n_negative_samples,
        method = "cpp_coo_neg_sampling"
      ),
      convergence = list(
        achieved = cpp_result$converged,
        error = cpp_result$final_mae,
        final_k = cpp_result$final_k
      )
    ),
    class = "topolow"
  )

  return(result)
}



# Copyright (c) 2024 Omid Arhami omid.arhami@uga.edu
# R/incremental.R
# Incremental embedding functionality for topolow

#' Incremental Euclidean Embedding
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Adds new points to an existing Euclidean embedding while keeping the positions
#' of existing points fixed. This is useful for incrementally updating antigenic
#' maps as new measurements become available, without re-optimizing the entire map.
#'
#' @details
#' The goal of incremental embedding is to produce results equivalent to 
#' re-optimizing the entire map from scratch, while being computationally 
#' efficient by only optimizing new points.
#'
#' ## Algorithm Steps
#' 1. Takes fixed positions of existing points and new measurements
#' 2. Identifies truly new points (those not in existing positions)
#' 3. Constructs a merged dissimilarity matrix:
#'    - Fixed-to-fixed distances come from existing positions (NEVER from new_measurements)
#'    - New-to-fixed and new-to-new distances come from new_measurements
#' 4. Reorders the matrix for optimal convergence (spectral ordering)
#' 5. Initializes new point positions near their nearest fixed anchor
#' 6. Optimizes only new point positions via C++ backend
#' 7. Only processes edges involving at least one new point (efficiency)
#' 8. Applies extra repulsion to compensate for reduced negative sampling
#'
#' ## Key Design Decisions
#' 
#' **Decision 1: Fixed Point Mask**
#' Points are marked as fixed (1) or optimizable (0). Fixed points never move;
#' only new points are optimized based on new measurements.
#'
#' **Decision 2: Edge Filtering**
#' Only edges involving at least one new point are processed. For 1000 fixed + 
#' 10 new points, this means ~20,000 edges instead of ~500,000.
#'
#' **Decision 3: Force Compensation**
#' In standard optimization, spring forces move BOTH endpoints. When one 
#' endpoint is fixed, only one point moves, causing ~half the distance change.
#' To compensate, forces are scaled by factor (1 + norm_movable/norm_fixed).
#' This ensures equivalent convergence behavior to full optimization.
#'
#' **Decision 4: Extra Repulsion Pass**
#' Edge filtering reduces negative sampling opportunities from fixed points.
#' To compensate, explicit repulsion is applied from randomly sampled fixed 
#' points to each new point after each iteration.
#'
#' **Decision 5: Contradictory Data Prevention**
#' Fixed-fixed distances in the matrix come ONLY from existing positions.
#' Any fixed-fixed measurements in new_measurements are ignored to prevent
#' impossible geometric constraints.
#'
#' **Decision 6: Immediate (Gauss-Seidel) Updates**
#' Positions are updated immediately when visited, matching the original
#' topolow algorithm's behavior.
#'
#' The function leverages the same C++ optimization backend as [euclidean_embedding()]
#' but with modifications for fixed-point handling and force compensation.
#'
#' @param fixed_positions Matrix. Coordinates of existing points that should remain
#'        fixed during optimization. Row names must identify the points.
#' @param new_measurements Data frame in long format containing new measurements.
#'        Must have columns for object identifiers, reference identifiers, and
#'        dissimilarity values.
#' @param object_col Character. Name of the column containing object identifiers.
#'        Default is "object".
#' @param reference_col Character. Name of the column containing reference identifiers.
#'        Default is "reference".
#' @param value_col Character. Name of the column containing dissimilarity values.
#'        Can include threshold indicators (< or >). Default is "value".
#' @param mapping_max_iter Integer. Maximum number of optimization iterations.
#'        Default is 1000.
#' @param k0 Numeric. Initial spring constant controlling spring forces. Default is 1.0.
#' @param cooling_rate Numeric. Rate of spring constant decay per iteration
#'        (0 < cooling_rate < 1). Default is 0.01.
#' @param c_repulsion Numeric. Repulsion constant controlling repulsive forces.
#'        Default is 0.01.
#' @param relative_epsilon Numeric. Convergence threshold for relative change in error.
#'        Default is 1e-4.
#' @param convergence_counter Integer. Number of consecutive iterations below threshold
#'        before declaring convergence. Default is 5.
#' @param n_negative_samples Integer. Number of negative samples per edge endpoint.
#'        Default is 5.
#' @param convergence_check_freq Integer. How often to check for convergence.
#'        Default is 10.
#' @param verbose Logical. Whether to print progress messages. Default is FALSE.
#'
#' @return A `list` object of class `topolow` containing:
#' \itemize{
#'   \item `positions`: Matrix of all point coordinates (fixed + new points)
#'   \item `est_distances`: Matrix of Euclidean distances in the final configuration
#'   \item `mae`: Mean Absolute Error on edges involving new points
#'   \item `iter`: Number of iterations performed
#'   \item `parameters`: List of optimization parameters used
#'   \item `convergence`: List with convergence status and final error
#'   \item `incremental_info`: List with details about the incremental update:
#'     \itemize{
#'       \item `n_fixed_points`: Number of points kept fixed
#'       \item `n_new_points`: Number of newly added points
#'       \item `n_new_edges`: Number of edges processed (involving new points)
#'       \item `new_point_names`: Names of the newly added points
#'       \item `warnings`: Any warnings generated during processing
#'     }
#' }
#'
#' @examples
#' # Create an initial map
#' initial_dist <- matrix(c(0, 2, 3, 2, 0, 4, 3, 4, 0), nrow = 3)
#' rownames(initial_dist) <- colnames(initial_dist) <- c("A", "B", "C")
#'
#' initial_map <- euclidean_embedding(
#'   dissimilarity_matrix = initial_dist,
#'   ndim = 2, mapping_max_iter = 100,
#'   k0 = 1.0, cooling_rate = 0.01, c_repulsion = 0.01
#' )
#'
#' # Add new points with new measurements
#' new_data <- data.frame(
#'   object = c("D", "D", "E"),
#'   reference = c("A", "B", "C"),
#'   value = c(2.5, 3.0, 1.5)
#' )
#'
#' updated_map <- incremental_embedding(
#'   fixed_positions = initial_map$positions,
#'   new_measurements = new_data,
#'   mapping_max_iter = 100,
#'   k0 = 1.0, cooling_rate = 0.01, c_repulsion = 0.01
#' )
#'
#' @seealso [euclidean_embedding()] for creating initial embeddings,
#'          [list_to_matrix()] for converting list data to matrix format.
#'
#' @importFrom stats runif dist
#' @export
incremental_embedding <- function(fixed_positions,
                                   new_measurements,
                                   object_col = "object",
                                   reference_col = "reference",
                                   value_col = "value",
                                   mapping_max_iter = 1000,
                                   k0 = 1.0,
                                   cooling_rate = 0.01,
                                   c_repulsion = 0.01,
                                   relative_epsilon = 1e-4,
                                   convergence_counter = 5,
                                   n_negative_samples = 5,
                                   convergence_check_freq = 10,
                                   verbose = FALSE) {

  # ===========================================================================
  # INPUT VALIDATION
  # ===========================================================================
  if (!is.matrix(fixed_positions)) {
    stop("fixed_positions must be a matrix")
  }
  if (is.null(rownames(fixed_positions))) {
    stop("fixed_positions must have row names identifying the points")
  }
  if (!is.data.frame(new_measurements)) {
    stop("new_measurements must be a data frame")
  }
  
  required_cols <- c(object_col, reference_col, value_col)
  missing_cols <- setdiff(required_cols, names(new_measurements))
  if (length(missing_cols) > 0) {
    stop("new_measurements is missing required columns: ", 
         paste(missing_cols, collapse = ", "))
  }
  
  if (!is.numeric(mapping_max_iter) || mapping_max_iter < 1) {
    stop("mapping_max_iter must be a positive integer")
  }
  if (!is.numeric(k0) || k0 <= 0) {
    stop("k0 must be a positive number")
  }
  if (!is.numeric(cooling_rate) || cooling_rate <= 0 || cooling_rate >= 1) {
    stop("cooling_rate must be between 0 and 1")
  }
  if (!is.numeric(c_repulsion) || c_repulsion <= 0) {
    stop("c_repulsion must be a positive number")
  }
  
  ndim <- ncol(fixed_positions)
  n_fixed <- nrow(fixed_positions)
  fixed_names <- rownames(fixed_positions)
  
  # Initialize warnings collector
  incremental_warnings <- character(0)
  
  # ===========================================================================
  # IDENTIFY POINT SETS
  # ===========================================================================
  # Get all unique point names from new measurements
  measurement_objects <- as.character(new_measurements[[object_col]])
  measurement_references <- as.character(new_measurements[[reference_col]])
  measurement_names <- unique(c(measurement_objects, measurement_references))
  
  # Identify truly new points (not in fixed positions)
  new_point_names <- setdiff(measurement_names, fixed_names)
  n_new <- length(new_point_names)
  
  # Check for edge case: no new points
  if (n_new == 0) {
    warning("No new points found in measurements. All points already exist in fixed_positions. ",
            "Returning original positions unchanged.")
    incremental_warnings <- c(incremental_warnings,
                              "No new points - measurements only between existing points")
    
    # Return the fixed positions as-is with minimal processing
    p_dist_mat <- as.matrix(stats::dist(fixed_positions))
    rownames(p_dist_mat) <- colnames(p_dist_mat) <- fixed_names
    
    result <- structure(
      list(
        positions = fixed_positions,
        est_distances = p_dist_mat,
        mae = 0,
        iter = 0,
        parameters = list(
          ndim = ndim,
          k0 = k0,
          cooling_rate = cooling_rate,
          c_repulsion = c_repulsion,
          n_negative_samples = n_negative_samples,
          method = "incremental_cpp"
        ),
        convergence = list(
          achieved = TRUE,
          error = 0,
          final_k = k0
        ),
        incremental_info = list(
          n_fixed_points = n_fixed,
          n_new_points = 0,
          n_new_edges = 0,
          new_point_names = character(0),
          warnings = incremental_warnings
        )
      ),
      class = "topolow"
    )
    return(result)
  }
  
  # All point names (fixed first, then new)
  all_names <- c(fixed_names, new_point_names)
  n_total <- length(all_names)
  
  if (verbose) {
    cat(sprintf("Incremental embedding: %d fixed points, %d new points\n", 
                n_fixed, n_new))
  }
  
  # ===========================================================================
  # CHECK FOR NEW-NEW ONLY MEASUREMENTS (NO ANCHOR TO FIXED POINTS)
  # ===========================================================================
  # Check if any new point has measurements to fixed points
  has_anchor <- sapply(new_point_names, function(np) {
    # Check if this new point appears with any fixed point
    obj_matches <- measurement_objects == np & measurement_references %in% fixed_names
    ref_matches <- measurement_references == np & measurement_objects %in% fixed_names
    any(obj_matches | ref_matches)
  })
  
  unanchored_points <- new_point_names[!has_anchor]
  if (length(unanchored_points) > 0) {
    warning("The following new points have no measurements to fixed points and cannot be ",
            "reliably positioned: ", paste(unanchored_points, collapse = ", "))
    incremental_warnings <- c(incremental_warnings,
                              paste("Unanchored new points:", 
                                    paste(unanchored_points, collapse = ", ")))
  }
  
  # ===========================================================================
  # CHECK FOR INSUFFICIENT MEASUREMENTS PER NEW POINT
  # ===========================================================================
  # Count measurements per new point
  measurements_per_new_point <- sapply(new_point_names, function(np) {
    sum(measurement_objects == np | measurement_references == np)
  })
  
  underconstrained_points <- new_point_names[measurements_per_new_point < ndim]
  if (length(underconstrained_points) > 0) {
    warning("The following new points have fewer measurements (", 
            paste(measurements_per_new_point[measurements_per_new_point < ndim], collapse = ", "),
            ") than dimensions (", ndim, "). Their positions may have low accuracy: ",
            paste(underconstrained_points, collapse = ", "))
    incremental_warnings <- c(incremental_warnings,
                              paste("Underconstrained points (measurements < ndim):",
                                    paste(underconstrained_points, collapse = ", ")))
  }
  
  # ===========================================================================
  # BUILD MERGED DISSIMILARITY MATRIX
  # 
  # DECISION 5: CONTRADICTORY DATA PREVENTION
  # Fixed-fixed distances come ONLY from existing positions (computed via 
  # Euclidean distance from fixed_positions). Any fixed-fixed measurements
  # in new_measurements are IGNORED to prevent impossible constraints.
  # 
  # Example: If A and B are fixed at distance 3.0 apart, but new_measurements

  # contains A-B = 5.0, using 5.0 would create an impossible constraint since
  # A and B cannot move. We use 3.0 (the actual geometric distance).
  # ===========================================================================
  # Step 1: Convert new measurements to partial matrix using existing function
  partial_matrix <- list_to_matrix(
    data = new_measurements,
    object_col = object_col,
    reference_col = reference_col,
    value_col = value_col,
    is_similarity = FALSE
  )
  
  # Step 2: Get distances from fixed positions
  fixed_dist_mat <- as.matrix(stats::dist(fixed_positions))
  rownames(fixed_dist_mat) <- colnames(fixed_dist_mat) <- fixed_names
  
  # Step 3: Create full merged matrix
  # Use character matrix to preserve threshold indicators
  merged_matrix <- matrix(NA, nrow = n_total, ncol = n_total)
  rownames(merged_matrix) <- colnames(merged_matrix) <- all_names
  
  # Fill diagonal with zeros
  diag(merged_matrix) <- 0
  
  # Fill in fixed-to-fixed distances (these are known from existing positions)
  # These are AUTHORITATIVE and will NOT be overwritten by new_measurements
  # Vectorized assignment using index matching
  fixed_row_idx <- match(fixed_names, all_names)
  fixed_col_idx <- match(fixed_names, all_names)
  idx_grid <- expand.grid(row = fixed_row_idx, col = fixed_col_idx)
  merged_matrix[as.matrix(idx_grid)] <- as.vector(fixed_dist_mat)
  
  # Step 4: Fill in new measurements ONLY for edges involving at least one new point
  # This prevents contradictory data from overwriting fixed-fixed distances (FIX FOR ISSUE B)
  partial_rownames <- rownames(partial_matrix)
  partial_colnames <- colnames(partial_matrix)
  
  # Get indices for partial matrix elements in merged matrix
  row_map <- match(partial_rownames, all_names)
  col_map <- match(partial_colnames, all_names)
  
  # Track if any fixed-fixed measurements were skipped
  skipped_fixed_fixed <- 0
  
  # For each non-NA element in partial_matrix, set it in merged_matrix
  # ONLY if at least one endpoint is a new point
  for (i in seq_along(partial_rownames)) {
    for (j in seq_along(partial_colnames)) {
      val <- partial_matrix[i, j]
      if (!is.na(val) && partial_rownames[i] != partial_colnames[j]) {
        # Check if this is a fixed-fixed pair
        is_fixed_fixed <- (partial_rownames[i] %in% fixed_names) && 
                          (partial_colnames[j] %in% fixed_names)
        
        if (is_fixed_fixed) {
          # Skip - do not overwrite fixed-fixed distances
          skipped_fixed_fixed <- skipped_fixed_fixed + 1
        } else {
          # At least one endpoint is new - use the measurement
          merged_matrix[row_map[i], col_map[j]] <- val
          merged_matrix[col_map[j], row_map[i]] <- val  # Symmetry
        }
      }
    }
  }
  
  if (skipped_fixed_fixed > 0 && verbose) {
    cat(sprintf("Note: Skipped %d measurements between fixed points (using geometric distances instead)\n",
                skipped_fixed_fixed))
  }
  
  # Convert to character to handle threshold indicators properly
  merged_matrix <- matrix(as.character(merged_matrix), 
                          nrow = n_total, ncol = n_total)
  merged_matrix[merged_matrix == "NA"] <- NA
  rownames(merged_matrix) <- colnames(merged_matrix) <- all_names
  diag(merged_matrix) <- "0"
  
  # ===========================================================================
  # MATRIX REORDERING (same as euclidean_embedding)
  # ===========================================================================
  # Extract numeric values for clustering, handling threshold indicators
  numeric_matrix <- matrix(NA_real_, n_total, n_total)
  non_na_mask <- !is.na(merged_matrix)
  numeric_matrix[non_na_mask] <- as.numeric(gsub("^[<>]", "", merged_matrix[non_na_mask]))
  
  # Create ordering to concentrate largest values in corners
  tryCatch({
    diag(numeric_matrix) <- NA
    row_means <- rowMeans(numeric_matrix, na.rm = TRUE)
    col_means <- colMeans(numeric_matrix, na.rm = TRUE)
    avg_dissim <- (row_means + col_means) / 2
    avg_dissim[is.nan(avg_dissim)] <- 0
    
    if (sum(avg_dissim > 0) > 1) {
      new_order <- order(avg_dissim)
      merged_matrix <- merged_matrix[new_order, new_order]
      all_names <- all_names[new_order]
      
      if (verbose) {
        cat("Matrix reordered for spectral pattern\n")
      }
    }
  }, error = function(e) {
    if (verbose) cat("Could not reorder matrix:", e$message, "\n")
  })
  
  # Update fixed mask after reordering
  fixed_mask <- as.integer(all_names %in% fixed_names)
  n <- n_total
  
  # Create indices of fixed and new points (0-based for C++)
  # These are used for the extra repulsion step (FIX FOR ISSUE A)
  fixed_indices_cpp <- as.integer(which(fixed_mask == 1) - 1L)
  new_indices_cpp <- as.integer(which(fixed_mask == 0) - 1L)
  
  # ===========================================================================
  # PREPROCESSING FOR C++ (same as euclidean_embedding)
  # ===========================================================================
  is_na_matrix <- is.na(merged_matrix)
  node_degrees <- rowSums(!is_na_matrix)
  
  # Parse to numeric values and threshold codes
  distances_numeric <- matrix(Inf, n, n)
  threshold_mask <- matrix(0L, n, n)
  
  non_na_mask <- !is_na_matrix
  
  # Handle character matrix with thresholds
  gt_mask <- non_na_mask & startsWith(merged_matrix, ">")
  lt_mask <- non_na_mask & startsWith(merged_matrix, "<")
  normal_mask <- non_na_mask & !gt_mask & !lt_mask
  
  threshold_mask[gt_mask] <- 1L
  threshold_mask[lt_mask] <- -1L
  
  if (any(gt_mask)) {
    distances_numeric[gt_mask] <- as.numeric(sub("^>", "", merged_matrix[gt_mask]))
  }
  if (any(lt_mask)) {
    distances_numeric[lt_mask] <- as.numeric(sub("^<", "", merged_matrix[lt_mask]))
  }
  if (any(normal_mask)) {
    distances_numeric[normal_mask] <- as.numeric(merged_matrix[normal_mask])
  }
  
  # ===========================================================================
  # BUILD COO EDGE LIST (UPPER TRIANGLE, FILTERED FOR NEW POINTS)
  # 
  # DECISION 2: EDGE FILTERING
  # Only edges involving at least one new point are passed to C++.
  # Fixed-Fixed edges are excluded since neither endpoint can move.
  # This provides significant speedup: for 1000 fixed + 10 new points,
  # we process ~20,000 edges instead of ~500,000.
  # ===========================================================================
  upper_tri_mask <- upper.tri(distances_numeric)
  valid_edges <- upper_tri_mask & (distances_numeric != Inf)
  
  # Get all edge indices
  edge_indices <- which(valid_edges, arr.ind = TRUE)
  
  if (nrow(edge_indices) == 0) {
    stop("No valid measurements found in new_measurements")
  }
  
  # Filter to only edges involving at least one new point
  edge_involves_new <- (fixed_mask[edge_indices[, 1]] == 0) | 
                       (fixed_mask[edge_indices[, 2]] == 0)
  
  # Keep only edges involving new points for the main optimization
  filtered_edge_indices <- edge_indices[edge_involves_new, , drop = FALSE]
  
  n_new_edges <- nrow(filtered_edge_indices)
  
  if (n_new_edges == 0) {
    stop("No edges involving new points found. Check that new_measurements ",
         "contains measurements between new and existing points.")
  }
  
  if (verbose) {
    cat(sprintf("Processing %d edges involving new points (out of %d total edges)\n",
                n_new_edges, nrow(edge_indices)))
  }
  
  # Convert to 0-based indexing for C++
  edge_i <- as.integer(filtered_edge_indices[, 1] - 1L)
  edge_j <- as.integer(filtered_edge_indices[, 2] - 1L)
  
  # Extract corresponding distances and thresholds
  linear_idx <- filtered_edge_indices[, 1] + (filtered_edge_indices[, 2] - 1L) * n
  edge_dist <- distances_numeric[linear_idx]
  edge_thresh <- threshold_mask[linear_idx]
  
  # Create flattened has_measurement vector for negative sampling
  has_measurement <- !is_na_matrix
  has_measurement_flat <- as.integer(as.vector(t(has_measurement)))
  
  # ===========================================================================
  # INITIALIZATION OF POINT POSITIONS
  # ===========================================================================
  # Initialize with fixed positions for known points, random walk for new points
  
  # Calculate init_step similar to euclidean_embedding
  dissim_numeric <- suppressWarnings(as.numeric(as.character(merged_matrix)))
  dissim_numeric[is.na(dissim_numeric)] <- NA
  init_step <- max(dissim_numeric, na.rm = TRUE) / n
  
  # Create initial positions matrix
  initial_positions <- matrix(0, nrow = n, ncol = ndim)
  rownames(initial_positions) <- all_names
  
  # For fixed points: use their known positions (reordered to match merged_matrix)
  fixed_indices_in_merged <- which(fixed_mask == 1)
  for (idx in fixed_indices_in_merged) {
    point_name <- all_names[idx]
    initial_positions[idx, ] <- fixed_positions[point_name, ]
  }
  
  # For new points: initialize using incremental random walk from nearest fixed point
  new_indices_in_merged <- which(fixed_mask == 0)
  
  if (length(new_indices_in_merged) > 0) {
    for (idx in new_indices_in_merged) {
      point_name <- all_names[idx]
      row_vals <- distances_numeric[idx, ]
      col_vals <- distances_numeric[, idx]
      
      # Get measured distances to fixed points
      fixed_point_dists <- rep(Inf, n)
      for (f_idx in fixed_indices_in_merged) {
        if (row_vals[f_idx] != Inf) {
          fixed_point_dists[f_idx] <- row_vals[f_idx]
        } else if (col_vals[f_idx] != Inf) {
          fixed_point_dists[f_idx] <- col_vals[f_idx]
        }
      }
      
      # Find nearest fixed point with a measurement
      nearest_fixed_idx <- which.min(fixed_point_dists)
      
      if (is.finite(fixed_point_dists[nearest_fixed_idx])) {
        # Initialize near the nearest fixed point with random offset
        nearest_pos <- initial_positions[nearest_fixed_idx, ]
        random_offset <- stats::runif(ndim, -init_step, init_step)
        initial_positions[idx, ] <- nearest_pos + random_offset
      } else {
        # No measurement to fixed points - use incremental random walk
        random_steps <- stats::runif(ndim, 0, 2 * init_step)
        if (idx > 1) {
          initial_positions[idx, ] <- initial_positions[idx - 1, ] + random_steps
        } else {
          initial_positions[idx, ] <- random_steps
        }
      }
    }
  }
  
  # ===========================================================================
  # CALL C++ OPTIMIZATION
  # ===========================================================================
  if (verbose) cat("Starting C++ incremental optimization...\n")
  
  start_time <- Sys.time()
  
  cpp_result <- optimize_layout_incremental_cpp(
    initial_positions = initial_positions,
    edge_i = edge_i,
    edge_j = edge_j,
    edge_dist = edge_dist,
    edge_thresh = edge_thresh,
    has_measurement_flat = has_measurement_flat,
    degrees = as.integer(node_degrees),
    fixed_point_mask = as.integer(fixed_mask),
    fixed_indices = fixed_indices_cpp,
    new_indices = new_indices_cpp,
    n_points = as.integer(n),
    n_iter = as.integer(mapping_max_iter),
    k0 = k0,
    cooling_rate = cooling_rate,
    c_repulsion = c_repulsion,
    n_negative_samples = as.integer(n_negative_samples),
    relative_epsilon = relative_epsilon,
    convergence_window = as.integer(convergence_counter),
    convergence_check_freq = as.integer(convergence_check_freq),
    verbose = verbose
  )
  
  end_time <- Sys.time()
  
  if (verbose) {
    cat(sprintf("Optimization finished in %.2f seconds.\n",
                as.numeric(difftime(end_time, start_time, units = "secs"))))
  }
  
  # ===========================================================================
  # POST-PROCESSING
  # ===========================================================================
  positions <- cpp_result$positions
  rownames(positions) <- all_names
  
  # Calculate distances and metrics
  p_dist_mat <- as.matrix(stats::dist(positions))
  rownames(p_dist_mat) <- colnames(p_dist_mat) <- all_names
  
  # Calculate MAE only on edges involving new points (what we actually optimized)
  mae_values <- numeric(n_new_edges)
  for (e in seq_len(n_new_edges)) {
    i <- filtered_edge_indices[e, 1]
    j <- filtered_edge_indices[e, 2]
    target <- edge_dist[e]
    actual <- p_dist_mat[i, j]
    mae_values[e] <- abs(target - actual)
  }
  mae <- mean(mae_values)
  
  # ===========================================================================
  # CREATE RESULT OBJECT
  # ===========================================================================
  result <- structure(
    list(
      positions = positions,
      est_distances = p_dist_mat,
      mae = mae,
      iter = cpp_result$iterations,
      parameters = list(
        ndim = ndim,
        k0 = k0,
        cooling_rate = cooling_rate,
        c_repulsion = c_repulsion,
        n_negative_samples = n_negative_samples,
        method = "incremental_cpp"
      ),
      convergence = list(
        achieved = cpp_result$converged,
        error = cpp_result$final_mae,
        final_k = cpp_result$final_k
      ),
      incremental_info = list(
        n_fixed_points = n_fixed,
        n_new_points = n_new,
        n_new_edges = n_new_edges,
        new_point_names = new_point_names,
        warnings = incremental_warnings
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
#' the entire workflow from parameter optimization to final embedding,  with comprehensive diagnostic tracking and visualization.
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
#'        Default: 150
#' @param max_cores Integer. Maximum number of cores to use. Default: NULL (auto-detect)
#' @param folds Integer. Number of cross-validation folds. Default: 20
#' @param mapping_max_iter Integer. Maximum iterations for final embedding. 
#'        Half this value is used for parameter search. Default: 500
#' @param opt_subsample Integer or NULL. If specified, uses subsampling during parameter
#'   optimization (both initial and adaptive) to reduce computational cost. Randomly
#'   samples this many points for each parameter evaluation. Final embedding always
#'   uses the full dataset. Default: NULL (no subsampling).
#'   Recommended for large datasets (>300 points). See \code{\link{initial_parameter_optimization}}
#'   for details.
#' @param clean_intermediate Logical. Whether to remove intermediate files. Default: TRUE
#' @param verbose Character. Verbosity level: "off" (no output), "standard" (progress updates),
#'        or "full" (detailed output including from internal functions). Default: "standard"
#' @param fallback_to_defaults Logical. Whether to use default parameters if 
#'        optimization fails. Default: FALSE
#' @param save_results Logical. Whether to save the final positions as CSV. Default: FALSE
#' @param create_diagnostic_plots Logical. Whether to create diagnostic and trace plots
#'        showing the parameter optimization process and embedding quality. Default: FALSE
#' @param diagnostic_plot_types Character vector. Which plot types to create.
#'        Options: "all", "parameter_search", "convergence", "quality", "cv_errors".
#'        Default: "all"
#' 
#' @return A list containing:
#'   \item{positions}{Matrix of optimized coordinates}
#'   \item{est_distances}{Matrix of estimated distances}
#'   \item{mae}{Mean absolute error}
#'   \item{optimal_params}{List of optimal parameters found, including cross-validation MAE during optimization}
#'   \item{optimization_summary}{Summary of the optimization process}
#'   \item{data_characteristics}{Summary of input data characteristics}
#'   \item{runtime}{Total runtime in seconds}
#'   \item{all_samples}{Data frame of all parameter evaluations (if create_diagnostic_plots=TRUE)}
#'   \item{diagnostic_plots}{List of ggplot objects (if create_diagnostic_plots=TRUE)}
#'   \item{dissimilarity_matrix}{Input dissimilarity matrix (if create_diagnostic_plots=TRUE)}
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
#' # Example 6: Basic usage
#' test_data <- data.frame(
#'   object = rep(paste0("Obj", 1:4), each = 4),
#'   reference = rep(paste0("Ref", 1:4), 4),
#'   score = sample(c(1, 2, 4, 8, 16, 32, 64, "<1", ">12"), 16, replace = TRUE)
#' )
#' dist_mat <- list_to_matrix(
#'   data = test_data,
#'   object_col = "object",
#'   reference_col = "reference",
#'   value_col = "score",
#'   is_similarity = TRUE
#' )
#' \dontrun{
#' # Basic usage with diagnostics
#' result <- Euclidify(
#'   dissimilarity_matrix = dist_mat,
#'   output_dir = tempdir(),
#'   create_diagnostic_plots = TRUE
#' )
#' 
#' # View diagnostic report
#' report <- create_diagnostic_report(result)
#' cat(report, sep = "\n")
#' 
#' # Access specific diagnostic plots
#' print(result$diagnostic_plots$parameter_search)
#' print(result$diagnostic_plots$convergence)
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
                      opt_subsample = NULL,
                      clean_intermediate = TRUE,
                      verbose = "standard",
                      fallback_to_defaults = FALSE,
                      save_results = FALSE,
                      create_diagnostic_plots = FALSE,
                      diagnostic_plot_types = "all") {
  
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
  # Convert to numeric matrix for eigen decomposition
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
  
  # Suggest optimal dimensions based on eigenvalues
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
      opt_subsample = opt_subsample,
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
  
  # Create diagnostic plots if requested
  if (create_diagnostic_plots) {
    if (verbose_main) cat("\nGenerating diagnostic plots...\n")
    
    # Read the parameter samples that were already saved
    samples_file <- file.path(optimization_dir, "model_parameters", 
                              "auto_optimization_model_parameters.csv")
    
    if (!file.exists(samples_file)) {
      if (verbose_main) {
        cat("  ERROR: Parameter samples file not found at:", samples_file, "\n")
        cat("  Checking directory contents:\n")
        param_dir <- file.path(optimization_dir, "model_parameters")
        if (dir.exists(param_dir)) {
          cat("  Files in", param_dir, ":\n")
          print(list.files(param_dir))
        } else {
          cat("  Directory doesn't exist:", param_dir, "\n")
        }
      }
    } else {
      all_samples <- read.csv(samples_file, stringsAsFactors = FALSE)
      
      # Create diagnostic plots directory
      diag_dir <- file.path(output_dir, "diagnostics")
      dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)
      
      # Prepare data for diagnostic plots
      results$all_samples <- all_samples  # Store for later use
      results$dissimilarity_matrix <- dissimilarity_matrix  # Needed for quality plots
      
      # Call diagnostic plotting function
      tryCatch({
        results$diagnostic_plots <- plot_euclidify_diagnostics(
          euclidify_result = results,
          plot_types = diagnostic_plot_types,
          save_plots = TRUE,
          output_dir = diag_dir,
          return_plots = TRUE
        )
        
        # Also create text report
        report_file <- file.path(diag_dir, "diagnostic_report.txt")
        create_diagnostic_report(results, output_file = report_file)
        
        if (verbose_main) {
          cat("  Diagnostic plots saved to:", diag_dir, "\n")
          cat("  Diagnostic report saved to:", report_file, "\n")
        }
      }, error = function(e) {
        if (verbose_main) {
          cat("  Warning: Could not create diagnostic plots:", e$message, "\n")
        }
      })
    }
  }
  
  # Clean up intermediate files
  if (clean_intermediate) {
    unlink(optimization_dir, recursive = TRUE)
    if (verbose_main) cat("\nCleaned up intermediate files\n")
  }
  
  return(results)
}
