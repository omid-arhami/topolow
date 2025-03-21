# Copyright (c) 2024 Omid Arhami omid.arhami@uga.edu
# License: free of charge access granted to any academic researcher to use this software for non-commercial, academic research purposes **only**.  Nobody may modify, distribute, sublicense, or publicly share the Software or any derivative works, until the paper is published by the original authors.  The Software is provided "as is" without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose and noninfringement.  In no event shall the authors or copyright holders be liable for any claim, damages or other liability, whether in an action of contract, tort or otherwise, arising from, out of or in connection with the Software or the use or other dealings in the Software.

# R/core.R

#' TopoLow Core Functions
#' 
#' Core implementations of the TopoLow algorithm for mapping distances in high dimensions.
#' This file contains the main optimization functions and their variants.
#' 
#' @keywords internal
#' @name topolow-package
NULL


#' Process ideal distance values
#' 
#' Helper function that processes distance values that may include threshold indicators
#' (< or >). Used internally by the optimization functions to handle thresholded
#' measurements.
#' 
#' @param reported_distance Character or numeric. The reported distance value, possibly 
#'        with < or > prefix
#' @param distance Numeric. The calculated distance to compare against
#' @return Numeric or NULL. Returns the processed distance value if the threshold 
#'         condition is met, NULL otherwise
#' @keywords internal
process_ideal_distance <- function(reported_distance, distance) {
  if (!is.numeric(distance)) {
    stop("'distance' must be numeric")
  }
  
  if (is.character(reported_distance) && startsWith(reported_distance, ">")) {
    numeric_part <- as.numeric(sub(">", "", reported_distance))
    if (distance < numeric_part) {
      return(numeric_part)
    } else {
      return(NULL)
    }
  } else if (is.character(reported_distance) && startsWith(reported_distance, "<")) {
    numeric_part <- as.numeric(sub("<", "", reported_distance))
    if (distance > numeric_part) {
      return(numeric_part)
    } else {
      return(NULL)
    }
  } else {
    tryCatch({
      return(as.numeric(reported_distance))
    }, error = function(e) {
      stop("'reported_distance' must be numeric or a character starting with > or <")
    })
  }
}


#' Process distance matrix for convergence error calculations
#' 
#' Helper function that processes elements of the distance matrix for calculating 
#' convergence error. Handles threshold indicators and NA values.
#' 
#' @param reported_distance Character or numeric. The reported distance value
#' @param distance Numeric. The calculated distance to compare against
#' @return Numeric or NA. Returns processed distance if valid, NA otherwise
#' @keywords internal
process_distance_matrix <- function(reported_distance, distance) {
  if (is.na(reported_distance) || is.na(distance)) {
    return(NA)
  }
  
  if (is.character(reported_distance) && startsWith(reported_distance, ">")) {
    numeric_part <- as.numeric(sub(">", "", reported_distance))
    if (distance < numeric_part) {
      return(numeric_part)
    } else {
      return(NA)
    }
  } else if (is.character(reported_distance) && startsWith(reported_distance, "<")) {
    numeric_part <- as.numeric(sub("<", "", reported_distance))
    if (distance > numeric_part) {
      return(numeric_part)
    } else {
      return(NA)
    }
  } else {
    tryCatch({
      return(as.numeric(reported_distance))
    }, error = function(e) {
      return(NA)
    })
  }
}

# Vectorized function to process distance matrix for convergence error calculations
vectorized_process_distance_matrix <- function(distance_matrix, p_dist_mat) {
  # Create result matrix
  result <- matrix(NA, nrow = nrow(distance_matrix), ncol = ncol(distance_matrix))
  
  # Convert to character matrix for string operations
  char_matrix <- matrix(as.character(distance_matrix), nrow = nrow(distance_matrix))
  
  # Process greater than cases (">x")
  gt_mask <- grepl("^>", char_matrix)
  if (any(gt_mask)) {
    gt_indices <- which(gt_mask)
    numeric_parts <- suppressWarnings(as.numeric(gsub("^>", "", char_matrix[gt_indices])))
    valid <- !is.na(numeric_parts) & !is.na(p_dist_mat[gt_indices]) & 
      p_dist_mat[gt_indices] < numeric_parts
    result[gt_indices[valid]] <- numeric_parts[valid]
  }
  
  # Process less than cases ("<x")
  lt_mask <- grepl("^<", char_matrix)
  if (any(lt_mask)) {
    lt_indices <- which(lt_mask)
    numeric_parts <- suppressWarnings(as.numeric(gsub("^<", "", char_matrix[lt_indices])))
    valid <- !is.na(numeric_parts) & !is.na(p_dist_mat[lt_indices]) & 
      p_dist_mat[lt_indices] > numeric_parts
    result[lt_indices[valid]] <- numeric_parts[valid]
  }
  
  # Process other values (convert to numeric)
  other_mask <- !is.na(char_matrix) & !gt_mask & !lt_mask
  if (any(other_mask)) {
    other_indices <- which(other_mask)
    result[other_indices] <- suppressWarnings(as.numeric(char_matrix[other_indices]))
  }
  
  return(result)
}


#' Main TopoLow algorithm implementation
#' 
#' TopoLow (Topological Optimization for Low-Dimensional Mapping) optimizes point positions in n-dimensional 
#' space to match a target distance matrix. The algorithm uses a physics-inspired approach with 
#' spring and repulsive forces to find optimal point configurations while handling missing 
#' and thresholded measurements.
#'
#' @details
#' The algorithm iteratively updates point positions using:
#' * Spring forces between points with measured distances
#' * Repulsive forces between points without measurements
#' * Modified forces for thresholded measurements (< or >)
#' * Adaptive spring constant that decays over iterations
#' * Convergence monitoring based on relative error change
#'
#' Valid parameter ranges and constraints:
#' * ndim: Positive integer, typically 2-20.
#' * k0: Initial spring constant, positive numeric > 0. Typical range: 0.1-30
#'      Controls initial force strength
#' * cooling_rate: Spring and repulsion decay rate, numeric between 0 and 1. Typical range: 0.0001-0.1
#'          Controls how quickly spring forces weaken
#' * c_repulsion: Repulsion constant, positive numeric > 0. Typical range: 0.00001-0.1
#'      Controls strength of repulsive forces
#' * relative_epsilon: Positive numeric, typically 1e-9 to 1e-3
#'                    Smaller values require more iterations but give higher precision
#' * convergence_counter: Positive integer, typically 5-20
#'                       Higher values do not necessarily lead to a better convergence
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
#'        Default is TRUE.
#' @param verbose Logical. Whether to print progress messages. Default is TRUE.
#'
#' @return A list with class "topolow" containing:
#' \itemize{
#'   \item positions: Matrix of optimized point coordinates
#'   \item est_distances: Matrix of distances in the optimized configuration
#'   \item mae: Mean absolute error between target and optimized distances
#'   \item r: Pearson correlation between target and optimized distances
#'   \item iter: Number of iterations performed
#'   \item parameters: List of input parameters used
#'   \item convergence: List with convergence status and final error
#' }
#'
#' @examples
#' # Create a simple distance matrix
#' dist_mat <- matrix(c(0, 2, 3, 2, 0, 4, 3, 4, 0), nrow=3)
#' 
#' # Run TopoLow in 2D
#' result <- create_topolow_map(dist_mat, ndim=2, mapping_max_iter=1000, 
#'                       k0=1.0, cooling_rate=0.001, c_repulsion=0.1)
#'                       
#' # Plot results
#' plot(result$positions)
#'
#' @export
create_topolow_map <- function(distance_matrix, 
                         ndim, 
                         mapping_max_iter, 
                         k0, 
                         cooling_rate, 
                         c_repulsion, 
                         relative_epsilon = 1e-4,
                         convergence_counter = 5,
                         initial_positions = NULL,
                         write_positions_to_csv = TRUE,
                         verbose = FALSE) {
  
  # Input validation with informative messages
  if (!is.matrix(distance_matrix)) {
    stop("distance_matrix must be a matrix")
  }
  if (nrow(distance_matrix) != ncol(distance_matrix)) {
    stop("distance_matrix must be square")
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
    if (nrow(initial_positions) != nrow(distance_matrix)) {
      stop("initial_positions must have same number of rows as distance_matrix")
    }
    if (ncol(initial_positions) != ndim) {
      stop("initial_positions must have ndim columns")
    }
  }
  
  # Initialize variables
  n <- nrow(distance_matrix)
  
  # OPTIMIZATION 1: Pre-compute NA status once
  is_na_matrix <- is.na(distance_matrix)
  node_degrees <- rowSums(!is_na_matrix)
  node_degrees_1 <- node_degrees + 1
  
  # OPTIMIZATION 2: Replace ifelse with direct indexing
  distances <- matrix(Inf, n, n)
  distances[!is_na_matrix] <- distance_matrix[!is_na_matrix]
  
  # Initialize positions (unchanged)
  positions <- if(is.null(initial_positions)) {
    distance_matrix_numeric <- suppressWarnings(as.numeric(as.character(distance_matrix)))
    distance_matrix_numeric[is.na(distance_matrix_numeric)] <- NA
    init_step <- max(distance_matrix_numeric, na.rm=TRUE) / n
    pos <- matrix(0, n, ndim)
    for (i in 2:n) {
      pos[i, ] <- pos[i-1, ] + stats::runif(ndim, 0, 2*init_step)
    }
    pos
  } else {
    initial_positions
  }
  
  rownames(positions) <- rownames(distance_matrix)
  
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
  
  # OPTIMIZATION 4: Pre-compute which pairs have measured distances
  has_measurement <- matrix(FALSE, nrow(point_pairs), 1)
  for(idx in 1:nrow(point_pairs)) {
    i <- point_pairs[idx, 1]
    j <- point_pairs[idx, 2]
    has_measurement[idx] <- distances[i, j] != Inf
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
      
      # OPTIMIZATION 6: Use pre-computed measurement status
      if(shuffled_has_measurement[pair_idx]) {
        ideal_distance <- distances[i, j]
        ideal_distance_processed <- process_ideal_distance(ideal_distance, distance)
        
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
    
    # OPTIMIZATION 8: Check for numerical stability less frequently
    if(iter %% 5 == 0) {
      # Vectorized check for NA/Inf values
      if(any(is.na(positions)) || any(is.infinite(positions))) {
        warning("Numerical instability detected at iteration ", iter)
        stop <- TRUE
        break
      }
    }
    
    if(stop) break
    
    # Calculate current distances and convergence error
    p_dist_mat <- coordinates_to_matrix(positions)
    
    # Process thresholded elements for convergence error
    # distance_matrix_convergence_error <- mapply(
    #   process_distance_matrix, 
    #   distance_matrix, 
    #   p_dist_mat
    # )
    distance_matrix_convergence_error <- vectorized_process_distance_matrix(distance_matrix, p_dist_mat)
    
    distance_matrix_convergence_error <- matrix(
      distance_matrix_convergence_error, 
      nrow = nrow(distance_matrix)
    )
    
    # OPTIMIZATION 9: Vectorized convergence error calculation
    distance_vec <- as.vector(as.numeric(distance_matrix_convergence_error))
    p_dist_vec <- as.vector(p_dist_mat)
    valid_indices <- !is.na(distance_vec)
    convergence_error <- mean(
      abs(distance_vec[valid_indices] - p_dist_vec[valid_indices])
    )
    
    # Calculate evaluation metrics similarly
    distance_vec_raw <- suppressWarnings(as.vector(as.numeric(distance_matrix)))
    valid_indices_raw <- !is.na(distance_vec_raw)
    mae <- mean(abs(distance_vec_raw[valid_indices_raw] - p_dist_vec[valid_indices_raw]))
    
    if (verbose) {
      cat(sprintf(
        "Iteration=%d, MAE=%.4f, convergence_error=%.4f\n", 
        iter, mae, convergence_error
      ))
    }
    
    # Check convergence (relatively unchanged)
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
  
  # Save positions if requested (unchanged)
  if(write_positions_to_csv) {
    csv_filename <- sprintf(
      "Positions_dim_%d_k0_%.4f_cooling_%.4f_c_repulsion_%.4f.csv",
      ndim, k0, cooling_rate, c_repulsion
    )
    utils::write.csv(positions, file = csv_filename, row.names = TRUE)
  }
  
  # Create result object (unchanged)
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


#' Print method for topolow objects
#' 
#' Provides a concise display of key optimization results including dimensions,
#' iterations, error metrics and convergence status.
#' 
#' @param x A topolow object returned by create_topolow_map()
#' @param ... Additional arguments passed to print (not used)
#' 
#' @examples
#' dist_mat <- matrix(c(0, 2, 3, 2, 0, 4, 3, 4, 0), nrow=3)
#' result <- create_topolow_map(dist_mat, ndim=2, mapping_max_iter=100, k0=1.0, cooling_rate=0.001, c_repulsion=0.1)
#' print(result)
#' 
#' @export
print.topolow <- function(x, ...) {
  cat("TopoLow optimization result:\n")
  cat(sprintf("Dimensions: %d\n", x$parameters$ndim))
  cat(sprintf("Iterations: %d\n", x$iter))
  cat(sprintf("MAE: %.4f\n", x$mae))
  cat(sprintf("Convergence achieved: %s\n", x$convergence$achieved))
  cat(sprintf("Final convergence error: %.4f\n", x$convergence$error))
}


#' Summary method for topolow objects
#' 
#' Provides a detailed summary of the optimization results including parameters,
#' convergence and performance metrics.
#' 
#' @param object A topolow object returned by create_topolow_map()
#' @param ... Additional arguments passed to summary (not used)
#' 
#' @examples
#' dist_mat <- matrix(c(0, 2, 3, 2, 0, 4, 3, 4, 0), nrow=3)
#' result <- create_topolow_map(dist_mat, ndim=2, mapping_max_iter=100, k0=1.0, cooling_rate=0.001, c_repulsion=0.1)
#' summary(result)
#' 
#' @export
summary.topolow <- function(object, ...) {
  print(object)
  cat("\nParameters:\n")
  cat(sprintf("k0: %.4f\n", object$parameters$k0))
  cat(sprintf("cooling_rate: %.4f\n", object$parameters$cooling_rate))
  cat(sprintf("c_repulsion: %.4f\n", object$parameters$c_repulsion))
}


#' Sigmoid transform function for threshold handling
#' 
#' Helper function that implements a sigmoid transformation used in the Smith variant
#' for smooth threshold handling.
#' 
#' @param x Numeric. Input value
#' @return Numeric. Transformed value in (0,1) range
#' @keywords internal
gmultiple <- function(x) {
  if (!is.numeric(x)) {
    stop("Input must be numeric")
  }
  1/(1 + exp(-10*x))
}

