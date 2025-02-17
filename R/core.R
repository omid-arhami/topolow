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
#' * cooling_rate: Spring decay rate, numeric between 0 and 1. Typical range: 0.0001-0.1
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
#' @param max_iter Integer. Maximum number of optimization iterations.
#' @param k0 Numeric. Initial spring constant controlling spring forces.
#' @param cooling_rate Numeric. Rate of spring constant decay per iteration (0 < cooling_rate < 1).
#' @param c_repulsion Numeric. Repulsion constant controlling repulsive forces.
#' @param relative_epsilon Numeric. Convergence threshold for relative change in error.
#'        Default is 1e-4.
#' @param convergence_counter Integer. Number of iterations below threshold before declaring
#'        convergence. Default is 10.
#' @param initial_positions Matrix or NULL. Optional starting coordinates. If NULL,
#'        random initialization is used. Matrix should have nrow = nrow(distance_matrix)
#'        and ncol = ndim.
#' @param write_positions_to_csv Logical. Whether to save point positions to CSV file.
#'        Default is TRUE.
#' @param verbose Logical. Whether to print progress messages. Default is TRUE.
#' @param trace_sse Logical. Whether to track sum of squared errors. Default is TRUE.
#'
#' @return A list with class "topolow" containing:
#' \itemize{
#'   \item positions: Matrix of optimized point coordinates
#'   \item est_distances: Matrix of distances in the optimized configuration
#'   \item mae: Mean absolute error between target and optimized distances
#'   \item r: Pearson correlation between target and optimized distances
#'   \item iter: Number of iterations performed
#'   \item trace_convergence_error_df: Data frame tracking convergence
#'   \item parameters: List of input parameters used
#'   \item convergence: List with convergence status and final error
#' }
#'
#' @seealso 
#' \code{\link{topolow_Smith_obj}} for a variant based on the squishing function in Smith et al 2004 for HI assay data
#'
#' @examples
#' # Create a simple distance matrix
#' dist_mat <- matrix(c(0, 2, 3, 2, 0, 4, 3, 4, 0), nrow=3)
#' 
#' # Run TopoLow in 2D
#' result <- topolow_full(dist_mat, ndim=2, max_iter=1000, 
#'                       k0=1.0, cooling_rate=0.001, c_repulsion=0.1)
#'                       
#' # Plot results
#' plot(result$positions)
#'
#' @export
topolow_full <- function(distance_matrix, 
                         ndim, 
                         max_iter, 
                         k0, 
                         cooling_rate, 
                         c_repulsion, 
                         relative_epsilon = 1e-4,
                         convergence_counter = 5,
                         initial_positions = NULL,
                         write_positions_to_csv = TRUE,
                         verbose = FALSE,
                         trace_sse = TRUE) {
  
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
  if (!is.numeric(max_iter) || max_iter < 1 || max_iter != round(max_iter)) {
    stop("max_iter must be a positive integer")
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

  # Sort input matrix
  #distance_matrix <- sort_matrix_by_year(distance_matrix)
  
  # Initialize variables
  n <- nrow(distance_matrix)
  node_degrees <- rowSums(!is.na(distance_matrix))
  node_degrees_1 <- node_degrees + 1
  distances <- ifelse(is.na(distance_matrix), Inf, distance_matrix)
  
  # Initialize positions
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
  
  trace_convergence_error_df <- data.table::data.table(
    iteration = numeric(max_iter), 
    convergence_error = NULL
  )

  # Create a list of all valid (i,j) pairs maintaining i < j
  point_pairs <- list()
  for(i in 1:(n-1)) {
    for(j in (i+1):n) {
      point_pairs[[length(point_pairs) + 1]] <- c(i,j)
    }
  }
  
  # Convert to matrix
  point_pairs <- do.call(rbind, point_pairs)
  
  ################## Main optimization loop
  for(iter in 1:max_iter) {
    k_2 <- k / 2
    stop <- FALSE
    
    # Randomize point pair ordering for this iteration
    random_order <- sample(nrow(point_pairs))
    shuffled_pairs <- point_pairs[random_order,]
    
    # Calculate forces between pairs in random order
    for(pair_idx in 1:nrow(shuffled_pairs)) {
      i <- shuffled_pairs[pair_idx, 1]
      j <- shuffled_pairs[pair_idx, 2]

      ideal_distance <- distances[i, j]
      delta <- positions[j,] - positions[i,]
      distance <- sqrt(sum(delta^2))
      distance_01 <- distance + 0.01
      
      # Spring forces for measured distances
      if(ideal_distance != Inf) {
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
      
      # Check for numerical stability
      if(any(is.na(positions)) || any(is.infinite(positions))) {
        warning("Numerical instability detected at iteration ", iter)
        stop <- TRUE
        break
      }
      if(stop) break
    }
    
    if(stop) break
    
    # Calculate current distances and convergence error
    p_dist_mat <- coordinates_to_matrix(positions)
    
    # Process thresholded elements for convergence error
    distance_matrix_convergence_error <- mapply(
      process_distance_matrix, 
      distance_matrix, 
      p_dist_mat
    )
    distance_matrix_convergence_error <- matrix(
      distance_matrix_convergence_error, 
      nrow = nrow(distance_matrix)
    )
    
    # Calculate convergence error
    convergence_error_df <- data.table::data.table(
      distance_matrix_convergence_error = suppressWarnings(
        as.vector(as.numeric(distance_matrix_convergence_error))
      ),
      p_dist_mat = as.vector(p_dist_mat)
    )
    convergence_error_df <- na.omit(convergence_error_df)
    convergence_error <- mean(
      abs(convergence_error_df$distance_matrix_convergence_error - 
            convergence_error_df$p_dist_mat), 
      na.rm = TRUE
    )
    
    # Calculate evaluation metrics
    evaldf <- data.table::data.table(
      distance_matrix = suppressWarnings(
        as.vector(as.numeric(distance_matrix))
      ),
      p_dist_mat = as.vector(p_dist_mat)
    )
    evaldf <- na.omit(evaldf)
    mae <- mean(abs(evaldf$distance_matrix - evaldf$p_dist_mat), na.rm = TRUE)
    
    # Record trace
    trace_convergence_error_df[iter, "iteration"] <- iter
    trace_convergence_error_df[iter, "convergence_error"] <- convergence_error
    
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
            "Please check model parameters k, decay, and c_repulsion.\n"
          ))
          cat(sprintf(
            "dim=%d, k0=%.4f, decay=%.4f, c_repulsion=%.4f, MAE=%.4f, convergence_error=%.4f\n",
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
    csv_filename <- sprintf(
      "Positions_dim_%d_k0_%.4f_decay_%.4f_c_repulsion_%.4f.csv",
      ndim, k0, cooling_rate, c_repulsion
    )
    utils::write.csv(positions, file = csv_filename, row.names = TRUE)
  }
  
  # Calculate final correlation
  pearson_corr <- stats::cor(
    evaldf$p_dist_mat, 
    evaldf$distance_matrix, 
    method = "pearson"
  )
  
  # Create result object
  result <- structure(
    list(
      positions = positions,
      est_distances = p_dist_mat,
      mae = mae,
      r = pearson_corr,
      iter = iter,
      trace_convergence_error_df = trace_convergence_error_df,
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
#' @param x A topolow object returned by topolow_full() or topolow_Smith_obj()
#' @param ... Additional arguments passed to print (not used)
#' 
#' @examples
#' dist_mat <- matrix(c(0, 2, 3, 2, 0, 4, 3, 4, 0), nrow=3)
#' result <- topolow_full(dist_mat, ndim=2, max_iter=100, k0=1.0, cooling_rate=0.001, c_repulsion=0.1)
#' print(result)
#' 
#' @export
print.topolow <- function(x, ...) {
  cat("TopoLow optimization result:\n")
  cat(sprintf("Dimensions: %d\n", x$parameters$ndim))
  cat(sprintf("Iterations: %d\n", x$iter))
  cat(sprintf("MAE: %.4f\n", x$mae))
  cat(sprintf("Correlation: %.4f\n", x$r))
  cat(sprintf("Convergence achieved: %s\n", x$convergence$achieved))
  cat(sprintf("Final convergence error: %.4f\n", x$convergence$error))
}


#' Summary method for topolow objects
#' 
#' Provides a detailed summary of the optimization results including parameters,
#' convergence trace, and performance metrics.
#' 
#' @param object A topolow object returned by topolow_full() or topolow_Smith_obj()
#' @param ... Additional arguments passed to summary (not used)
#' 
#' @examples
#' dist_mat <- matrix(c(0, 2, 3, 2, 0, 4, 3, 4, 0), nrow=3)
#' result <- topolow_full(dist_mat, ndim=2, max_iter=100, k0=1.0, cooling_rate=0.001, c_repulsion=0.1)
#' summary(result)
#' 
#' @export
summary.topolow <- function(object, ...) {
  print(object)
  cat("\nParameters:\n")
  cat(sprintf("k0: %.4f\n", object$parameters$k0))
  cat(sprintf("cooling_rate: %.4f\n", object$parameters$cooling_rate))
  cat(sprintf("c_repulsion: %.4f\n", object$parameters$c_repulsion))
  
  # Add convergence trace summary
  trace <- object$trace_convergence_error_df
  cat("\nConvergence trace summary:\n")
  cat(sprintf("Initial error: %.4f\n", trace$convergence_error[1]))
  cat(sprintf("Final error: %.4f\n", tail(trace$convergence_error, 1)))
  cat(sprintf("Error reduction: %.2f%%\n", 
              100 * (1 - tail(trace$convergence_error, 1) / trace$convergence_error[1])))
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


#' Smith variant of TopoLow algorithm
#' 
#' A variant of the TopoLow algorithm specifically designed for HI assay data, using
#' modified force calculations with sigmoid thresholding as described in Smith et al.
#' This version handles threshold measurements differently from the standard algorithm.
#'
#' @details
#' The key differences from topolow_full are:
#' * Modified force calculation for threshold measurements using sigmoid function
#' * Offset parameter for threshold calculations
#' * Specialized handling of HI assay-style measurements
#' 
#' Like Topolow algorithm, this variant is particularly suited for data where:
#' * Measurements represent binding affinities
#' * Many measurements are thresholded (< or >)
#' * True distances follow certain biological constraints
#'
#' Valid parameter ranges and constraints:
#' * ndim: Positive integer, typically 2-20. Higher dimensions increase computational cost
#' * k0: Initial spring constant, positive numeric > 0. Typical range: 0.1-30
#'      Controls initial force strength
#' * cooling_rate: Spring decay rate, numeric between 0 and 1. Typical range: 0.0001-0.1
#'          Controls how quickly spring forces weaken
#' * c_repulsion: Repulsion constant, positive numeric > 0. Typical range: 0.00001-0.2
#'      Controls repulsive force strength
#' * relative_epsilon: Positive numeric, typically 1e-6 to 1e-2
#'                    Smaller values require more iterations but give higher precision
#' * convergence_counter: Positive integer, typically 5-20
#'                     Higher values ensure more stable convergence
#' * ofs: Numeric offset parameter > 0. Typical range: 0.5-2
#'      Controls threshold sensitivity
#'
#' @inheritParams topolow_full
#' @param ofs Numeric. Offset parameter for threshold calculations. Default is 1.
#'
#' @return A list of class "topolow" containing the same elements as topolow_full() plus:
#' \itemize{
#'   \item variant: Character string "smith" indicating the algorithm variant
#'   \item parameters$ofs: The offset parameter used
#' }
#'
#' @references 
#' Smith, D. J., et al. (2004) "Mapping the Antigenic and Genetic Evolution of 
#' Influenza Virus" Science, 305(5682), 371-376.
#'
#' @seealso 
#' \code{\link{topolow_full}} for the standard algorithm version
#'
#' @examples
#' # Create a simple distance matrix with thresholds
#' dist_mat <- matrix(c(0, ">2", 3, ">2", 0, 4, 3, 4, 0), nrow=3)
#' 
#' # Run Smith variant
#' result <- topolow_Smith_obj(dist_mat, ndim=2, max_iter=1000, 
#'                            k0=1.0, cooling_rate=0.001, c_repulsion=0.1)
#'                            
#' @export
topolow_Smith_obj <- function(distance_matrix, 
                              ndim, 
                              max_iter, 
                              k0, 
                              cooling_rate, 
                              c_repulsion, 
                              relative_epsilon = 1e-4,
                              convergence_counter = 10,
                              initial_positions = NULL,
                              write_positions_to_csv = TRUE,
                              verbose = TRUE,
                              trace_sse = TRUE,
                              ofs = 1) {
  
  # Input validation
  if (!is.matrix(distance_matrix)) {
    stop("distance_matrix must be a matrix")
  }
  if (nrow(distance_matrix) != ncol(distance_matrix)) {
    stop("distance_matrix must be square")
  }
  if (!is.numeric(ndim) || ndim < 1 || ndim != round(ndim)) {
    stop("ndim must be a positive integer")  
  }
  if (ndim > 20) {
    warning("High dimensionality (ndim > 20) may lead to increased computational cost")
  }
  if (!is.numeric(max_iter) || max_iter < 1 || max_iter != round(max_iter)) {
    stop("max_iter must be a positive integer")
  }
  if (!is.numeric(ofs) || ofs <= 0) {
    stop("ofs must be a positive number")
  }
  if (ofs < 0.5 || ofs > 2) {
    warning("ofs value outside typical range [0.5, 2] may affect threshold handling")
  }
  
  # Initialize variables
  n <- nrow(distance_matrix)
  node_degrees <- rowSums(!is.na(distance_matrix))
  node_degrees_1 <- node_degrees + 1
  distances <- ifelse(is.na(distance_matrix), Inf, distance_matrix)
  
  # Initialize positions
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
  
  trace_convergence_error_df <- data.table::data.table(
    iteration = numeric(max_iter), 
    convergence_error = NULL
  )
  
  # Main optimization loop
  for(iter in 1:max_iter) {
    # Pre-compute invariant terms
    k_2 <- k / 2
    
    # Calculate forces between all pairs of points
    for(i in 1:(n-1)) {
      for(j in (i+1):n) {
        ideal_distance <- distances[i, j]
        delta <- positions[j,] - positions[i,]
        distance <- sqrt(sum(delta^2))
        distance_01 <- distance + 0.01
        
        # Spring forces:
        if(ideal_distance != Inf) {
          if (is.character(ideal_distance) && startsWith(ideal_distance, ">")) {
            threshold_value <- as.numeric(sub(">", "", ideal_distance))
            # Key difference: Modified force calculation for thresholds
            adjustment_factor <- 2*k*(threshold_value-distance+ofs)*delta/distance_01
            adjustment_factor <- adjustment_factor * gmultiple(threshold_value-distance+ofs)
          } else {
            ideal_distance <- as.numeric(ideal_distance)
            adjustment_factor <- 2*k*(ideal_distance-distance)*delta/distance_01
          }
          positions[i,] <- positions[i,] - adjustment_factor/(4*node_degrees_1[i]+k)
          positions[j,] <- positions[j,] + adjustment_factor/(4*node_degrees_1[j]+k)
        } else {
          # Repulsive force for missing measurements
          force <- c_repulsion/(2*distance_01^2)*(delta/distance_01)
          positions[i,] <- positions[i,] - force/node_degrees_1[i]
          positions[j,] <- positions[j,] + force/node_degrees_1[j]
        }
      }
    }
    
    # Calculate current distances and convergence error
    p_dist_mat <- coordinates_to_matrix(positions)
    
    # Process thresholded elements for convergence error
    distance_matrix_convergence_error <- mapply(
      process_distance_matrix, 
      distance_matrix, 
      p_dist_mat
    )
    distance_matrix_convergence_error <- matrix(
      distance_matrix_convergence_error, 
      nrow = nrow(distance_matrix)
    )
    
    # Calculate convergence error
    convergence_error_df <- data.table::data.table(
      distance_matrix_convergence_error = suppressWarnings(
        as.vector(as.numeric(distance_matrix_convergence_error))
      ),
      p_dist_mat = as.vector(p_dist_mat)
    )
    convergence_error_df <- na.omit(convergence_error_df)
    convergence_error <- mean(
      abs(convergence_error_df$distance_matrix_convergence_error - 
            convergence_error_df$p_dist_mat), 
      na.rm = TRUE
    )
    
    # Calculate evaluation metrics
    evaldf <- data.table::data.table(
      distance_matrix = suppressWarnings(
        as.vector(as.numeric(distance_matrix))
      ),
      p_dist_mat = as.vector(p_dist_mat)
    )
    evaldf <- na.omit(evaldf)
    mae <- mean(abs(evaldf$distance_matrix - evaldf$p_dist_mat), na.rm = TRUE)
    
    # Record trace
    trace_convergence_error_df[iter, "iteration"] <- iter
    trace_convergence_error_df[iter, "convergence_error"] <- convergence_error
    
    if (verbose) {
      cat(sprintf(
        "Iteration=%d, MAE=%.4f, convergence_error=%.4f\n", 
        iter, mae, convergence_error
      ))
    }
    
    # Check convergence with the same criteria as topolow_full
    if (iter > 1) {
      if (is.na(prev_convergence_error) || 
          is.na(convergence_error) || 
          convergence_error > convergence_error0 || 
          is.na((prev_convergence_error - convergence_error)/prev_convergence_error)) {
        
        if (verbose) {
          cat("! Skipping convergence check for this iteration.\n")
          cat(sprintf(
            "dim=%d, k0=%.4f, decay=%.4f, c_repulsion=%.4f, MAE=%.4f, convergence_error=%.4f\n",
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
    csv_filename <- sprintf(
      "Positions_dim_%d_k0_%.4f_decay_%.4f_c_repulsion_%.4f.csv",
      ndim, k0, cooling_rate, c_repulsion
    )
    utils::write.csv(positions, file = csv_filename, row.names = TRUE)
  }
  
  # Calculate final correlation
  pearson_corr <- stats::cor(
    evaldf$p_dist_mat, 
    evaldf$distance_matrix, 
    method = "pearson"
  )
  
  # Create result object
  result <- structure(
    list(
      positions = positions,
      est_distances = p_dist_mat,
      mae = mae,
      r = pearson_corr,
      iter = iter,
      trace_convergence_error_df = trace_convergence_error_df,
      parameters = list(
        ndim = ndim,
        k0 = k0,
        cooling_rate = cooling_rate,
        c_repulsion = c_repulsion,
        ofs = ofs
      ),
      convergence = list(
        achieved = convergence_counter <= 0,
        error = convergence_error
      ),
      variant = "smith"
    ),
    class = "topolow"
  )
  
  return(result)
}