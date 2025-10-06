# Copyright (c) 2024 Omid Arhami omid.arhami@uga.edu
# R/utils.R

#' Utility functions for the topolow package

# New
#' Extract Numeric Values from Mixed Data
#'
#' @description
#' Extracts numeric values from data that may contain threshold indicators
#' (e.g., "<10", ">1280") or regular numeric values.
#'
#' @param x A vector that may contain numeric values, character strings with
#'   threshold indicators, or a mix of both.
#' @return A numeric vector with threshold indicators converted to their
#'   numeric equivalents.
#' @examples
#' # Mixed data with threshold indicators
#' mixed_data <- c(10, 20, "<5", ">100", 50)
#' extract_numeric_values(mixed_data)
#'
#' @export
extract_numeric_values <- function(x) {
  sapply(x, function(val) {
    if (is.character(val)) {
      if (grepl("^<", val)) {
        as.numeric(sub("<", "", val))
      } else if (grepl("^>", val)) {
        as.numeric(sub(">", "", val))
      } else {
        as.numeric(val)
      }
    } else {
      as.numeric(val)
    }
  }, USE.NAMES = FALSE)
}


#' Create Cross-Validation Folds for a Dissimilarity Matrix
#'
#' @description
#' Creates k-fold cross-validation splits from a dissimilarity matrix while maintaining
#' symmetry. Each fold in the output consists of a training matrix (with some
#' values masked as `NA`) and a corresponding ground truth matrix for validation.
#'
#' @param dissimilarity_matrix The input dissimilarity matrix, which may contain noise.
#' @param ground_truth_matrix An optional, noise-free dissimilarity matrix to be used as the ground truth for evaluation. If `NULL`, the input `dissimilarity_matrix` is used as the truth.
#' @param n_folds The integer number of folds to create.
#' @param random_seed An optional integer to set the random seed for reproducibility.
#'
#' @return A list of length `n_folds`. Each element of the list is itself a list
#'   containing two matrices: `truth` (the ground truth for that fold) and `train`
#'   (the training matrix with `NA` values for validation).
#' @note This function has breaking changes from previous versions:
#'   \itemize{
#'     \item Parameter \code{truth_matrix} renamed to \code{dissimilarity_matrix}
#'     \item Parameter \code{no_noise_truth} renamed to \code{ground_truth_matrix}  
#'     \item Return structure now uses named elements (\code{$truth}, \code{$train})
#'   }
#' @examples
#' # Create a sample dissimilarity matrix
#' d_mat <- matrix(runif(100), 10, 10)
#' diag(d_mat) <- 0
#'
#' # Create 5-fold cross-validation splits
#' folds <- create_cv_folds(d_mat, n_folds = 5, random_seed = 123)
#' @export
create_cv_folds <- function(dissimilarity_matrix, ground_truth_matrix = NULL,
                            n_folds = 10, random_seed = NULL) {
  # --- Input Validation and Setup ---
  if (!is.null(random_seed)) {
    if (!is.numeric(random_seed) || random_seed != round(random_seed)) {
      stop("`random_seed` must be an integer.")
    }
    set.seed(random_seed)
  }

  if (!is.matrix(dissimilarity_matrix)) {
    stop("`dissimilarity_matrix` must be a matrix.")
  }

  if (!is.null(ground_truth_matrix)) {
    if (!is.matrix(ground_truth_matrix)) {
      stop("`ground_truth_matrix` must be NULL or a matrix.")
    }
    if (!identical(dim(dissimilarity_matrix), dim(ground_truth_matrix))) {
      stop("`dissimilarity_matrix` and `ground_truth_matrix` must have the same dimensions.")
    }
  }

  if (!is.numeric(n_folds) || n_folds < 2 || n_folds != round(n_folds)) {
    stop("`n_folds` must be an integer greater than or equal to 2.")
  }

  if (n_folds > nrow(dissimilarity_matrix)) {
    stop("`n_folds` cannot be larger than the number of rows in the matrix.")
  }

  # Use the perfect and full matrix as the ground truth for evaluation if provided
  eval_truth <- if (!is.null(ground_truth_matrix)) ground_truth_matrix else dissimilarity_matrix

  # --- Fold Creation ---
  # Initialize a list to store the [truth, train] pairs for each fold
  matrix_list <- vector("list", n_folds)
  for(i in 1:n_folds) {
    matrix_list[[i]] <- list(truth = eval_truth, train = NULL)
  }

  # Determine the number of elements to hold out in each fold
  num_elements <- sum(!is.na(dissimilarity_matrix))
  holdout_size <- floor(num_elements / (n_folds * 2)) # Factor of 2 for symmetry

  # Create a temporary copy to track available elements for sampling
  sampling_pool <- dissimilarity_matrix

  for(i in 1:n_folds) {
    # Ensure there are enough elements left to sample for the holdout set
    if (sum(!is.na(sampling_pool)) < holdout_size) {
        warning("Could not create all requested folds due to data sparsity. Returning fewer folds.")
        return(matrix_list[1:(i-1)])
    }

    # Sample validation indices from the remaining available elements
    random_indices <- sample(which(!is.na(sampling_pool)), size = holdout_size)

    # Create the training matrix for this fold by masking the holdout set
    train_matrix <- dissimilarity_matrix
    for(index in random_indices) {
      # Convert the linear index to row/column coordinates
      row <- (index - 1) %/% nrow(dissimilarity_matrix) + 1
      col <- (index - 1) %% ncol(dissimilarity_matrix) + 1

      # Mask the validation entries symmetrically
      train_matrix[row, col] <- NA
      train_matrix[col, row] <- NA

      # Remove the selected elements from the sampling pool for subsequent folds
      sampling_pool[row, col] <- NA
      # Ensure symmetry is maintained in the sampling pool as well
      sampling_pool[col, row] <- NA
    }

    matrix_list[[i]]$train <- train_matrix


  }

  return(matrix_list)
}


#' Check Dissimilarity Matrix Connectivity
#'
#' @description
#' Checks whether a dissimilarity matrix forms a connected graph, meaning all
#' points can be reached from any other point through a path of observed
#' dissimilarities. This is critical for ensuring that subsampled data will
#' allow proper embedding optimization.
#'
#' @param dissimilarity_matrix Square symmetric matrix of dissimilarities.
#'   Can contain NA values for missing measurements.
#' @param min_completeness Numeric. Minimum network completeness (fraction of
#'   possible edges that are observed) to consider acceptable. Default: 0.10.
#'   This is used as a warning threshold, not a hard requirement.
#'
#' @return A list containing connectivity diagnostics:
#'   \item{is_connected}{Logical. TRUE if the graph forms a single connected component.}
#'   \item{n_components}{Integer. Number of separate connected components (islands).}
#'   \item{completeness}{Numeric. Fraction of possible edges that are observed (0-1).}
#'   \item{n_points}{Integer. Number of points in the matrix.}
#'   \item{n_measurements}{Integer. Number of observed dissimilarities.}
#'
#' @details
#' A connected graph means there are no isolated groups of points. If the graph
#' has multiple components (islands), optimization will fail because points in
#' different islands have no observed relationships to constrain their relative
#' positions.
#'
#' The function uses \code{\link{analyze_network_structure}} to build an
#' adjacency matrix, then uses igraph to identify connected components.
#'
#' @examples
#' # Create a connected matrix
#' connected_mat <- matrix(c(0, 1, 2, 1, 0, 1.5, 2, 1.5, 0), nrow = 3)
#' result <- check_matrix_connectivity(connected_mat)
#' print(result$is_connected)  # TRUE
#'
#' # Create a disconnected matrix (two islands)
#' disconnected_mat <- matrix(NA, nrow = 4, ncol = 4)
#' diag(disconnected_mat) <- 0
#' disconnected_mat[1, 2] <- disconnected_mat[2, 1] <- 1
#' disconnected_mat[3, 4] <- disconnected_mat[4, 3] <- 1
#' result <- check_matrix_connectivity(disconnected_mat)
#' print(result$is_connected)  # FALSE
#' print(result$n_components)  # 2
#'
#' @export
check_matrix_connectivity <- function(dissimilarity_matrix,
                                      min_completeness = 0.1) {
  # Input validation
  if (!is.matrix(dissimilarity_matrix)) {
    stop("dissimilarity_matrix must be a matrix")
  }
  if (nrow(dissimilarity_matrix) != ncol(dissimilarity_matrix)) {
    stop("dissimilarity_matrix must be square")
  }
  if (nrow(dissimilarity_matrix) < 2) {
    stop("dissimilarity_matrix must have at least 2 points")
  }
  
  # Check if igraph is available
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for connectivity checking. ",
         "Please install it with: install.packages('igraph')")
  }
  
  # Get network structure using existing topolow function
  network_structure <- analyze_network_structure(dissimilarity_matrix)
  
  # Convert adjacency matrix to igraph object
  adj_matrix <- network_structure$adjacency
  graph <- igraph::graph_from_adjacency_matrix(
    adj_matrix,
    mode = "undirected"
  )
  
  # Check connectivity
  components <- igraph::components(graph)
  is_connected <- (components$no == 1)
  n_components <- components$no
  
  # Get completeness and measurement counts
  completeness <- network_structure$summary$completeness
  n_points <- network_structure$summary$n_points
  n_measurements <- network_structure$summary$n_measurements
  
  # Issue warning if completeness is low but graph is connected
  if (is_connected && completeness < min_completeness) {
    warning(sprintf(
      "Network is connected but sparse (%.1f%% complete). This may lead to poor optimization. Consider using more data points.",
      completeness * 100
    ))
  }
  
  return(list(
    is_connected = is_connected,
    n_components = n_components,
    completeness = completeness,
    n_points = n_points,
    n_measurements = n_measurements
  ))
}


#' Subsample Dissimilarity Matrix with Connectivity Guarantee
#'
#' @description
#' Randomly samples a subset of points from a dissimilarity matrix, ensuring
#' the resulting submatrix forms a connected graph. If connectivity cannot be
#' achieved, the function adaptively increases the subsample size.
#'
#' @param dissimilarity_matrix Square symmetric dissimilarity matrix.
#' @param sample_size Integer. Target number of points to sample.
#' @param max_attempts Integer. Maximum number of sampling attempts before
#'   giving up. Default: 5
#' @param min_completeness Numeric. Minimum acceptable network completeness
#'   (0-1). Default: 0.10. Used for warnings, not hard requirement.
#' @param random_seed Integer or NULL. Random seed for reproducibility. If NULL,
#'   no seed is set. Default: NULL.
#' @param verbose Logical. Print progress messages. Default: FALSE.
#'
#' @return A list containing:
#'   \item{subsampled_matrix}{Matrix. The subsampled dissimilarity matrix.}
#'   \item{selected_indices}{Integer vector. Indices of selected points from
#'     the original matrix.}
#'   \item{selected_names}{Character vector or NULL. Names of selected points
#'     (if original matrix had row/column names).}
#'   \item{is_connected}{Logical. Whether the final subsample is connected.}
#'   \item{n_components}{Integer. Number of connected components in final subsample.}
#'   \item{completeness}{Numeric. Network completeness of final subsample (0-1).}
#'   \item{attempt_number}{Integer. Which attempt succeeded.}
#'
#' @details
#' The function performs the following steps:
#' 1. Validates inputs and checks if subsampling is necessary
#' 2. Attempts to sample \code{sample_size} points randomly
#' 3. Checks if the subsample forms a connected graph using
#'    \code{\link{check_matrix_connectivity}}
#' 4. If not connected, retries
#' 5. Returns the connected subsample or stops with an error if unsuccessful
#'
#' **Why connectivity matters**: A disconnected graph has isolated groups of
#' points with no observed dissimilarities between groups. The algorithm
#' cannot determine the relative positions of disconnected components,
#' leading to arbitrary and meaningless results.
#'
#' @examples
#' # Create a well-connected matrix
#' n <- 100
#' coords <- matrix(rnorm(n * 3), ncol = 3)
#' dist_mat <- as.matrix(dist(coords))
#'
#' # Subsample to 30 points
#' result <- subsample_dissimilarity_matrix(
#'   dissimilarity_matrix = dist_mat,
#'   sample_size = 30,
#'   random_seed = 123,
#'   verbose = TRUE
#' )
#'
#' print(result$is_connected)
#' print(dim(result$subsampled_matrix))
#'
#' @seealso
#' \code{\link{check_matrix_connectivity}} for connectivity checking,
#'
#' @export
subsample_dissimilarity_matrix <- function(dissimilarity_matrix,
                                           sample_size,
                                           max_attempts = 5,
                                           min_completeness = 0.1,
                                           random_seed = NULL,
                                           verbose = FALSE) {
  # ============================================================================
  # Input Validation
  # ============================================================================
  if (!is.matrix(dissimilarity_matrix)) {
    stop("dissimilarity_matrix must be a matrix")
  }
  
  if (nrow(dissimilarity_matrix) != ncol(dissimilarity_matrix)) {
    stop("dissimilarity_matrix must be square")
  }
  
  matrix_size <- nrow(dissimilarity_matrix)
  
  if (!is.numeric(sample_size) || sample_size < 2) {
    stop("sample_size must be a numeric value >= 2")
  }
  
  sample_size <- as.integer(floor(sample_size))
  
  if (sample_size >= matrix_size) {
    if (verbose) {
      cat(sprintf("Requested sample_size (%d) >= matrix size (%d). ",
                  "Using full matrix.\n", sample_size, matrix_size))
    }
    # Return full matrix with appropriate metadata
    connectivity <- check_matrix_connectivity(dissimilarity_matrix,
                                              min_completeness)
    return(list(
      subsampled_matrix = dissimilarity_matrix,
      selected_indices = seq_len(matrix_size),
      selected_names = rownames(dissimilarity_matrix),
      is_connected = connectivity$is_connected,
      n_components = connectivity$n_components,
      completeness = connectivity$completeness,
      attempt_number = 1
    ))
  }
  
  if (!is.null(random_seed)) {
    set.seed(random_seed)
  }
  
  # ============================================================================
  # Sampling Loop
  # ============================================================================
  current_sample_size <- sample_size
  
  for (attempt in 1:max_attempts) {
    if (verbose && attempt > 1) {
      cat(sprintf("  Attempt %d/%d", attempt, max_attempts))
      cat("...\n")
    }
    
    # Sample indices
    selected_indices <- sample(seq_len(matrix_size),
                               size = current_sample_size,
                               replace = FALSE)
    
    # Extract submatrix
    subsampled_matrix <- dissimilarity_matrix[selected_indices,
                                              selected_indices,
                                              drop = FALSE]
    
    # Preserve names if they exist
    if (!is.null(rownames(dissimilarity_matrix))) {
      selected_names <- rownames(dissimilarity_matrix)[selected_indices]
      rownames(subsampled_matrix) <- selected_names
      colnames(subsampled_matrix) <- selected_names
    } else {
      selected_names <- NULL
    }
    
    # Check connectivity
    connectivity <- tryCatch({
      check_matrix_connectivity(subsampled_matrix, min_completeness)
    }, error = function(e) {
      if (verbose) {
        cat(sprintf("  Warning: Connectivity check failed: %s\n", e$message))
      }
      list(is_connected = FALSE, n_components = -1,
           completeness = 0, n_points = current_sample_size,
           n_measurements = 0)
    })
    
    # Success - return connected subsample
    if (connectivity$is_connected) {
      if (verbose) {
        cat(sprintf("  [OK] Connected subsample obtained (attempt %d, size %d, %.1f%% complete)\n",
                    attempt, current_sample_size,
                    connectivity$completeness * 100))
      }
      
      return(list(
        subsampled_matrix = subsampled_matrix,
        selected_indices = selected_indices,
        selected_names = selected_names,
        is_connected = TRUE,
        n_components = connectivity$n_components,
        completeness = connectivity$completeness,
        attempt_number = attempt
      ))
    }
    
    # Failed this attempt
    if (verbose) {
      cat(sprintf("  X Not connected (%d components, %.1f%% complete)\n",
                  connectivity$n_components,
                  connectivity$completeness * 100))
    }
  }
  
  # ============================================================================
  # Failed to find connected subsample
  # ============================================================================
  stop(sprintf(
    paste0(
      "Failed to obtain a connected subsample after %d attempts.\n",
      "  Final sample size tried: %d (started at %d)\n",
      "  Original matrix size: %d\n",
      "  Last attempt had %d components with %.1f%% completeness\n\n",
      "Possible solutions:\n",
      "  1. Increase opt_subsample (current: %d)\n",
      "  2. Reduce number of CV folds\n",
      "  3. Use full dataset (opt_subsample = NULL)\n",
      "  4. Check if your data has inherent disconnected groups"
    ),
    max_attempts, current_sample_size, sample_size, matrix_size,
    connectivity$n_components, connectivity$completeness * 100,
    sample_size
  ))
}


#' Sanity Check for Subsampled Data Before Optimization
#'
#' @description
#' Validates that a subsampled dissimilarity matrix has sufficient data for
#' reliable parameter optimization with cross-validation. We need enough ground
#' truth in each fold to be able to compare the predictions with them with high
#' confidence.
#'
#' @param subsampled_matrix Matrix. The subsampled dissimilarity matrix to check.
#' @param folds Integer. Number of cross-validation folds that will be used.
#' @param min_points_per_fold Integer. Minimum number of data points expected
#'   per fold. Default: 5.
#' @param min_measurements_per_fold Integer. Minimum number of measurements
#'   per fold. Default: 10.
#' @param verbose Logical. Print detailed diagnostic messages. Default: TRUE.
#'
#' @return A list with check results:
#'   \item{all_checks_passed}{Logical. TRUE if all checks passed.}
#'   \item{checks}{List of individual check results (logical values).}
#'   \item{diagnostics}{List of diagnostic values calculated.}
#'   \item{warnings}{Character vector of warning messages (empty if no warnings).}
#'
#' @details
#' The function performs these checks:
#' \itemize{
#'   \item \strong{Sufficient points}: At least 2x folds points in total
#'   \item \strong{Sufficient measurements}: At least folds x min_measurements_per_fold
#'   \item \strong{Adequate points per fold}: Average points per fold >= min_points_per_fold
#'   \item \strong{Adequate measurements per fold}: Average measurements per fold >= min_measurements_per_fold
#'   \item \strong{Not too sparse}: Overall sparsity < 95%
#' }
#'
#' If checks fail, warnings are issued but execution continues. This allows the
#' user to proceed with awareness of potential issues.
#'
#' @examples
#' # Create a small matrix
#' small_mat <- matrix(c(0, 1, 2, 1, 0, 1.5, 2, 1.5, 0), nrow = 3)
#'
#' # Check with 5 folds (will fail - too few points)
#' result <- sanity_check_subsample(small_mat, folds = 5)
#' print(result$all_checks_passed)
#' print(result$warnings)
#'
#' # Create a larger connected matrix
#' n <- 50
#' coords <- matrix(rnorm(n * 2), ncol = 2)
#' large_mat <- as.matrix(dist(coords))
#'
#' # Check with 10 folds (should pass)
#' result <- sanity_check_subsample(large_mat, folds = 10)
#' print(result$all_checks_passed)
#'
#' @export
sanity_check_subsample <- function(subsampled_matrix,
                                   folds = 20,
                                   min_points_per_fold = 3,
                                   min_measurements_per_fold = 3,
                                   verbose = TRUE) {
  # Calculate diagnostics
  n_points <- nrow(subsampled_matrix)
  n_measurements <- as.integer(sum(!is.na(subsampled_matrix)) / 2) # Divide by 2 for symmetry
  total_possible <- n_points * (n_points - 1) / 2
  sparsity <- 1 - (n_measurements / total_possible)
  
  # Each fold will have approximately this many measurements
  # In k-fold CV, roughly 1/k of data is held out per fold
  est_measurements_per_fold <- n_measurements / folds
  est_points_per_fold <- n_points  # All points are used, some connections hidden
  
  # Initialize checks
  checks <- list()
  warnings_list <- character(0)
  
  # Check 1: Sufficient total points (at least 2x folds for meaningful CV)
  checks$sufficient_points <- n_points >= (2 * folds)
  if (!checks$sufficient_points) {
    msg <- sprintf(
      "Very few points (%d) for %d-fold CV. Consider reducing folds or increasing subsample size.",
      n_points, folds
    )
    warnings_list <- c(warnings_list, msg)
    if (verbose) warning(msg, call. = FALSE)
  }
  
  # Check 2: Sufficient total measurements
  min_total_measurements <- folds * min_measurements_per_fold
  checks$sufficient_measurements <- n_measurements >= min_total_measurements
  if (!checks$sufficient_measurements) {
    msg <- sprintf(
      "Insufficient measurements (%d) for %d-fold CV. Expected at least %d.",
      n_measurements, folds, min_total_measurements
    )
    warnings_list <- c(warnings_list, msg)
    if (verbose) warning(msg, call. = FALSE)
  }
  
  # Check 3: Adequate measurements per fold
  checks$adequate_measurements_per_fold <- est_measurements_per_fold >= min_measurements_per_fold
  if (!checks$adequate_measurements_per_fold) {
    msg <- sprintf(
      "Only ~%.1f measurements per fold (expected >= %d). Results may be unreliable.",
      est_measurements_per_fold, min_measurements_per_fold
    )
    warnings_list <- c(warnings_list, msg)
    if (verbose) warning(msg, call. = FALSE)
  }
  
  # Check 4: Not extremely sparse
  checks$not_too_sparse <- sparsity < 0.95
  if (!checks$not_too_sparse) {
    msg <- sprintf(
      "Matrix is %.1f%% sparse. Such extreme sparsity may cause optimization issues.",
      sparsity * 100
    )
    warnings_list <- c(warnings_list, msg)
    if (verbose) warning(msg, call. = FALSE)
  }
  
  # Check 5: Has any measurements at all
  checks$has_measurements <- n_measurements > 0
  if (!checks$has_measurements) {
    msg <- "No measurements found in subsampled matrix!"
    warnings_list <- c(warnings_list, msg)
    if (verbose) warning(msg, call. = FALSE)
  }
  
  # Overall assessment
  all_checks_passed <- all(unlist(checks))
  
  # Compile diagnostics
  diagnostics <- list(
    n_points = n_points,
    n_measurements = n_measurements,
    sparsity = sparsity,
    folds = folds,
    est_measurements_per_fold = est_measurements_per_fold,
    est_points_per_fold = est_points_per_fold
  )
  
  if (verbose && all_checks_passed) {
    cat(sprintf(
      "[OK] Subsample sanity checks passed: %d points, %d measurements (%.1f%% sparse), %d folds\n",
      n_points, n_measurements, sparsity * 100, folds
    ))
  }
  
  return(list(
    all_checks_passed = all_checks_passed,
    checks = checks,
    diagnostics = diagnostics,
    warnings = warnings_list
  ))
}


#' Save ggplot with white background
#'
#' @description
#' Wrapper around `ggplot2::ggsave` that ensures a white background by default.
#' @importFrom ggplot2 ggsave
#' @inheritParams ggplot2::ggsave
#' @return No return value, called for side effects.
#' @export
ggsave_white_bg <- function(..., bg = 'white') {
  ggplot2::ggsave(..., bg = bg)
}


#' Color Palettes
#'
#' @description
#' Predefined color palettes optimized for visualization.
#'
#' @name color_palettes
NULL

#' @rdname color_palettes
#' @export
c25 <- c(
  "purple", "green1", "blue1", "gold1", "red",
  "darkturquoise", "darkorange", "skyblue2", "green4", "maroon",
  "yellow3", "gray40", "hotpink", "darkorange4", "deeppink1",
  "khaki2", "palegreen2", "dodgerblue2", "brown", "orchid1"
)
