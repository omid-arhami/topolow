# Copyright (c) 2024 Omid Arhami omid.arhami@uga.edu
# R/diagnostics.R

#' Model Diagnostics and Convergence Testing

# Newed:
#' Check Multivariate Gaussian Convergence
#'
#' @description
#' Assesses the convergence of multivariate samples by monitoring the stability of the
#' mean vector and covariance matrix over a sliding window. This is useful for checking
#' if a set of parameter samples has stabilized.
#'
#' @param data Matrix or Data Frame. A matrix of samples where columns are parameters.
#' @param window_size Integer. The size of the sliding window used to compute statistics.
#' @param tolerance Numeric. The convergence threshold for the relative change in the
#'        mean and covariance.
#' @return An object of class `topolow_convergence` containing diagnostics about the
#'   convergence of the multivariate samples. This list includes logical flags for
#'   convergence (`converged`, `mean_converged`, `cov_converged`) and the history
#'   of the mean and covariance changes.
#' @examples
#' # Create sample data for the example
#' chain_data <- as.data.frame(matrix(rnorm(500 * 4), ncol = 4))
#' colnames(chain_data) <- c("param1", "param2", "param3", "param4")
#'
#' # Run the convergence check
#' conv_results <- check_gaussian_convergence(chain_data)
#' print(conv_results)
#'
#' # The plot method for this object can be used to create convergence plots.
#' # plot(conv_results)
#' @export
check_gaussian_convergence <- function(data, window_size = 300, tolerance = 0.01) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("data must be a matrix or data frame")
  }

  data[] <- lapply(data, as.numeric)
  data <- stats::na.omit(data)

  if(is.data.frame(data)) {
    data <- as.matrix(data)
  }

  n_samples <- nrow(data)
  n_dims <- ncol(data)

  if (window_size >= n_samples) {
    stop("window_size must be less than the number of samples")
  }
  if (!is.numeric(tolerance) || tolerance <= 0) {
    stop("tolerance must be a positive number")
  }
  if (!all(apply(data, 2, is.numeric))) {
    stop("All columns must be numeric")
  }

  # Initialize storage for tracking statistics
  mean_history <- matrix(NA, nrow = n_samples, ncol = n_dims)
  cov_norm_history <- numeric(n_samples)

  # Calculate running statistics
  for(i in window_size:n_samples) {
    current_data <- data[1:i,]
    mean_history[i,] <- colMeans(current_data)
    cov_norm_history[i] <- norm(stats::cov(current_data), "F")  # Frobenius norm
  }

  # Fill in earlier values to avoid NA warnings in plots
  if (window_size > 1) {
    for(i in 1:(window_size-1)) {
      mean_history[i,] <- mean_history[window_size,]
    }
  }

  # Calculate relative changes in mean (Euclidean distance)
  mean_changes <- sapply((window_size + 1):n_samples, function(i) {
      sqrt(sum(((mean_history[i,] - mean_history[i-1,]) / mean_history[i-1,])^2))
  })

  # Calculate relative changes in covariance norm
  cov_changes <- abs(cov_norm_history[(window_size + 1):n_samples] - cov_norm_history[window_size:(n_samples - 1)]) /
                 cov_norm_history[window_size:(n_samples - 1)]


  # Check convergence based on the last window
  mean_converged <- all(utils::tail(mean_changes, window_size) < tolerance, na.rm = TRUE)
  cov_converged <- all(utils::tail(cov_changes, window_size) < tolerance, na.rm = TRUE)

  # Return results with the new S3 class
  structure(list(
    converged = mean_converged && cov_converged,
    mean_converged = mean_converged,
    cov_converged = cov_converged,
    final_mean = colMeans(utils::tail(data, window_size)),
    final_cov = stats::cov(utils::tail(data, window_size)),
    mean_history = mean_history,
    cov_changes = cov_changes,
    param_names = colnames(data)
  ), class = "topolow_convergence")
}


# Newed:
#' Calculate MCMC-style Diagnostics for Sampling Chains
#'
#' @description
#' Calculates standard MCMC-style convergence diagnostics for multiple chains from an
#' optimization or sampling run. It computes the R-hat (potential scale reduction factor)
#' and effective sample size (ESS) to help assess if the chains have converged to a
#' stable distribution.
#'
#' @param chain_files Character vector. Paths to CSV files, where each file represents a chain of samples.
#' @param mutual_size Integer. Number of samples to use from the end of each chain for calculations.
#' @return A list object of class `topolow_diagnostics` containing convergence diagnostics for the MCMC chains.
#'   \item{rhat}{A numeric vector of the R-hat (potential scale reduction factor) statistic for each parameter. Values close to 1 indicate convergence.}
#'   \item{ess}{A numeric vector of the effective sample size for each parameter.}
#'   \item{chains}{A list of data frames, where each data frame is a cleaned and trimmed MCMC chain.}
#'   \item{param_names}{A character vector of the parameter names being analyzed.}
#'   \item{mutual_size}{The integer number of samples used from the end of each chain for calculations.}
#' @examples
#' # This example demonstrates how to use the function with temporary files.
#' # Create dummy chain files in a temporary directory
#' temp_dir <- tempdir()
#' chain_files <- character(3)
#' par_names <- c("log_N", "log_k0", "log_cooling_rate", "log_c_repulsion")
#' sample_data <- data.frame(
#'   log_N = rnorm(100), log_k0 = rnorm(100),
#'   log_cooling_rate = rnorm(100), log_c_repulsion = rnorm(100),
#'   NLL = runif(100), Holdout_MAE = runif(100)
#' )
#' for (i in 1:3) {
#'   chain_files[i] <- file.path(temp_dir, paste0("chain", i, ".csv"))
#'   write.csv(sample_data, chain_files[i], row.names = FALSE)
#' }
#'
#' # Calculate diagnostics
#' diag_results <- calculate_diagnostics(chain_files, mutual_size = 50)
#' print(diag_results)
#'
#' # Clean up the temporary files and directory
#' unlink(chain_files)
#' unlink(temp_dir, recursive = TRUE)
#' @export
calculate_diagnostics <- function(chain_files, mutual_size=500) {
  # Check if coda package is available
  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("coda package is required for diagnostics. Please install with install.packages('coda').")
  }
  # Validate inputs
  if (!is.character(chain_files)) {
    stop("chain_files must be a character vector")
  }
  if (!all(file.exists(chain_files))) {
    missing <- chain_files[!file.exists(chain_files)]
    stop("Missing chain files: ", paste(missing, collapse = ", "))
  }
  if (!is.numeric(mutual_size) || mutual_size <= 0 || mutual_size != round(mutual_size)) {
    stop("mutual_size must be a positive integer")
  }

  # Define parameter names to analyze
  par_names <- c("log_N", "log_k0", "log_cooling_rate", "log_c_repulsion")

  # Read and clean chains
  chains <- lapply(chain_files, function(file) {
    df <- utils::read.csv(file, na.strings = c("NA", ""), stringsAsFactors = FALSE)
    required_cols <- c(par_names, "NLL", "Holdout_MAE")
    if(!all(required_cols %in% names(df))) {
        stop("File ", file, " missing required columns.")
    }
    for(col in required_cols) df[[col]] <- as.numeric(df[[col]])

    # Filter out rows with NAs or Infs in any required column
    valid_rows <- complete.cases(df[, required_cols]) &
                  apply(df[, required_cols], 1,
                        function(row) all(is.finite(row)))

    clean_df <- df[valid_rows, ]

    # Apply clean_data to all columns of the dataframe
    clean_df <- as.data.frame(lapply(clean_df, clean_data, k = 3))
    # Remove rows with NAs introduced by clean_data
    clean_df <- na.omit(clean_df)
    return(clean_df[, par_names, drop=FALSE])
  })

  # Check if chains have enough rows after cleaning
  if(any(sapply(chains, nrow) < mutual_size)) {
    stop("One or more chains have fewer valid rows than mutual_size.")
  }

  # Extract the last `mutual_size` rows from each chain
  chains <- lapply(chains, function(chain) utils::tail(chain, mutual_size))

  # Convert to mcmc.list for coda diagnostics
  mcmc_list <- tryCatch({
    coda::mcmc.list(lapply(chains, coda::mcmc))
  }, error = function(e) {
    warning("Error in Gelman diagnostics: ", e$message)
      # Return NA values if calculation fails
      matrix(NA, nrow=length(par_names), ncol=2,
             dimnames=list(par_names, c("Point est.", "Upper C.I.")))
  })

  # Calculate Gelman-Rubin diagnostic (R-hat)
  gelman_result <- tryCatch(coda::gelman.diag(mcmc_list), error = function(e) {
    warning("Could not calculate Gelman-Rubin diagnostic: ", e$message); NULL
  })
  rhat <- if (!is.null(gelman_result)) gelman_result$psrf[,1] else rep(NA, length(par_names))

  # Calculate Effective Sample Size (ESS)
  ess <- tryCatch(coda::effectiveSize(mcmc_list), error = function(e) {
    warning("Could not calculate Effective Sample Size: ", e$message); rep(NA, length(par_names))
  })

  # Add names to avoid missing value issues:
  names(rhat) <- par_names
  names(ess) <- par_names

  # Return results with the new S3 class
  structure(list(
    rhat = rhat,
    ess = ess,
    chains = chains,
    param_names = par_names,
    mutual_size = mutual_size
  ), class = "topolow_diagnostics")
}


# Newed:
#' Analyze Network Structure
#'
#' @description
#' Analyzes the connectivity of a dissimilarity matrix, returning node degrees and
#' overall completeness.
#'
#' @param dissimilarity_matrix Square symmetric matrix of dissimilarities.
#' @return A list containing the network analysis results:
#'   \item{adjacency}{A logical `matrix` where `TRUE` indicates a measured dissimilarity.}
#'   \item{connectivity}{A `data.frame` with node-level metrics, including the `completeness` (degree) for each point.}
#'   \item{summary}{A list of overall network statistics, including `n_points`, `n_measurements`, and total `completeness`.}
#' @examples
#' # Create a sample dissimilarity matrix
#' dist_mat <- matrix(runif(25), 5, 5)
#' rownames(dist_mat) <- colnames(dist_mat) <- paste0("Point", 1:5)
#' dist_mat[lower.tri(dist_mat)] <- t(dist_mat)[lower.tri(dist_mat)]
#' diag(dist_mat) <- 0
#' dist_mat[1, 3] <- NA; dist_mat[3, 1] <- NA
#'
#' # Analyze the network structure
#' metrics <- analyze_network_structure(dist_mat)
#' print(metrics$summary$completeness)
#' @export
analyze_network_structure <- function(dissimilarity_matrix) {
  if (!is.matrix(dissimilarity_matrix) || nrow(dissimilarity_matrix) != ncol(dissimilarity_matrix)) {
    stop("Input must be a square matrix")
  }
  if (nrow(dissimilarity_matrix) < 2) {
    stop("Matrix must have at least 2 rows/columns")
  }

  # Create adjacency matrix (TRUE where a measurement exists)
  adjacency <- !is.na(dissimilarity_matrix)
  diag(adjacency) <- FALSE # No self-loops

  # Get or create node names
  node_names <- rownames(dissimilarity_matrix)
  if (is.null(node_names)) {
    node_names <- paste0("Node", seq_len(nrow(dissimilarity_matrix)))
  }

  # Set the rownames and colnames of the adjacency matrix
  rownames(adjacency) <- node_names
  colnames(adjacency) <- node_names

  # Calculate node-level metrics
  connectivity <- data.frame(
    node = node_names,
    degree = rowSums(adjacency),
    completeness = rowSums(adjacency) / (ncol(dissimilarity_matrix) - 1)
  )

  # Calculate overall network summary statistics
  summary_stats <- list(
    n_points = nrow(dissimilarity_matrix),
    n_measurements = sum(adjacency) / 2, # Each edge counted twice in sum
    completeness = sum(adjacency) / (nrow(dissimilarity_matrix) * (nrow(dissimilarity_matrix)-1))
  )

  return(list(
    adjacency = adjacency,
    connectivity = connectivity,
    summary = summary_stats
  ))
}

