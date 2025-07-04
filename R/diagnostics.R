# Copyright (c) 2024 Omid Arhami omid.arhami@uga.edu
# R/diagnostics.R

#' Model Diagnostics and Convergence Testing
#' 
#' @description
#' This file contains functions for assessing model convergence, analyzing chains,
#' and performing diagnostic tests. Functions are designed to be general-purpose
#' and usable with any iterative optimization or sampling procedure.
#'
#' Functions handle:
#' - Convergence testing
#' - Chain analysis
#' - Statistical diagnostics
#' - Parameter distribution analysis
#'
#' @importFrom coda gelman.diag effectiveSize mcmc mcmc.list
#' @importFrom stats var sd quantile
#' @importFrom parallel mclapply detectCores
#' @keywords internal
"_PACKAGE"


#' Check Multivariate Gaussian Convergence
#'
#' @description
#' Assesses convergence of multivariate samples by monitoring changes in mean 
#' vector and covariance matrix over a sliding window. Useful for checking
#' stability of parameter distributions in optimization or sampling.
#'
#' @param data Matrix or data frame of samples where columns are parameters
#' @param window_size Integer size of sliding window for statistics
#' @param tolerance Numeric convergence threshold for relative changes
#' @return An object of class `topolow_convergence` containing diagnostics about the convergence of the multivariate samples. This list includes:
#'   \item{converged}{A logical flag, `TRUE` if both mean and covariance have converged.}
#'   \item{mean_converged}{A logical flag, `TRUE` if the mean vector has converged.}
#'   \item{cov_converged}{A logical flag, `TRUE` if the covariance matrix has converged.}
#'   \item{final_mean}{The mean vector calculated from the last `window_size` samples.}
#'   \item{final_cov}{The covariance matrix calculated from the last `window_size` samples.}
#'   \item{mean_history}{A matrix tracking the history of the running mean of each parameter.}
#'   \item{cov_changes}{A numeric vector of the relative changes in the Frobenius norm of the covariance matrix.}
#'   \item{param_names}{The names of the parameters (columns) from the input data.}
#' @examples
#' # Assuming 'chain_data' is a data frame with samples
#' chain_data <- as.data.frame(matrix(rnorm(500 * 4), ncol = 4))
#' colnames(chain_data) <- c("log_N", "log_k0", "log_cooling_rate", "log_c_repulsion")
#' conv_results <- check_gaussian_convergence(chain_data)
#' print(conv_results)  # Shows summary
#' # The plot method for this object would create convergence plots.
#' # plot(conv_results)
#' @importFrom stats na.omit cov
#' @importFrom utils tail
#' @export
check_gaussian_convergence <- function(data, window_size = 300, tolerance = 1e-2) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("data must be a matrix or data frame")
  }
  
  data[] <- lapply(data, as.numeric)
  data <- na.omit(data)
  
  # Convert to matrix if needed
  if(is.data.frame(data)) {
    data <- as.matrix(data)
  }

  n_samples <- nrow(data)
  n_dims <- ncol(data)
  
  if (window_size >= n_samples) {
    stop("window_size must be less than number of samples")
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
    cov_norm_history[i] <- norm(cov(current_data), "F")  # Frobenius norm
  }
  
  # Calculate relative changes in mean
  mean_changes <- numeric(n_samples - 1)
  for(i in (window_size + 1):n_samples) {
    # Calculate relative Euclidean distance between consecutive means
    mean_changes[i-1] <- sqrt(sum(((mean_history[i,] - mean_history[i-1,])/
                                    mean_history[i-1,])^2))
  }
  
  # Calculate relative changes in covariance norm
  cov_changes <- numeric(n_samples - 1)
  for(i in (window_size + 1):n_samples) {
    cov_changes[i-1] <- abs(cov_norm_history[i] - cov_norm_history[i-1])/
      cov_norm_history[i-1]
  }
  
  # Check convergence
  recent_mean <- colMeans(tail(mean_history, window_size))
  recent_cov <- mean(tail(cov_norm_history, window_size))
  
  recent_mean_changes <- tail(mean_changes, window_size)
  recent_cov_changes <- tail(cov_changes, window_size)
  
  mean_converged <- all(recent_mean_changes < tolerance, na.rm = TRUE)
  cov_converged <- all(recent_cov_changes < tolerance, na.rm = TRUE)
  
  # Return results with class
  structure(list(
    converged = mean_converged && cov_converged,
    mean_converged = mean_converged,
    cov_converged = cov_converged,
    final_mean = recent_mean,
    final_cov = cov(tail(data, window_size)),
    mean_history = mean_history,
    cov_changes = cov_changes,
    param_names = colnames(data)
  ), class = "topolow_convergence")
}

#' Calculate Adaptive Monte Carlo Sampling Diagnostics
#'
#' @description
#' Calculates standard Adaptive Monte Carlo Sampling diagnostics including R-hat (potential scale reduction)
#' and effective sample size for multiple chains. Can be used with any iterative
#' sampling or optimization procedure that produces chain-like output.
#'
#' @param chain_files Character vector of paths to CSV files containing chains
#' @param mutual_size Integer number of samples to use from end of each chain
#' @return A list object of class `topolow_amcs_diagnostics` containing convergence diagnostics for the MCMC chains.
#'   \item{rhat}{A numeric vector of the R-hat (potential scale reduction factor) statistic for each parameter. Values close to 1 indicate convergence.}
#'   \item{ess}{A numeric vector of the effective sample size for each parameter.}
#'   \item{chains}{A list of data frames, where each data frame is a cleaned and trimmed MCMC chain.}
#'   \item{param_names}{A character vector of the parameter names being analyzed.}
#'   \item{mutual_size}{The integer number of samples used from the end of each chain for calculations.}
#' @examples
#' # This example demonstrates how to use the function with temporary files,
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
#' @importFrom utils read.csv
#' @importFrom stats complete.cases var
#' @importFrom coda mcmc.list mcmc gelman.diag effectiveSize
#' @export
calculate_diagnostics <- function(chain_files, mutual_size=500) {
  # Validate inputs
  if (!is.character(chain_files)) {
    stop("chain_files must be a character vector")
  }
  
  if (!all(file.exists(chain_files))) {
    missing <- chain_files[!file.exists(chain_files)]
    stop("Missing chain files: ", paste(missing, collapse = ", "))
  }
  
  if (!is.numeric(mutual_size) || mutual_size <= 0 ||
      mutual_size != round(mutual_size)) {
    stop("mutual_size must be a positive integer")
  }
  
  # Get parameter names we need to analyze
  par_names <- c("log_N", "log_k0", "log_cooling_rate", "log_c_repulsion")
  
  # Read and clean chains
  chains <- lapply(chain_files, function(file) {
    # Read the file with special handling for empty strings
    df <- read.csv(file, na.strings = c("NA", ""), stringsAsFactors = FALSE)
    
    # Check for required columns
    missing_cols <- setdiff(c(par_names, "NLL", "Holdout_MAE"), names(df))
    if(length(missing_cols) > 0) {
      stop("File ", file, " missing required columns: ", 
           paste(missing_cols, collapse = ", "))
    }
    
    # Ensure all columns are numeric
    for(col in c(par_names, "NLL", "Holdout_MAE")) {
      if(!is.numeric(df[[col]])) {
        df[[col]] <- as.numeric(df[[col]])
      }
    }
    
    # Filter out rows with NAs or Infs in any required column
    valid_rows <- complete.cases(df[, c(par_names, "NLL", "Holdout_MAE")]) & 
                  apply(df[, c(par_names, "NLL", "Holdout_MAE")], 1, 
                        function(row) all(is.finite(row)))
    
    clean_df <- df[valid_rows, ]
    # Apply clean_data to all columns of the dataframe
    clean_df <- as.data.frame(lapply(clean_df, clean_data, k = 3))
    
    # Return only parameter columns
    return(clean_df[, par_names, drop=FALSE])
  })
  
  # Check if chains have enough rows after cleaning
  chain_sizes <- sapply(chains, nrow)
  if(any(chain_sizes < mutual_size)) {
    problem_chains <- which(chain_sizes < mutual_size)
    sizes_info <- paste(problem_chains, ":", chain_sizes[problem_chains], collapse=", ")
    stop("Some chains have fewer valid rows than mutual_size (", mutual_size, 
         "). Chain sizes: ", sizes_info)
  }
  
  # Extract the last mutual_size rows from each chain
  for(i in seq_along(chains)) {
    n_rows <- nrow(chains[[i]])
    chains[[i]] <- chains[[i]][(n_rows - mutual_size + 1):n_rows, , drop=FALSE]
  }
  
  # Check for zero variance parameters that could cause problems
  for(i in seq_along(chains)) {
    variances <- apply(chains[[i]], 2, var)
    zero_var_params <- which(variances < 1e-10)
    if(length(zero_var_params) > 0) {
      warning("Chain ", i, " has near-zero variance in parameters: ", 
              paste(names(zero_var_params), collapse=", "), 
              ". This may cause problems with diagnostics.")
    }
  }
  
  # Convert to mcmc.list with error handling
  tryCatch({
    mcmc_list <- coda::mcmc.list(lapply(chains, function(chain) {
      coda::mcmc(as.matrix(chain))
    }))
    
    # Calculate diagnostics
    gelman_result <- tryCatch({
      coda::gelman.diag(mcmc_list)
    }, error = function(e) {
      warning("Error in Gelman diagnostics: ", e$message)
      # Return NA values if calculation fails
      matrix(NA, nrow=length(par_names), ncol=2, 
             dimnames=list(par_names, c("Point est.", "Upper C.I.")))
    })
    
    # Get R-hat values (if available)
    rhat <- if(is.null(gelman_result$psrf)) {
      rep(NA, length(par_names))
    } else {
      gelman_result$psrf[,1]
    }
    
    # Calculate effective sample size
    ess <- tryCatch({
      coda::effectiveSize(mcmc_list)
    }, error = function(e) {
      warning("Error calculating effective sample size: ", e$message)
      rep(NA, length(par_names))
    })
    
    # Return results
    structure(list(
      rhat = rhat,
      ess = ess,
      chains = chains,
      param_names = par_names,
      mutual_size = mutual_size
    ), class = "topolow_amcs_diagnostics")
    
  }, error = function(e) {
    stop("Error in MCMC diagnostics: ", e$message)
  })
}



#' Calculate Cumulative Distance Metrics
#'
#' @description
#' Calculates cumulative distance metrics either from a reference point or between
#' all pairs. Handles both seasonal and year-based analyses.
#'
#' @param df_coords Data frame containing:
#'        - V1...Vn coordinate columns
#'        - year: Numeric years
#'        - season: Character season identifiers.
#'        - cluster: Factor cluster assignments
#'        - color: Character color codes
#' @param ndim Number of coordinate dimensions
#' @param reference_row Integer index of reference row (or FALSE for all-pairs analysis)
#' @param na.rm Logical indicating whether to remove NA values
#' @return A list containing the calculated distance metrics. The content of the list depends 
#' on the `reference_row` parameter.
#' \itemize{
#'   \item If `reference_row` is specified, the list contains `summary_data`: a `data.frame` 
#'    with distances from the specified reference point to all other points, summarized by 
#'    season and cluster. The columns include `season_num`, `cluster`, `color`, `avg_euclidean_dist`, 
#'    `count`, `total_count`, and `fraction`.
#'   \item If `reference_row` is `FALSE`, the list contains `dist_data`: a `data.frame` with all unique 
#'    pairwise distances. The columns include `year_diff`, `euclidean_dist`, and `ref_year`.
#' }
#' @examples
#' # Create sample coordinate data
#' coords <- data.frame(V1 = rnorm(10), V2 = rnorm(10), year = rep(2000:2004, 2),
#'                      season = paste0(rep(2000:2004, 2), "-S"),
#'                      cluster = factor(rep(1:2, 5)), color = "blue")
#' 
#' # Calculate distances from reference point
#' ref_distances <- calculate_cumulative_distances(coords, ndim=2, reference_row=1)
#'
#' # Calculate all pairwise distances
#' all_distances <- calculate_cumulative_distances(coords, ndim=2, reference_row=FALSE)
#'
#' @importFrom stats na.omit
#' @importFrom dplyr %>% filter mutate group_by summarize
#' @importFrom grDevices colors
#' @export
calculate_cumulative_distances <- function(df_coords, ndim, 
                                        reference_row=FALSE,
                                        na.rm=TRUE) {
  if (!is.data.frame(df_coords)) {
    stop("df_coords must be a data frame")
  }
  
  # Validate coordinate columns
  coord_cols <- paste0("V", 1:ndim)
  if (!all(coord_cols %in% names(df_coords))) {
    stop("Missing coordinate columns V1 through V", ndim)
  }
  
  if (na.rm) {
    df_coords <- na.omit(df_coords)
  }
  
  if (reference_row) {
    if (!all(c("season", "cluster", "color") %in% names(df_coords))) {
      stop("Missing required columns: season, cluster, and color")
    }
    
    # Calculate reference-based distances
    dist_data <- expand.grid(
      sample1 = seq_len(nrow(df_coords)),
      sample2 = seq_len(nrow(df_coords))
    ) %>%
      filter(sample1 == reference_row) %>%
      mutate(
        euclidean_dist = sqrt(rowSums(
          (df_coords[sample1, coord_cols] - df_coords[sample2, coord_cols])^2)),
        year = df_coords$year[sample2],
        season = df_coords$season[sample2],
        season_num = as.numeric(sub("-.*", "", season)) + 0.5,
        ref_year = df_coords$year[sample1],
        cluster = df_coords$cluster[sample2],
        color = df_coords$color[sample2]
      )
    
    # Filter valid colors
    valid_colors <- colors()
    dist_data <- dist_data %>%
      filter(color %in% valid_colors)
    
    # Calculate summary statistics
    summary_data <- dist_data %>%
      group_by(season_num, cluster, color) %>%
      summarize(avg_euclidean_dist = mean(euclidean_dist),
                count = dplyr::n(), .groups = 'drop') %>%
      group_by(season_num) %>% 
      mutate(total_count = sum(count),
             fraction = count / total_count)
    
    return(list(summary_data = summary_data))
    
  } else {
    # Calculate all pairwise distances
    dist_data <- expand.grid(
      sample1 = seq_len(nrow(df_coords)),
      sample2 = seq_len(nrow(df_coords))
    ) %>%
      filter(sample1 < sample2) %>% 
      mutate(
        euclidean_dist = sqrt(rowSums(
          (df_coords[sample1, coord_cols] - df_coords[sample2, coord_cols])^2)),
        year_diff = abs(
          df_coords$year[sample2] - df_coords$year[sample1]),
        ref_year = df_coords$year[sample1]
      )
    
    return(list(dist_data = dist_data))
  }
}

#' Calculate Annual Distance Metrics
#'
#' @description
#' Calculates year-over-year antigenic distances and statistics. Compares each point
#' to the mean coordinates of the previous year.
#'
#' @param df_coords Data frame containing:
#'        - V1...Vn coordinate columns
#'        - year: Numeric years
#'        - name: Point identifiers (will use rownames if missing)
#' @param ndim Number of coordinate dimensions
#' @param na.rm Logical indicating whether to remove NA values
#' @return A list containing year-over-year antigenic distance metrics:
#'   \item{dist_data}{A `data.frame` where each row represents a point and its `distance` 
#'    to the mean coordinate of the previous year.}
#'   \item{summary}{A list containing the `overall_mean` and `overall_sd` (standard deviation) 
#'    of the annual distances across all years.}
#' @examples
#' # Create sample coordinate data
#' coords <- data.frame(V1 = rnorm(10), V2 = rnorm(10), year = rep(2000:2004, 2))
#' annual_stats <- calculate_annual_distances(coords, ndim=2)
#' print(annual_stats$summary$overall_mean)
#'
#' @importFrom stats na.omit sd
#' @importFrom dplyr %>% mutate left_join group_by summarise across all_of select
#' @export
calculate_annual_distances <- function(df_coords, ndim, na.rm=TRUE) {
  # Add name column if missing
  if (!("name" %in% names(df_coords))) {
    df_coords$name <- rownames(df_coords)
  }
  
  # Clean point identifiers
  df_coords$name <- sub("^(S/|V/)", "", df_coords$name)
  
  if (na.rm) {
    df_coords <- na.omit(df_coords)
  }
  
  coord_cols <- paste0("V", 1:ndim)
  
  # Calculate distances from previous year mean
  dist_data <- df_coords %>%
    mutate(ref_date = year - 1) %>%
    left_join(
      df_coords %>%
        group_by(year) %>%
        summarise(across(all_of(coord_cols), mean)),
      by = c("ref_date" = "year")
    ) %>%
    mutate(
      distance = sqrt(rowSums((
        dplyr::select(., paste0(coord_cols, ".x")) - 
          dplyr::select(., paste0(coord_cols, ".y")))^2))
    ) %>%
    dplyr::select(year, distance)
  
  # Calculate summary statistics
  overall_mean <- mean(dist_data$distance, na.rm=TRUE)
  overall_sd <- sd(dist_data$distance, na.rm=TRUE)
  
  return(list(
    dist_data = dist_data,
    summary = list(
      overall_mean = overall_mean,
      overall_sd = overall_sd
    )
  ))
}


#' Calculate Network Analysis Metrics
#'
#' @description
#' Analyzes the connectivity pattern in a distance matrix by converting it to
#' a network representation. Useful for assessing data completeness and structure.
#'
#' @param distance_matrix Square symmetric matrix of distances
#' @return A list containing the network analysis results:
#'   \item{adjacency}{A logical `matrix` where `TRUE` indicates a measured distance 
#'    between two points, representing the network's adjacency matrix.}
#'   \item{connectivity}{A `data.frame` with node-level metrics, including the `completeness` 
#'    (degree) for each point.}
#'   \item{summary}{A list of overall network statistics, including `n_points`, `n_measurements`, 
#'    and total `completeness`.}
#' @examples
#' # Create a sample distance matrix
#' dist_mat <- matrix(runif(25), 5, 5)
#' # Add row and column names
#' rownames(dist_mat) <- colnames(dist_mat) <- paste0("Point", 1:5)
#' dist_mat[lower.tri(dist_mat)] <- t(dist_mat)[lower.tri(dist_mat)]
#' diag(dist_mat) <- 0
#' dist_mat[1, 3] <- NA; dist_mat[3, 1] <- NA
#' 
#' # Analyze the network structure
#' metrics <- analyze_network_structure(dist_mat)
#' print(metrics$summary$completeness)
#'
#' @importFrom igraph graph_from_adjacency_matrix
#' @export
analyze_network_structure <- function(distance_matrix) {
  if (!is.matrix(distance_matrix)) {
    stop("Input must be a matrix")
  }
  
  if (nrow(distance_matrix) != ncol(distance_matrix)) {
    stop("Matrix must be square")
  }
  
  # Validate matrix properties
  if (!isSymmetric(unname(distance_matrix))) {
    stop("Matrix must be symmetric")
  }
  
  if (nrow(distance_matrix) < 2) {
    stop("Matrix must have at least 2 rows/columns")
  }

  # Validate connectivity
  if (all(is.na(distance_matrix))) {
    stop("Matrix contains no measurements")
  }
  
  # Create adjacency matrix
  adjacency <- !is.na(distance_matrix)
  
  # Create graph
  graph <- graph_from_adjacency_matrix(
    adjacency,
    mode = "undirected",
    weighted = TRUE
  )
  
  # Calculate node-level metrics
  connectivity <- data.frame(
    node = rownames(distance_matrix),
    completeness = rowSums(!is.na(distance_matrix)) / ncol(distance_matrix)
  )
  
  # Calculate overall metrics
  summary <- list(
    n_points = nrow(distance_matrix),
    n_measurements = sum(!is.na(distance_matrix)) / 2,
    completeness = sum(!is.na(distance_matrix)) / length(distance_matrix)
  )
  
  return(list(
    adjacency = adjacency,
    connectivity = connectivity,
    summary = summary
  ))
}


#' Generate Distance Matrix Heatmap Data
#'
#' @description
#' Prepares distance matrix data for heatmap visualization by handling missing values
#' and calculating relevant statistics.
#'
#' @param distance_matrix Square symmetric matrix of distances
#' @param cluster_rows Logical; whether to cluster rows
#' @param cluster_cols Logical; whether to cluster columns
#' @return A list of data prepared for generating a heatmap of the distance matrix:
#'   \item{matrix_data}{The distance `matrix`, potentially reordered by clustering.}
#'   \item{row_order}{An integer vector of the row indices after clustering. If `cluster_rows` 
#'    is `FALSE`, this is the original order.}
#'   \item{col_order}{An integer vector of the column indices after clustering. If `cluster_cols` 
#'    is `FALSE`, this is the original order.}
#'   \item{stats}{A list of summary statistics for the distance matrix, including `mean`, `sd`, 
#'    `min`, `max`, and `completeness`.}
#' @examples
#' # Create a sample distance matrix
#' dist_mat <- matrix(runif(25), 5, 5)
#' 
#' # Prepare data for a heatmap
#' heatmap_data <- prepare_heatmap_data(dist_mat)
#' print(heatmap_data$stats$completeness)
#'
#' @importFrom stats sd hclust dist
#' @export
prepare_heatmap_data <- function(distance_matrix, 
                                cluster_rows = FALSE,
                                cluster_cols = FALSE) {
  if (!is.matrix(distance_matrix)) {
    stop("Input must be a matrix")
  }
  
  # Validate sizes
  n_rows <- nrow(distance_matrix)
  n_cols <- ncol(distance_matrix)
  
  if (n_rows != n_cols) {
    stop("Matrix must be square")
  }
  
  if (n_rows < 2) {
    stop("Matrix must have at least 2 rows/columns")
  }
  
  # Validate data types
  if (!is.logical(cluster_rows)) {
    stop("cluster_rows must be logical")
  }
  
  if (!is.logical(cluster_cols)) {
    stop("cluster_cols must be logical")
  }
  
  # Check if clustering is possible
  if ((cluster_rows || cluster_cols) && 
      !requireNamespace("stats", quietly = TRUE)) {
    stop("stats package required for clustering")
  }
  
  # Calculate basic statistics
  stats <- list(
    mean = mean(distance_matrix, na.rm = TRUE),
    sd = sd(distance_matrix, na.rm = TRUE),
    min = min(distance_matrix, na.rm = TRUE),
    max = max(distance_matrix, na.rm = TRUE),
    completeness = sum(!is.na(distance_matrix)) / length(distance_matrix)
  )
  
  # Handle clustering if requested
  row_order <- seq_len(nrow(distance_matrix))
  col_order <- seq_len(ncol(distance_matrix))
  
  if (cluster_rows) {
    dist_mat <- dist(distance_matrix)
    row_order <- hclust(dist_mat)$order
  }
  
  if (cluster_cols) {
    dist_mat <- dist(t(distance_matrix))
    col_order <- hclust(dist_mat)$order
  }
  
  # Reorder matrix if needed
  matrix_data <- distance_matrix[row_order, col_order]
  
  return(list(
    matrix_data = matrix_data,
    row_order = row_order,
    col_order = col_order,
    stats = stats
  ))
}