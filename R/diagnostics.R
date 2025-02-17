# Copyright (c) 2024 Omid Arhami omid.arhami@uga.edu
# License: free of charge access granted to any academic researcher to use this software for non-commercial, academic research purposes **only**.  Nobody may modify, distribute, sublicense, or publicly share the Software or any derivative works, until the paper is published by the original authors.  The Software is provided "as is" without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose and noninfringement.  In no event shall the authors or copyright holders be liable for any claim, damages or other liability, whether in an action of contract, tort or otherwise, arising from, out of or in connection with the Software or the use or other dealings in the Software.

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
#' @importFrom stats var sd cor quantile
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
#' @return List containing:
#'   \item{converged}{Logical indicating if convergence achieved}
#'   \item{mean_converged}{Logical for mean convergence}
#'   \item{cov_converged}{Logical for covariance convergence}
#'   \item{final_mean}{Vector of final mean values}
#'   \item{final_cov}{Final covariance matrix}
#'   \item{mean_history}{Matrix of mean values over iterations}
#'   \item{cov_changes}{Vector of covariance changes}
#' @examples
#' \dontrun{
#' data <- read.csv("chain_data.csv")
#' conv_results <- check_gaussian_convergence(data)
#' print(conv_results)  # Shows summary
#' plot(conv_results)   # Creates convergence plots
#' }
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
#' @return List containing:
#'   \item{rhat}{R-hat statistic for each parameter}
#'   \item{ess}{Effective sample size for each parameter}
#' @examples
#' \dontrun{
#' chain_files <- c("chain1.csv", "chain2.csv", "chain3.csv")
#' diag <- calculate_diagnostics(chain_files, mutual_size = 1000)
#' print(diag)  # Shows R-hat and ESS
#' plot(diag)   # Creates trace and density plots
#' print(diag$rhat) # Should be close to 1
#' print(diag$ess)  # Should be large enough (>400) for reliable inference
#' }
#' @export
calculate_diagnostics <- function(chain_files, mutual_size=2000) {
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

  # Check file format
  lapply(chain_files, function(f) {
    tryCatch({
      df <- read.csv(f)
      required_cols <- c("NLL", "Holdout_MAE")
      if (!all(required_cols %in% names(df))) {
        stop("File ", f, " missing required columns: ",
             paste(setdiff(required_cols, names(df)), collapse = ", "))
      }
    }, error = function(e) {
      stop("Error reading file ", f, ": ", e$message)
    })
  })

  # Read chains
  chains <- lapply(chain_files, read.csv)
  
  # Get parameter names, excluding NLL/Holdout columns
  par_names <- c("log_N", "log_k0", "log_cooling_rate", "log_c_repulsion")
  
  # Process each chain
  for(i in 1:length(chains)) {
    chains[[i]] <- chains[[i]] %>% 
      filter(!is.na(NLL) & !is.na(Holdout_MAE) & 
               is.finite(NLL) & is.finite(Holdout_MAE))
    chains[[i]] <- na.omit(chains[[i]])
    chains[[i]] <- chains[[i]][
      (nrow(chains[[i]])-mutual_size-1):nrow(chains[[i]]),
      par_names]
  }
  
  # Check dimensions
  n_params <- unique(sapply(chains, ncol))
  if (length(n_params) != 1) {
    stop("All chains must have the same number of parameters")
  }
  
  # Convert to mcmc.list
  mcmc_list <- mcmc.list(lapply(chains, function(chain) {
    mcmc(as.matrix(chain))
  }))
  
  # Calculate diagnostics
  gelman_result <- coda::gelman.diag(mcmc_list)
  rhat <- gelman_result$psrf[,1]
  
  ess <- coda::effectiveSize(mcmc_list)
  
   structure(list(
    rhat = rhat,
    ess = ess,
    chains = chains,
    param_names = par_names,
    mutual_size = mutual_size
  ), class = "topolow_amcs_diagnostics")
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
#' @return List containing either:
#'   If reference_row provided:
#'   \item{summary_data}{Data frame with columns:
#'     \itemize{
#'       \item season_num: Numeric season identifier based on Influenza A.
#'       \item cluster: Cluster assignment  
#'       \item color: Point color
#'       \item avg_euclidean_dist: Mean distance to reference
#'       \item count: Points per cluster
#'       \item total_count: Total points per season
#'       \item fraction: Proportion of points in cluster
#'     }
#'   }
#'   If reference_row = FALSE:
#'   \item{dist_data}{Data frame with columns:
#'     \itemize{
#'       \item year_diff: Years between points
#'       \item euclidean_dist: Distance between points
#'       \item ref_year: Reference year
#'     }
#'   }
#' @examples
#' \dontrun{
#' # Calculate distances from reference point
#' ref_distances <- calculate_cumulative_distances(coords, ndim=2, reference_row=1)
#'
#' # Calculate all pairwise distances
#' all_distances <- calculate_cumulative_distances(coords, ndim=2, reference_row=FALSE)
#' }
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
                count = n(), .groups = 'drop') %>%
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
#' @return List containing:
#'   \item{dist_data}{Data frame with columns:
#'     \itemize{
#'       \item year: Collection year
#'       \item distance: Distance from previous year mean
#'     }
#'   }
#'   \item{summary}{List with:
#'     \itemize{
#'       \item overall_mean: Mean distance across all years
#'       \item overall_sd: Standard deviation of distances
#'     }
#'   }
#' @examples
#' \dontrun{
#' annual_stats <- calculate_annual_distances(coords, ndim=2)
#' print(annual_stats$summary$overall_mean)
#' }
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
#' @return List containing:
#'   \item{adjacency}{Logical matrix indicating presence of measurements}
#'   \item{connectivity}{Data frame with connectivity metrics per point}
#'   \item{summary}{List of overall network statistics}
#' @examples
#' \dontrun{
#' metrics <- analyze_network_structure(dist_mat)
#' print(metrics$summary$completeness)
#' }
#' @importFrom igraph graph_from_adjacency_matrix degree
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
#' @return List containing:
#'   \item{matrix_data}{Processed matrix for visualization}
#'   \item{row_order}{Optional row ordering from clustering}
#'   \item{col_order}{Optional column ordering from clustering}
#'   \item{stats}{List of matrix statistics}
#' @examples
#' \dontrun{
#' heatmap_data <- prepare_heatmap_data(dist_mat)
#' print(heatmap_data$stats$completeness)
#' }
#' @importFrom stats hclust dist
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