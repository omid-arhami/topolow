# Copyright (c) 2024 Omid Arhami omid.arhami@uga.edu
# R/diagnostics.R

#' Model Diagnostics and Convergence Testing


#' Create Diagnostic Plots for Multiple Sampling Chains
#'
#' @description
#' Creates trace and density plots for multiple sampling or optimization chains to help
#' assess convergence and mixing. It displays parameter trajectories and their
#' distributions across all chains.
#'
#' @param chain_files A character vector of paths to CSV files, where each file contains data for one chain.
#' @param mutual_size Integer. The number of samples to use from the end of each chain for plotting.
#' @param output_file Character. The path for saving the plot. Required if `save_plot` is TRUE.
#' @param output_dir Character. The directory for saving output files. Required if `save_plot` is TRUE.
#' @param save_plot Logical. If TRUE, saves the plot to a file. Default: FALSE.
#' @param width,height,dpi Numeric. The dimensions and resolution for the saved plot.
#' @return A `ggplot` object of the combined plots.
#' @examples
#' # This example uses sample data files that would be included with the package.
#' chain_files <- c(
#'   system.file("extdata", "diag_chain1.csv", package = "topolow"),
#'   system.file("extdata", "diag_chain2.csv", package = "topolow"),
#'   system.file("extdata", "diag_chain3.csv", package = "topolow")
#' )
#'
#' # Only run the example if the files are found
#' if (all(nzchar(chain_files))) {
#'   # Create diagnostic plot without saving to a file
#'   plot_mcmc_diagnostics(chain_files, mutual_size = 50, save_plot = FALSE)
#' }
#'
#' @export
plot_mcmc_diagnostics <- function(chain_files,
                                    mutual_size = 2000,
                                    output_file = "diagnostic_plots.png",
                                    output_dir,
                                    save_plot = FALSE,
                                    width = 3000, height = 3000, dpi = 300) {
                                      # Check if gridextra is available
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("gridExtra package is required for plotting. Please install with install.packages('gridExtra').")
  }
  if (save_plot && (missing(output_dir) || missing(output_file))) {
    stop("`output_dir` and `output_file` must be provided when `save_plot` is TRUE.", call. = FALSE)
  }

  # Read and process chain data
  par_names <- c("log_N", "log_k0", "log_cooling_rate", "log_c_repulsion")
  chains <- lapply(chain_files, function(file) {
    df <- utils::read.csv(file)
    # Take the last `mutual_size` samples
    tail(df[, par_names], mutual_size)
  })

  n_params <- ncol(chains[[1]])
  n_chains <- length(chains)

  # Create a list to hold all the individual plots
  plot_list <- list()

  for (i in seq_len(n_params)) {
    # Combine data from all chains for the current parameter
    trace_data <- do.call(rbind, lapply(seq_len(n_chains), function(j) {
      data.frame(
        Chain = as.factor(j),
        Iteration = seq_len(nrow(chains[[j]])),
        Value = chains[[j]][, i]
      )
    }))

    # Create Trace Plot
    p_trace <- ggplot2::ggplot(trace_data,
                               ggplot2::aes(x = .data$Iteration, y = .data$Value, color = .data$Chain)) +
      ggplot2::geom_line(linewidth = 0.5) +
      ggplot2::labs(title = paste("Trace Plot:", par_names[i]), x = "Iteration", y = "Value") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "none",
        panel.grid.minor = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(color = "black", fill = NA)
      )
    plot_list <- c(plot_list, list(p_trace))

    # Create Density Plot
    p_density <- ggplot2::ggplot(trace_data, ggplot2::aes(x = .data$Value, color = .data$Chain)) +
      ggplot2::geom_density(alpha = 0.3) +
      ggplot2::labs(title = paste("Density Plot:", par_names[i]), x = "Value", y = "Density") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "none",
        panel.grid.minor = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(color = "black", fill = NA)
      )
    plot_list <- c(plot_list, list(p_density))
  }

  # Arrange all plots into a grid
  combined_plot <- gridExtra::grid.arrange(grobs = plot_list, ncol = 2)

  # Optionally save the combined plot
  if (save_plot) {
    full_output_path <- file.path(output_dir, output_file)
    ggsave_white_bg(full_output_path, combined_plot,
                    width = width / dpi, height = height / dpi,
                    dpi = dpi, limitsize = FALSE)
  }

  return(combined_plot)
}


#' Plot Embedding Quality
#'
#' @description
#' Creates diagnostic plots for assessing embedding quality.
#'
#' @param true_dissim Matrix of input dissimilarities.
#' @param est_dissim Matrix of estimated distances from embedding.
#'
#' @return A ggplot object with quality diagnostic plots.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_abline geom_density2d
#'   labs theme_minimal stat_smooth
#' @keywords internal
plot_embedding_quality <- function(true_dissim, est_dissim) {
  
  # Get observed (non-NA, non-diagonal) pairs
  mask <- !is.na(true_dissim) & row(true_dissim) != col(true_dissim)
  
  plot_data <- data.frame(
    true_dist = as.vector(true_dissim[mask]),
    est_dist = as.vector(est_dissim[mask])
  )
  
  # Remove any character values (thresholds)
  plot_data$true_dist <- suppressWarnings(as.numeric(as.character(plot_data$true_dist)))
  plot_data <- plot_data[!is.na(plot_data$true_dist), ]
  
  p1 <- ggplot(plot_data, aes(x = .data$true_dist, y = .data$est_dist)) +
    geom_point(alpha = 0.3, size = 1.5) +
    geom_abline(slope = 1, intercept = 0, color = "red", 
                linetype = "dashed", linewidth = 1) +
    stat_smooth(method = "loess", color = "blue", se = FALSE) +
    labs(
      title = "Embedding Quality: Predicted vs True Distances",
      subtitle = "Red line = perfect fit, Blue line = LOESS smooth",
      x = "Input Dissimilarity",
      y = "Estimated Euclidean Distance"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10)
    )
  
  # Residuals
  plot_data$residual <- plot_data$est_dist - plot_data$true_dist
  
  p2 <- ggplot(plot_data, aes(x = .data$true_dist, y = .data$residual)) +
    geom_point(alpha = 0.3, size = 1.5) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    stat_smooth(method = "loess", color = "blue", se = TRUE) +
    labs(
      title = "Residual Plot",
      x = "Input Dissimilarity",
      y = "Residual (Estimated - True)"
    ) +
    theme_minimal()
  
  p3 <- ggplot(plot_data, aes(x = .data$residual)) +
    geom_density(fill = "steelblue", alpha = 0.5) +
    labs(
      title = "Distribution of Residuals",
      x = "Residual",
      y = "Density"
    ) +
    theme_minimal()
  
  gridExtra::grid.arrange(p1, p2, p3, ncol = 1, heights = c(2, 1, 1))
}


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


#' Plot Running Minimum Error to Show Sampling Convergence
#'
#' @description
#' Creates a diagnostic plot showing how the best (minimum) NLL or MAE found improves
#' over iterations and eventually plateaus. This is a standard way to demonstrate that
#' sampling/optimization has converged in terms of the objective function.
#'
#' @param chain_files Character vector. Paths to CSV files containing sampling chains.
#' @param metric Character. Which metric to plot: "NLL", "MAE", or "both". Default: "both"
#' @param combine_chains Logical. If TRUE, combines all chains into one sequence.
#'   If FALSE, plots each chain separately. Default: TRUE
#' @param sort_combined Logical. If TRUE and `combine_chains` is TRUE, sorts combined
#'   samples from worst to best for each metric independently before computing the
#'   running minimum. This produces a smooth monotonic improvement curve that is not
#'   biased by chain ordering. If FALSE, chains are concatenated sequentially (original
#'   behavior). Default: TRUE
#' @param show_raw Logical. If TRUE, shows raw values as points behind the running minimum.
#'   Default: TRUE
#' @param window_size Integer. Window size for computing rolling mean (optional smoothing).
#'   Set to NULL for no smoothing. Default: NULL
#' @param output_file Character. Optional path to save the plot. Default: NULL
#' @param width,height Numeric. Dimensions for saved plot. Default: 10, 6
#' @param dpi Numeric. Resolution for saved plot. Default: 300
#'
#' @return A ggplot object showing the error plateau diagnostic.
#'
#' @details
#' The "running minimum" (or cumulative minimum) shows the best value found so far
#' at each iteration. When this line flattens out (plateaus), it indicates that
#' the sampling has found the optimal region and additional samples are not
#' improving the objective.
#'
#' When `combine_chains = TRUE` and `sort_combined = TRUE`, samples from all chains
#' are sorted from worst (highest) to best (lowest) independently for each metric
#' before computing the running minimum. This avoids the problem where the first
#' chain in the list dominates the visual, and produces a smooth improvement curve
#' showing how quality improves across the full pool of samples. Each metric facet
#' has its own sort order, so the NLL panel sorts by NLL and the MAE panel sorts
#' by MAE.
#'
#' This is distinct from:
#' - Trace plots (which show parameter values, not objective function)
#' - R-hat (which measures between-chain variance)
#' - ESS (which measures autocorrelation)
#'
#' @examples
#' \dontrun{
#' # Basic usage (sorted by default)
#' chain_files <- c("chain1.csv", "chain2.csv", "chain3.csv")
#' plot_performance_trace(chain_files)
#'
#' # Original sequential concatenation behavior
#' plot_performance_trace(chain_files, sort_combined = FALSE)
#'
#' # Plot only NLL with separate chains
#' plot_performance_trace(chain_files, metric = "NLL", combine_chains = FALSE)
#'
#' # Save the plot
#' plot_performance_trace(chain_files, output_file = "fig/convergence_plateau.png")
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_step labs theme_minimal
#'   scale_color_manual facet_wrap theme element_text ggsave
#' @importFrom dplyr mutate group_by ungroup bind_rows
#' @export
plot_performance_trace <- function(chain_files,
                                metric = "both",
                                combine_chains = TRUE,
                                sort_combined = TRUE,
                                show_raw = TRUE,
                                window_size = NULL,
                                output_file = NULL,
                                width = 10,
                                height = 6,
                                dpi = 300) {
  
  # Validate inputs
  if (!all(file.exists(chain_files))) {
    missing <- chain_files[!file.exists(chain_files)]
    stop("Missing chain files: ", paste(missing, collapse = ", "))
  }
  metric <- match.arg(metric, c("NLL", "MAE", "both"))
  
  # Read and combine chain data
  all_data <- lapply(seq_along(chain_files), function(i) {
    df <- utils::read.csv(chain_files[i])
    df$chain <- paste0("Chain ", i)
    df$sample_in_chain <- seq_len(nrow(df))
    return(df)
  })
  
  if (combine_chains) {
    combined <- do.call(rbind, all_data)
    combined$chain <- "Combined"
    
    if (sort_combined) {
      # Sort independently per metric for a smooth improvement curve
      # ----------------------------------------------------------
      # NLL: sort worst (highest) to best (lowest)
      nll_order <- order(-combined$NLL)
      combined$iteration_nll <- NA_integer_
      combined$iteration_nll[nll_order] <- seq_len(nrow(combined))
      
      # MAE: sort worst (highest) to best (lowest)
      mae_order <- order(-combined$Holdout_MAE)
      combined$iteration_mae <- NA_integer_
      combined$iteration_mae[mae_order] <- seq_len(nrow(combined))
      
      # Compute running minimums along the sorted orders
      combined$running_min_NLL <- NA_real_
      combined$running_min_NLL[nll_order] <- cummin(combined$NLL[nll_order])
      
      combined$running_min_MAE <- NA_real_
      combined$running_min_MAE[mae_order] <- cummin(combined$Holdout_MAE[mae_order])
      
      # Optional: rolling mean along sorted orders
      if (!is.null(window_size) && window_size > 1) {
        combined$rolling_mean_NLL <- NA_real_
        combined$rolling_mean_NLL[nll_order] <- zoo::rollmean(
          combined$NLL[nll_order], k = window_size, fill = NA, align = "right")
        
        combined$rolling_mean_MAE <- NA_real_
        combined$rolling_mean_MAE[mae_order] <- zoo::rollmean(
          combined$Holdout_MAE[mae_order], k = window_size, fill = NA, align = "right")
      }
      
      # For single-metric mode, set a default iteration column
      if (metric == "NLL") {
        combined$iteration <- combined$iteration_nll
      } else if (metric == "MAE") {
        combined$iteration <- combined$iteration_mae
      }
      # For "both", iteration_nll and iteration_mae are used directly below
      
      plot_data <- combined
      
    } else {
      # Original behavior: sequential concatenation
      combined$iteration <- seq_len(nrow(combined))
      plot_data <- combined
    }
    
  } else {
    # Keep chains separate but add global iteration counter
    plot_data <- do.call(rbind, all_data)
    plot_data$iteration <- plot_data$sample_in_chain
  }
  
  # Filter out invalid rows
  plot_data <- plot_data[!is.na(plot_data$NLL) & !is.na(plot_data$Holdout_MAE) &
                           is.finite(plot_data$NLL) & is.finite(plot_data$Holdout_MAE), ]
  
  # Calculate running minimum for each chain (non-sorted path)
  if (!(combine_chains && sort_combined)) {
    plot_data <- do.call(rbind, lapply(split(plot_data, plot_data$chain), function(chain_df) {
      chain_df <- chain_df[order(chain_df$iteration), ]
      chain_df$running_min_NLL <- cummin(chain_df$NLL)
      chain_df$running_min_MAE <- cummin(chain_df$Holdout_MAE)
      
      # Optional: rolling mean for smoothing
      if (!is.null(window_size) && window_size > 1) {
        chain_df$rolling_mean_NLL <- zoo::rollmean(chain_df$NLL, k = window_size, 
                                                     fill = NA, align = "right")
        chain_df$rolling_mean_MAE <- zoo::rollmean(chain_df$Holdout_MAE, k = window_size, 
                                                     fill = NA, align = "right")
      }
      return(chain_df)
    }))
  }
  
  # Prepare data for plotting
  if (metric == "both") {
    
    if (combine_chains && sort_combined) {
      # Each metric gets its own iteration axis (sorted independently)
      plot_long <- rbind(
        data.frame(
          iteration = plot_data$iteration_nll,
          chain = plot_data$chain,
          value = plot_data$NLL,
          running_min = plot_data$running_min_NLL,
          metric = "Log-Likelihood"
        ),
        data.frame(
          iteration = plot_data$iteration_mae,
          chain = plot_data$chain,
          value = plot_data$Holdout_MAE,
          running_min = plot_data$running_min_MAE,
          metric = "Mean Absolute Error (MAE)"
        )
      )
    } else {
      # Shared iteration axis (original behavior)
      plot_long <- rbind(
        data.frame(
          iteration = plot_data$iteration,
          chain = plot_data$chain,
          value = plot_data$NLL,
          running_min = plot_data$running_min_NLL,
          metric = "Log-Likelihood"
        ),
        data.frame(
          iteration = plot_data$iteration,
          chain = plot_data$chain,
          value = plot_data$Holdout_MAE,
          running_min = plot_data$running_min_MAE,
          metric = "Mean Absolute Error (MAE)"
        )
      )
    }
    
    p <- ggplot2::ggplot(plot_long, ggplot2::aes(x = iteration))
    
    if (show_raw) {
      p <- p + ggplot2::geom_point(ggplot2::aes(y = value, color = chain), 
                                    alpha = 0.3, size = 1)
    }
    
    subtitle_text <- if (combine_chains && sort_combined) {
      "Samples sorted worst-to-best per metric; plateau indicates optimal region"
    } else {
      "Plateau indicates optimal region has been found"
    }
    
    p <- p + 
      ggplot2::geom_step(ggplot2::aes(y = running_min, color = chain), 
                          linewidth = 1.2, direction = "hv") +
      ggplot2::facet_wrap(~metric, scales = "free", ncol = 1) +
      ggplot2::labs(
        title = "Sampling Performance: Running Minimum Error",
        subtitle = subtitle_text,
        x = "Sample (sorted by metric)",
        y = "Value",
        color = "Chain"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 14, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 10, color = "gray40"),
        strip.text = ggplot2::element_text(size = 11, face = "bold"),
        panel.grid.minor = ggplot2::element_blank()
      )
    
  } else {
    # Single metric
    value_col <- if (metric == "NLL") "NLL" else "Holdout_MAE"
    running_min_col <- if (metric == "NLL") "running_min_NLL" else "running_min_MAE"
    y_label <- if (metric == "NLL") "Negative Log-Likelihood" else "Mean Absolute Error"
    
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = iteration))
    
    if (show_raw) {
      p <- p + ggplot2::geom_point(ggplot2::aes_string(y = value_col, color = "chain"), 
                                    alpha = 0.3, size = 1)
    }
    
    subtitle_text <- if (combine_chains && sort_combined) {
      "Samples sorted worst-to-best; plateau indicates optimal region"
    } else {
      "Plateau indicates optimal region has been found"
    }
    
    p <- p + 
      ggplot2::geom_step(ggplot2::aes_string(y = running_min_col, color = "chain"), 
                          linewidth = 1.2, direction = "hv") +
      ggplot2::labs(
        title = paste("Sampling Performance: Running Minimum", metric),
        subtitle = subtitle_text,
        x = "Sample (sorted by metric)",
        y = y_label,
        color = "Chain"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 14, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 10, color = "gray40"),
        panel.grid.minor = ggplot2::element_blank()
      )
  }
  
  # Remove legend if only one chain
  if (length(unique(plot_data$chain)) == 1) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  
  # Save if requested
  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, p, width = width, height = height, dpi = dpi, bg = "white")
    message("Plot saved to: ", output_file)
  }
  
  return(p)
}


#' Plot Log-Likelihood Improvement Over Iterations
#'
#' @description
#' A visualization showing the improvement in log-likelihood 
#' relative to the first sample, making it easy to see when gains diminish.
#'
#' @param chain_files Character vector. Paths to CSV files containing sampling chains.
#' @param combine_chains Logical. If TRUE, combines all chains. Default: TRUE
#'
#' @return A ggplot object.
#' @export
plot_ll_improvement <- function(chain_files, combine_chains = TRUE) {
  
  # Read and combine chain data
  all_data <- lapply(seq_along(chain_files), function(i) {
    df <- utils::read.csv(chain_files[i])
    df$chain <- paste0("Chain ", i)
    df$sample_in_chain <- seq_len(nrow(df))
    return(df)
  })
  
  if (combine_chains) {
    plot_data <- do.call(rbind, all_data)
    plot_data$iteration <- seq_len(nrow(plot_data))
    plot_data$chain <- "Combined"
  } else {
    plot_data <- do.call(rbind, all_data)
    plot_data$iteration <- plot_data$sample_in_chain
  }
  
  # Filter out invalid rows
  plot_data <- plot_data[!is.na(plot_data$NLL) & is.finite(plot_data$NLL), ]
  
  # Calculate LL improvement (since we minimize NLL, improvement = reduction)
  plot_data <- do.call(rbind, lapply(split(plot_data, plot_data$chain), function(chain_df) {
    chain_df <- chain_df[order(chain_df$iteration), ]
    initial_best <- chain_df$NLL[1]
    chain_df$running_min_NLL <- cummin(chain_df$NLL)
    chain_df$improvement <- initial_best - chain_df$running_min_NLL
    return(chain_df)
  }))
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = iteration, y = improvement, color = chain)) +
    ggplot2::geom_step(linewidth = 1.2) +
    ggplot2::labs(
      title = "Cumulative NLL Improvement Over Sampling",
      subtitle = "Flattening curve indicates convergence to optimal region",
      x = "Sample Iteration",
      y = "NLL Improvement (relative to initial best)",
      color = "Chain"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10, color = "gray40")
    )
  
  if (length(unique(plot_data$chain)) == 1) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  
  return(p)
}