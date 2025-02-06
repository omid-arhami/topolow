# Copyright (c) 2024 Omid Arhami omid.arhami@uga.edu
# Licensed under MIT License
# R/utils.R

#' Utility functions for the topolow package
#' 
#' @description
#' This file contains utility functions used throughout the topolow package.
#' Functions include data manipulation, format conversion, and optional SLURM support.
#'
#' @keywords internal
"_PACKAGE"

#' (Not used anywhere) Scale matrix columns to 0-10 range
#'
#' @param mat Numeric matrix to scale
#' @return Matrix with columns scaled to 0-10 range
#' @keywords internal
scale_columns_0_to_10 <- function(mat) {
  if (!is.matrix(mat)) {
    stop("Input must be a matrix")
  }
  
  scale_vector <- function(x) {
    non_na_x <- x[!is.na(x)]
    if (length(non_na_x) == 0) return(x)
    
    min_val <- min(non_na_x)
    max_val <- max(non_na_x)
    
    if (min_val == max_val) {
      return(ifelse(is.na(x), NA, 0))
    }
    
    scaled_x <- (x - min_val) / (max_val - min_val) * 10
    return(scaled_x)
  }
  
  scaled_mat <- apply(mat, 2, scale_vector)
  return(scaled_mat)
}


#' Convert 2-digit to 4-digit year
#'
#' @param dataframe Data frame containing year column
#' @param column_name Name of year column
#' @return Data frame with converted year column
#' @keywords internal
yy_to_yyyy <- function(dataframe, column_name) {
  if (!is.data.frame(dataframe)) {
    stop("First argument must be a data frame")
  }
  if (!column_name %in% names(dataframe)) {
    stop(paste("Column", column_name, "not found in data frame"))
  }
 
  # Validate year values
  years <- dataframe[[column_name]]
  if (!is.numeric(years)) {
    stop("Year column must be numeric")
  }
  
  if (any(years < 0, na.rm = TRUE)) {
    stop("Year values cannot be negative")
  }
  
  dataframe[[column_name]] <- sapply(dataframe[[column_name]], function(x) {
    year <- as.integer(x)
    if (year < 25) {
      return(year + 2000)
    } else if (year < 100) {
      return(year + 1900) 
    }
    return(year)
  })
  
  return(dataframe)
}


#' Extract year from sequence name
#'
#' @param seq_name Character string containing sequence name
#' @return Extracted year as integer
#' @keywords internal
extract_year <- function(seq_name) {
  if (!is.character(seq_name)) {
    stop("seq_name must be a character string")
  }
  
  # Check for year pattern
  if (!grepl("\\d{4}", seq_name)) {
    stop("No 4-digit year found in sequence name")
  }
  
  year <- as.integer(sub(".*\\/(\\d{4}).*", "\\1", seq_name))
  
  # Validate extracted year
  if (is.na(year) || year < 1900 || year > 2100) {
    stop("Invalid year value extracted: ", year)
  }
  
  return(year)
}


#' Sort symmetric matrix by year
#'
#' @param matrix symmetric matrix
#' @return Sorted symmetric matrix
#' @keywords internal
sort_matrix_by_year <- function(matrix) {
  # Verify matrix is symmetric
  if (!isSymmetric(matrix)) {
    stop("Input matrix must be symmetric")
  }
  
  # Get row names and extract years
  names <- rownames(matrix)
  years <- sapply(names, extract_year)
  
  # Create ordering based on years
  order_idx <- order(years)
  
  # Create new sorted matrix
  sorted_matrix <- matrix[order_idx, order_idx]
  
  # Verify that the sorting preserved symmetry
  if (!isSymmetric(sorted_matrix)) {
    warning("Resulting matrix is not symmetric - something went wrong")
  }
  
  # Return sorted matrix
  return(sorted_matrix)
}


#' Generate unique string identifiers with year suffix
#'
#' @param n Number of strings to generate
#' @param length Length of random part of string (default: 8)
#' @param lower_bound Lower bound for year suffix (default: 1)
#' @param upper_bound Upper bound for year suffix (default: 20)
#' @return Character vector of unique strings with year suffixes
#' @export
generate_unique_string <- function(n, length = 8, lower_bound = 1, upper_bound = 20) {
  if (n <= 0) stop("n must be positive")
  if (length <= 0) stop("length must be positive")
  if (lower_bound > upper_bound) {
    stop("lower_bound must be less than or equal to upper_bound")
  }
  
  # Validate numeric parameters
  if (!all(sapply(list(n, length, lower_bound, upper_bound), is.numeric))) {
    stop("All numeric parameters must be numeric")
  }
  
  if (!all(sapply(list(n, length), function(x) x == round(x)))) {
    stop("n and length must be integers")
  }

  # Generate base strings
  base_strings <- sapply(1:n, function(x) {
    paste0(sample(c(LETTERS, 0:9), length, replace = TRUE), collapse = "")
  })
  
  # Calculate strings per year
  num_per_group <- n / (upper_bound - lower_bound + 1)
  
  # Generate year suffixes
  suffixes <- rep(lower_bound:upper_bound, each = ceiling(num_per_group))[1:n]
  
  # Combine strings with suffixes
  unique_strings <- paste0(base_strings, "/", suffixes)
  
  return(unique_strings)
}


#' Extract last part of a path or name after final slash 
#'
#' @param name Character string containing path or name
#' @return Last part of string after final slash
#' @keywords internal
extract_last_part <- function(name) {
  if (!is.character(name)) {
    stop("name must be a character string")
  }
  parts <- strsplit(name, "/")[[1]]
  return(tail(parts, n=1))
}


#' Generate Complex High-Dimensional Data for Testing
#'
#' @description 
#' Generates synthetic high-dimensional data with clusters and trends for testing 
#' dimensionality reduction methods. Creates data with specified properties:
#' - Multiple clusters along a trend line
#' - Variable density regions 
#' - Controllable noise levels
#' - Optional visualization
#'
#' The function generates cluster centers along a trend line, adds points around those
#' centers with specified spread, and incorporates random noise to create high and 
#' low density areas. The data is useful for testing dimensionality reduction and
#' visualization methods.
#'
#' @param n_points Integer number of points to generate
#' @param n_dim Integer number of dimensions 
#' @param n_clusters Integer number of clusters
#' @param cluster_spread Numeric controlling cluster variance
#' @param fig_name Character path to save visualization (optional)
#' @return Data frame with generated coordinates in n_dim dimensions. Column names 
#'         are "Dim1" through "DimN" where N is n_dim.
#' @examples
#' \dontrun{
#' # Generate basic dataset
#' data <- generate_complex_data(n_points = 500, n_dim = 10, 
#'                              n_clusters = 4, cluster_spread = 1)
#'                              
#' # Generate and visualize dataset
#' data <- generate_complex_data(n_points = 500, n_dim = 10,
#'                              n_clusters = 4, cluster_spread = 1,
#'                              fig_name = "cluster_viz.png")
#' }
#' @importFrom MASS mvrnorm
#' @importFrom stats runif rnorm prcomp
#' @importFrom ggplot2 ggplot aes geom_point theme_minimal labs coord_fixed
#' @export
generate_complex_data <- function(n_points = 500, n_dim = 10, n_clusters = 4, 
                                cluster_spread = 1, fig_name = NA) {
  # Input validation
  if (!is.numeric(n_points) || n_points <= 0 || n_points %% 1 != 0) {
    stop("n_points must be a positive integer")
  }
  if (!is.numeric(n_dim) || n_dim <= 0 || n_dim %% 1 != 0) {
    stop("n_dim must be a positive integer")
  }
  if (!is.numeric(n_clusters) || n_clusters <= 0 || n_clusters %% 1 != 0) {
    stop("n_clusters must be a positive integer")
  }
  if (!is.numeric(cluster_spread) || cluster_spread <= 0) {
    stop("cluster_spread must be a positive number")
  }
  if (!is.na(fig_name) && !is.character(fig_name)) {
    stop("fig_name must be NA or a character string ending in .png")
  }
  
  if (n_clusters > n_points) {
    stop("Number of clusters cannot exceed number of points")
  }
  
  # Generate a trend vector
  trend <- seq(-10, 10, length.out = n_points)
  
  # Generate cluster centers along the trend
  cluster_centers <- matrix(0, nrow = n_clusters, ncol = n_dim)
  for (i in 1:n_clusters) {
    cluster_centers[i, ] <- 0.5 * runif(n_dim, min = -10, max = 10) + 
      trend[floor(i * n_points / n_clusters)]
  }
  
  # Generate points around each cluster center
  points_per_cluster <- floor(n_points / n_clusters)
  data <- do.call(rbind, lapply(1:n_clusters, function(i) {
    mvrnorm(n = points_per_cluster, 
            mu = cluster_centers[i, ], 
            Sigma = diag(cluster_spread, n_dim))
  }))
  
  # Add random noise to create high and low density areas
  noise <- matrix(rnorm(nrow(data) * n_dim, mean = 0, sd = 3.3), 
                 ncol = n_dim)
  data <- data + noise
  
  # Add trend to multiple dimensions
  for (j in 1:n_dim) {
    data[, j] <- data[, j] + trend[1:nrow(data)]
  }
  
  # Convert to data frame with proper column names
  data <- as.data.frame(data)
  colnames(data) <- paste0("Dim", 1:n_dim)
  
  # Create visualization if requested
  if (!is.na(fig_name)) {
    # Perform PCA for visualization
    pca_result <- prcomp(data, scale. = FALSE)
    pca_data <- as.data.frame(pca_result$x[, 1:2])
    colnames(pca_data) <- c("PC1", "PC2")
    
    # Create plot
    p <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
      geom_point(alpha = 0.5) +
      theme_minimal() +
      labs(title = "PCA of Simulated High-Dimensional Data",
           x = "PC 1", y = "PC 2") +
      coord_fixed()
    
    # Save plot
    ggsave(filename = fig_name, plot = p, dpi = 300)
  }
  
  return(data)
}



#' Increase Missing Values in a Matrix
#'
#' @description
#' Strategically introduces NA values into a distance matrix while maintaining symmetry.
#' New NA values are added preferentially farther from the diagonal to simulate real-world 
#' measurement patterns where distant pairs are more likely to be unmeasured.
#'
#' @details
#' The function:
#' 1. Calculates needed additional NAs to reach target percentage
#' 2. Creates probability matrix favoring off-diagonal elements
#' 3. Randomly selects positions weighted by distance from diagonal
#' 4. Maintains matrix symmetry by mirroring NAs
#'
#' @param mat Matrix to modify
#' @param target_na_percentage Numeric between 0 and 1 specifying desired proportion of NAs
#' @return Matrix with increased NA values, maintaining symmetry
#' @examples
#' \dontrun{
#' # Create sample distance matrix
#' dist_mat <- matrix(runif(100), 10, 10)
#' dist_mat[lower.tri(dist_mat)] <- t(dist_mat)[lower.tri(dist_mat)]
#' diag(dist_mat) <- 0
#' 
#' # Increase NAs to 70%
#' sparse_mat <- increase_na_percentage(dist_mat, 0.7)
#' }
#' @export
increase_na_percentage <- function(mat, target_na_percentage) {
  # Input validation
  if (!is.matrix(mat)) {
    stop("Input must be a matrix")
  }
  if (!is.numeric(target_na_percentage) || 
      target_na_percentage <= 0 || 
      target_na_percentage >= 1) {
    stop("target_na_percentage must be between 0 and 1")
  }
  
  # Calculate current and target NA counts
  total_elements <- length(mat)
  current_na_count <- sum(is.na(mat))
  target_na_count <- ceiling(target_na_percentage * total_elements)
  additional_na_count <- target_na_count - current_na_count
  
  if (additional_na_count <= 0) {
    warning("The matrix already has more NA values than the target percentage.")
    return(mat)
  }
  
  # Get indices of non-NA elements
  non_na_indices <- which(!is.na(mat))
  
  # Calculate distance from diagonal for each element
  n <- nrow(mat)
  distances <- outer(1:n, 1:n, FUN = function(i, j) abs(i - j))
  
  # Create probability matrix favoring off-diagonal elements
  max_distance <- max(distances)
  # Add padding to avoid strictly diagonal matrix
  distances <- distances + 0.15 * max_distance
  
  # Normalize probabilities
  prob_matrix <- distances / max_distance
  prob_vector <- as.vector(prob_matrix)
  prob_vector <- prob_vector / sum(prob_vector)
  
  # Randomly sample indices based on probability vector
  random_indices <- sample(non_na_indices, 
                         size = additional_na_count,
                         prob = prob_vector[non_na_indices], 
                         replace = FALSE)
  
  # Replace selected elements with NA while maintaining symmetry
  for (index in random_indices) {
    row <- (index - 1) %/% nrow(mat) + 1
    col <- (index - 1) %% ncol(mat) + 1
    if (row != col && !is.na(mat[row, col])) {
      mat[row, col] <- NA
      mat[col, row] <- NA
    }
  }
  
  return(mat)
}

#' Add Noise and Bias to Matrix Data
#'
#' @description
#' Creates noisy versions of a distance matrix by adding random noise and/or systematic bias.
#' Useful for testing robustness of algorithms to measurement errors and systematic biases.
#'
#' @details
#' The function generates three variants of the input matrix:
#' 1. n1: Matrix with random Gaussian noise
#' 2. n2: Different realization of random noise
#' 3. nb: Matrix with both random noise and systematic negative bias
#'
#' The noise level is scaled relative to the data mean to maintain realistic error magnitudes.
#'
#' @param matrix_data Numeric matrix to add noise to
#' @return List containing three matrices:
#'   \item{n1}{Matrix with first noise realization}
#'   \item{n2}{Matrix with second noise realization}
#'   \item{nb}{Matrix with noise and negative bias}
#' @examples
#' \dontrun{
#' # Create sample distance matrix
#' dist_mat <- matrix(runif(100), 10, 10)
#' dist_mat[lower.tri(dist_mat)] <- t(dist_mat)[lower.tri(dist_mat)]
#' diag(dist_mat) <- 0
#'
#' # Generate noisy versions
#' noisy_variants <- add_noise_bias(dist_mat)
#' }
#' @importFrom stats rnorm
#' @export
add_noise_bias <- function(matrix_data) {
  # Input validation
  if (!is.matrix(matrix_data) || !is.numeric(matrix_data)) {
    stop("Input must be a numeric matrix")
  }
  
  # Calculate noise scale based on data mean
  data_mean <- mean(matrix_data, na.rm = TRUE)
  measurement_er <- 0.05 * data_mean
  
  # Generate first noisy variant
  noise1 <- matrix(rnorm(length(matrix_data), 
                        mean = 0, 
                        sd = measurement_er), 
                  nrow = nrow(matrix_data), 
                  ncol = ncol(matrix_data))
  noisy_matrix1 <- matrix_data + noise1
  
  # Enforce non-negative values and zero diagonal
  noisy_matrix1[noisy_matrix1 < 0] <- 0
  diag(noisy_matrix1) <- 0
  
  # Generate second noisy variant
  noise2 <- matrix(rnorm(length(matrix_data), 
                        mean = 0, 
                        sd = measurement_er),
                  nrow = nrow(matrix_data), 
                  ncol = ncol(matrix_data))
  noisy_matrix2 <- matrix_data + noise2
  noisy_matrix2[noisy_matrix2 < 0] <- 0
  diag(noisy_matrix2) <- 0
  
  # Generate biased variant
  bias <- matrix(measurement_er, 
                nrow = nrow(matrix_data), 
                ncol = ncol(matrix_data))
  noisy_biased_matrix <- matrix_data + noise2 - bias
  noisy_biased_matrix[noisy_biased_matrix < 0] <- 0
  diag(noisy_biased_matrix) <- 0

  # Return all variants
  return(list(
    n1 = noisy_matrix1,
    n2 = noisy_matrix2,
    nb = noisy_biased_matrix
  ))
}


#' Save ggplot with white background
#' 
#' @description
#' Wrapper around ggplot2::ggsave that ensures white background.
#' This function masks ggplot2::ggsave.
#' @inheritParams ggplot2::ggsave
#' @export
ggsave <- function(..., bg = 'white') {
  ggplot2::ggsave(..., bg = bg)
}

#' Calculate Euclidean Distance Between Colors
#' 
#' @param x,y Numeric vectors representing coordinates in LAB color space
#' @return Numeric distance between colors
#' @keywords internal
euclidean_dist <- function(x, y) {
  sqrt(sum((x - y)^2))
}


#' Color Palettes
#' 
#' @description
#' Predefined color palettes optimized for visualization
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

#' @rdname color_palettes
#' @export
c25_claud <- c(
  "red", "green1", "blue1", "yellow3", "purple",
  "darkorange", "hotpink", "darkturquoise", "gold1",
  "maroon", "palegreen2", "dodgerblue2", "brown", "orchid1",
  "green4", "sandybrown", "steelblue4", "yellow4", "plum3",
  "darkorange4", "skyblue2", "deeppink1", "khaki2", "gray70"
)

#' @rdname color_palettes
#' @export
c25_old <- c(
  "purple", "darkturquoise", "gold1", "darkorange4",
  "green1", "blue1", "red", "plum3", "green4", "dodgerblue2",
  "darkorange", "black", "skyblue2", "hotpink", "palegreen2",
  "sandybrown", "gray70", "khaki2", "maroon", "orchid1",
  "yellow4", "deeppink1", "steelblue4", "yellow3", "brown"
)

#' @rdname color_palettes
#' @export
c25_older <- c(
  "dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00",
  "black", "gold1", "skyblue2", "#FB9A99", "palegreen2",
  "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon",
  "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise",
  "green1", "yellow4", "yellow3", "darkorange4", "brown"
)