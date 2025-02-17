# Copyright (c) 2024 Omid Arhami omid.arhami@uga.edu
# License: free of charge access granted to any academic researcher to use this software for non-commercial, academic research purposes **only**.  Nobody may modify, distribute, sublicense, or publicly share the Software or any derivative works, until the paper is published by the original authors.  The Software is provided "as is" without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose and noninfringement.  In no event shall the authors or copyright holders be liable for any claim, damages or other liability, whether in an action of contract, tort or otherwise, arising from, out of or in connection with the Software or the use or other dealings in the Software.

# R/experiments.R

#' Experiment Running Functions
#' 
#' @description
#' Functions for running parameter optimization, comparison experiments,
#' and other computational experiments either locally or via SLURM.
#'
#' @keywords internal
"_PACKAGE"


#' Create Cross-validation Folds for Distance Matrix
#'
#' @description
#' Creates k-fold cross-validation splits of a distance matrix while maintaining 
#' symmetry. Each fold has a training matrix with some values masked for validation.
#'
#' @param truth_matrix Matrix of true distances 
#' @param no_noise_truth Optional matrix of noise-free distances. If provided, used as truth.
#' @param n_folds Integer number of folds to create
#' @param random_seed Integer random seed for reproducibility
#' @return List of lists, each containing:
#'   \item{truth}{Truth matrix for this fold} 
#'   \item{train}{Training matrix with masked validation entries}
#' @examples
#' \dontrun{
#' # Create 5-fold CV splits
#' folds <- create_cv_folds(dist_matrix, n_folds = 5, random_seed = 123)
#' }
#' @export
create_cv_folds <- function(truth_matrix, no_noise_truth = NULL, 
                          n_folds = 10, random_seed = NULL) {
  if (!is.null(random_seed)) {
    if (!is.numeric(random_seed) || random_seed != round(random_seed)) {
      stop("random_seed must be an integer")
    }
    set.seed(random_seed) 
  }

  # Validate truth matrices
  if (!is.matrix(truth_matrix)) {
    stop("truth_matrix must be a matrix")
  }
  
  if (!is.null(no_noise_truth)) {
    if (!is.matrix(no_noise_truth)) {
      stop("no_noise_truth must be NULL or a matrix")
    }
    if (!identical(dim(truth_matrix), dim(no_noise_truth))) {
      stop("truth_matrix and no_noise_truth must have same dimensions")
    }
  }
  
  # Validate n_folds
  if (!is.numeric(n_folds) || n_folds < 2 || n_folds != round(n_folds)) {
    stop("n_folds must be an integer >= 2")
  }
  
  if (n_folds > nrow(truth_matrix)) {
    stop("n_folds cannot be larger than number of rows in matrix")
  }

  # Use noise-free matrix as truth if provided
  if (!is.null(no_noise_truth)) {
    eval_truth <- no_noise_truth
  } else {
    eval_truth <- truth_matrix
  }

  # Initialize list to store folds
  matrix_list <- vector("list", n_folds)
  for(i in 1:n_folds) {
    matrix_list[[i]] <- list(eval_truth, NULL)
  }

  num_elements <- sum(!is.na(truth_matrix))
  holdout_size <- floor(num_elements/(n_folds * 2)) # Factor of 2 for symmetry
  
  # Create training copy to track available elements
  D_train <- truth_matrix
  
  # Create folds
  for(i in 1:n_folds) {
    # Sample validation indices from remaining elements
    random_indices <- sample(which(!is.na(D_train)), size = holdout_size)
    
    # Create training matrix for this fold
    input_matrix <- truth_matrix
    for(index in random_indices) {
      # Convert linear index to row/col
      row <- (index - 1) %/% nrow(truth_matrix) + 1  
      col <- (index - 1) %% ncol(truth_matrix) + 1
      
      # Mask validation entries
      input_matrix[row, col] <- NA
      input_matrix[col, row] <- NA # Maintain symmetry
    }
    
    # Store training matrix
    matrix_list[[i]][[2]] <- input_matrix
    
    # Update available elements for next fold
    for(index in random_indices) {
      row <- (index - 1) %/% nrow(D_train) + 1
      col <- (index - 1) %% ncol(D_train) + 1
      D_train[row, col] <- NA 
      D_train[col, row] <- NA
    }
  }

  return(matrix_list)
}



#' Create and Optimize RACMACS Map
#'
#' @description
#' Creates and optimizes an antigenic map using RACMACS and keeps the best optimization. 
#' This function wraps RACMACS functionality to provide a simplified interface for
#' map creation and optimization.
#'
#' @param titer_table Matrix or data frame of titer measurements
#' @param dim Integer number of dimensions for the map (default: 2)
#' @param optimization_number Integer number of optimization runs (default: 400)
#' @param scenario_name Character string for output file naming
#' @return RACMACS map object containing optimized coordinates
#' @examples
#' \dontrun{
#' # Create and optimize map from titer data
#' map <- create_and_optimize_RACMACS_map(titer_table)
#' 
#' # Create map with specific settings
#' map <- create_and_optimize_RACMACS_map(
#'   titer_table, 
#'   dim = 3,
#'   optimization_number = 1000,
#'   scenario_name = "example_map"
#' )
#' }
#' @importFrom Racmacs acmap optimizeMap keepBestOptimization
#' @importFrom parallel detectCores
#' @export
create_and_optimize_RACMACS_map <- function(titer_table, dim = 2,
                                          optimization_number = 400,
                                          scenario_name) {
  # Set number of cores for optimization
  options(RacOptimizer.num_cores = detectCores() - 4)
  
  # Create map object
  map <- acmap(titer_table = titer_table)
  
  # Optimize map
  map <- optimizeMap(
    options = list(ignore_disconnected = TRUE),
    map = map,
    number_of_dimensions = dim,
    number_of_optimizations = optimization_number,
    minimum_column_basis = "none"
  )
  
  # Keep best optimization
  map <- keepBestOptimization(map)
  
  # Save coordinates if scenario name provided
  if(!missing(scenario_name)) {
    save.coords(map, 
               filename = paste0(scenario_name, "_RACMACS_coords.csv"),
               optimization_number = 1,
               antigens = TRUE, 
               sera = TRUE)
  }
  
  return(map)
}



#' Generate Synthetic Distance Matrices with Missing Data
#'
#' @description
#' Creates synthetic distance matrices with controlled levels of missingness and noise
#' for testing and validating mapping algorithms. Generates multiple datasets with 
#' different dimensionalities and missingness patterns.
#'
#' @param n_dims_list Numeric vector of dimensions to generate data for
#' @param seeds Integer vector of random seeds (same length as n_dims_list)
#' @param n_points Integer number of points to generate
#' @param missingness_levels Named list of missingness percentages (default: list(S=0.67, M=0.77, L=0.87))
#' @param output_dir Character path to directory for saving outputs (optional)
#' @param prefix Character string to prefix output files (optional)
#' @param save_plots Logical whether to save network visualization plots
#' @return List containing:
#'   \item{matrices}{List of generated distance matrices}
#'   \item{panels}{List of generated assay panels}
#'   \item{metadata}{Data frame with generation parameters}
#' @examples
#' \dontrun{
#' # Generate datasets with different dimensions
#' results <- generate_synthetic_datasets(
#'   n_dims_list = c(2, 5, 10),
#'   seeds = c(123, 456, 789),
#'   n_points = 250,
#'   output_dir = "sim_data"
#' )
#' 
#' # Custom missingness levels
#' results <- generate_synthetic_datasets(
#'   n_dims_list = c(2, 5),
#'   seeds = c(123, 456),
#'   n_points = 200,
#'   missingness_levels = list(low=0.5, high=0.8)
#' )
#' }
#' @export
generate_synthetic_datasets <- function(n_dims_list, seeds, n_points,
                                     missingness_levels = list(S=0.67, M=0.77, L=0.87),
                                     output_dir = NULL,
                                     prefix = "sim",
                                     save_plots = FALSE) {
  
  # Input validation
  if (length(n_dims_list) != length(seeds)) {
    stop("n_dims_list and seeds must have the same length")
  }
  
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
  }
  
  # Validate dimensions
  if (!is.numeric(n_dims_list) || any(n_dims_list < 1) || 
      any(n_dims_list != round(n_dims_list))) {
    stop("All dimensions must be positive integers")
  }
  
  # Validate seeds
  if (!is.numeric(seeds) || any(seeds != round(seeds))) {
    stop("All seeds must be integers")
  }
  
  # Validate n_points
  if (!is.numeric(n_points) || n_points < 2 || n_points != round(n_points)) {
    stop("n_points must be an integer >= 2")
  }
  
  # Validate missingness levels
  if (!is.list(missingness_levels)) {
    stop("missingness_levels must be a list")
  }
  
  if (!all(unlist(missingness_levels) > 0 & unlist(missingness_levels) < 1)) {
    stop("All missingness levels must be between 0 and 1")
  }
  
  # Initialize return lists
  coordinates_list <- list()
  matrices_list <- list()
  panels_list <- list()
  metadata <- data.frame()
  
  # Process each dimension
  for (i in seq_along(n_dims_list)) {
    set.seed(seeds[i])
    n_dims <- n_dims_list[i]
    
    # Generate base dataset
    dataset <- generate_complex_data(n_points, n_dims, 
                                   n_clusters = 5, 
                                   cluster_spread = 0.1)
    
    # Add unique IDs
    dataset$unique_id <- generate_unique_string(nrow(dataset))
    
    # Generate distance matrix
    dist_matrix <- as.matrix(dist(dataset[, 1:n_dims]))
    rownames(dist_matrix) <- dataset$unique_id
    colnames(dist_matrix) <- dataset$unique_id
    
    # Select references
    n_ref <- floor(0.4 * n_points)
    selected_indices <- floor(seq(1, nrow(dataset), length.out = n_ref + 1))
    selected_indices <- selected_indices[1:n_ref]
    selected_names <- dataset$unique_id[selected_indices]
    
    # Label references and viruses
    virus_mask <- !(rownames(dist_matrix) %in% selected_names)
    dataset$unique_id[virus_mask] <- paste0("V/", dataset$unique_id[virus_mask])
    rownames(dist_matrix)[virus_mask] <- paste0("V/", rownames(dist_matrix)[virus_mask])
    colnames(dist_matrix)[virus_mask] <- paste0("V/", colnames(dist_matrix)[virus_mask])
    
    serum_mask <- rownames(dist_matrix) %in% selected_names
    dataset$unique_id[serum_mask] <- paste0("S/", dataset$unique_id[serum_mask])
    rownames(dist_matrix)[serum_mask] <- paste0("S/", rownames(dist_matrix)[serum_mask])
    colnames(dist_matrix)[serum_mask] <- paste0("S/", colnames(dist_matrix)[serum_mask])
    
    selected_names <- paste0("S/", selected_names)
    
    # Save the coordinates to return
    coordinates_list[[i]] <- dataset
    
    # Create base matrices
    matrices <- list()
    matrices$complete <- dist_matrix
    matrices$labelled <- only_virus_vs_as(dist_matrix, selected_names)
    
    # Generate matrices with different missingness levels
    for (level_name in names(missingness_levels)) {
      level_pct <- missingness_levels[[level_name]]
      missing_mat <- increase_na_percentage(matrices$labelled, level_pct)
      matrices[[paste0(level_name, "_missing")]] <- missing_mat
      
      # Add noise variations
      noisy <- add_noise_bias(missing_mat)
      matrices[[paste0(level_name, "_noise1")]] <- noisy$n1
      matrices[[paste0(level_name, "_noise2")]] <- noisy$n2
      matrices[[paste0(level_name, "_noise_bias")]] <- noisy$nb
    }
    
    # Generate panels
    panels <- list()
    panels$full <- symmetric_to_nonsymmetric_matrix(matrices$labelled, selected_names)
    
    for (level_name in names(missingness_levels)) {
      missing_mat <- matrices[[paste0(level_name, "_missing")]]
      panel <- symmetric_to_nonsymmetric_matrix(missing_mat, selected_names)
      panel <- panel[rowSums(!is.na(panel)) > 0, ]
      panels[[level_name]] <- panel
    }
    
    # Store results
    matrices_list[[i]] <- matrices
    panels_list[[i]] <- panels
    
    # Record metadata
    metadata <- rbind(metadata, data.frame(
      dimension = n_dims,
      seed = seeds[i],
      n_points = n_points,
      n_references = n_ref
    ))
    
    # Save files if output directory provided
    if (!is.null(output_dir)) {
      base_name <- file.path(output_dir, 
                            sprintf("%s_dim%d", prefix, n_dims))
      
      # Save matrices
      for (mat_name in names(matrices)) {
        saveRDS(matrices[[mat_name]], 
                paste0(base_name, "_", mat_name, ".rds"))
      }
      
      # Save panels
      for (panel_name in names(panels)) {
        saveRDS(panels[[panel_name]], 
                paste0(base_name, "_panel_", panel_name, ".rds"))
      }
      
      # Create network plots if requested
      if (save_plots) {
        for (level_name in names(missingness_levels)) {
          mat <- matrices[[paste0(level_name, "_missing")]]
          net <- analyze_network_structure(mat)
          plot_network_structure(net, 
            paste0(base_name, "_network_", level_name))
        }
      }
    }
  }
  
  # Return results
  return(list(
    coordinates = coordinates_list,
    matrices = matrices_list,
    panels = panels_list,
    metadata = metadata
  ))
}