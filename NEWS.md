# topolow 0.3.2

*  Initial release to CRAN (revised per CRAN reviewr's instructions).
*  Introduces the Topolow algorithm, a physics-inspired method for antigenic cartography.
*  Provides robust mapping and complete positioning of all antigens, even with highly sparse datasets (>95% missing values).
*  Implements automatic, likelihood-based estimation to determine the optimal dimensionality of the antigenic map.
*  Includes functionality to calculate "antigenic velocity" vectors to quantify the rate and direction of antigenic drift.
*  Features tools for handling and processing cross-reactivity and binding affinity assay data, including those with thresholded values.
*  Demonstrates improved prediction accuracy and run-to-run stability compared to traditional MDS methods.

# topolow 1.0.0 (2025-07-01)

* All exported methods now include `\value` documentation describing the output's class, structure, and meaning.
* Examples for unexported functions have been omitted, and `\dontrun{}` wrappers have been removed5. Slower examples are now wrapped in `\donttest{}` as appropriate.
* Functions no longer write to user directories by default. Functions where writing a file is the main purpose now require the user to specify an output directory.
* The complex distributed processing functionality has been removed, as it was not essential for typical use cases.
* The link to our paper and citation information have been updated.


# topolow 2.0.0 (2025-07-30)

The wizard function `Euclidify` was added to run all the workflow needed to get the main output automatically. 

## Deprecations

* `create_topolow_map()` is now deprecated in favor of `euclidean_embedding()`. The old function will be removed in version 3.0.0.
  - Parameter name changed: `distance_matrix` --> `dissimilarity_matrix`  
  - Function name changed: `create_topolow_map()` --> `euclidean_embedding()`

## Breaking Changes

* **`initial_parameter_optimization()`**: Parameter `distance_matrix` renamed to `dissimilarity_matrix` 
  - **Migration**: Replace `distance_matrix = your_matrix` with `dissimilarity_matrix = your_matrix`
* **`run_adaptive_sampling()`**: Parameter `distance_matrix` renamed to `dissimilarity_matrix`
  - **Migration**: Replace `distance_matrix = your_matrix` with `dissimilarity_matrix = your_matrix`
* **`adaptive_MC_sampling()`**: 
  - Parameter `distance_matrix` renamed to `dissimilarity_matrix` 
  - Removed parameter `batch_size` from `adaptive_MC_sampling()`; its value had no effect in the processes anyway
  - Removed parameter `num_parallel_jobs` from `run_adaptive_sampling`; set `max_cores` to define the number of cores and parallel jobs
  - **Migration**: Replace `distance_matrix = your_matrix` with `dissimilarity_matrix = your_matrix` and remove `batch_size` arguments
* **`create_cv_folds()`**: Parameter names and return structure changed
  - **Parameter changes**: `truth_matrix` --> `dissimilarity_matrix`, `no_noise_truth` --> `ground_truth_matrix`
  - **Return structure**: Now returns named list elements (`$truth`, `$train`) instead of indexed elements
  - **Migration**: Update parameter names and change `result[[1]][[1]]` to `result[[1]]$truth`, `result[[1]][[2]]` to `result[[1]]$train`
* `take_log` parameter in `clean_data()` is deprecated
  - Perform log transformation before calling these functions instead
  - Parameter will be removed in next major version
* **`analyze_network_structure()`**: Parameter `distance_matrix` renamed to `dissimilarity_matrix` for consistency with other functions
* **`calculate_diagnostics()`**: Return class changed from `topolow_amcs_diagnostics` to `topolow_diagnostics` for naming consistency
* **`plot_network_structure()`**: Removed `aesthetic_config` and `layout_config` parameters
  - **Migration**: Replace with `width`, `height`, `dpi` parameters
  - Fixed aesthetic values improve consistency but reduce customization
  - Added better handling for empty network cases
* **`scatterplot_fitted_vs_true()`**: Parameter names updated for consistency
  - **Migration**: `distance_matrix` --> `dissimilarity_matrix`, `p_dist_mat` --> `p_dissimilarity_mat`
  - **Migration**: Default `save_plot` changed from `TRUE` to `FALSE`
  - Improved modern ggplot2 syntax using `linewidth` instead of deprecated `size`
* **`error_calculator_comparison()`**: Parameter names changed for consistency
  - `p_dist_mat` --> `predicted_dissimilarities`
  - `truth_matrix` --> `true_dissimilarities` 
  - `input_matrix` --> `input_dissimilarities` (now optional, defaults to `NULL`)
  - **Migration**: Update all parameter names in function calls
* **`calculate_prediction_interval()`**: Parameter names changed for consistency  
  - `distance_matrix` --> `dissimilarity_matrix`
  - `p_dist_mat` --> `predicted_dissimilarity_matrix`
  - **Migration**: Update parameter names in function calls
* **`long_to_matrix`** was renamed to `titers_list_to_matrix` since it is specific to viral titer data processing.
* Function `process_antigenic_data` accepts a data frame as input, instead of the previous form of a file path.
* In `process_antigenic_data`, `is_titer` became `is_similarity` for clearity for broader audience. Parameter `id_prefix` was removed.

## New Features

* Added `euclidean_embedding()` function with enhanced performance and features:
  - **Matrix reordering**: Automatic spectral ordering concentrates largest dissimilarity values in corners for improved optimization
  - **Enhanced validation**: Better input data quality checks with informative warnings
  - **Improved documentation**: More detailed examples and parameter guidance

## Improvements

* Package dependencies where reduced from 20 to 13
* Enhanced algorithm documentation with clearer physics-inspired approach description
* Better handling of edge cases in dissimilarity matrix processing
* Improved error messages for parameter validation
* Updated `parameter_sensitivity` function to use modern ggplot2 syntax
* Improved input validation and error handling in sensitivity analysis
* Enhanced MLE calculation algorithm for better robustness
* Replaced deprecated `size` parameter with `linewidth` in plots
* Enhanced input validation and error messages in `create_cv_folds()`
* `input_dissimilarities` parameter now optional in `error_calculator_comparison()`
* `initial_parameter_optimization` saves/returns the parameters in log scale, consistent with other function
* A vignette was added

## Deprecation Timeline

* Version 2.0.0: `create_topolow_map()` deprecated, issues warning
* Version 3.0.0 (planned): `create_topolow_map()` will be removed

## Migration Guide

To update your code:

```r
# Old (deprecated):
result <- create_topolow_map(distance_matrix = my_matrix, 
  # ... other parameters
)

# New (recommended):
result <- euclidean_embedding(dissimilarity_matrix = my_matrix,  # parameter name changed
  # ... other parameters (unchanged)
)
```

# topolow 2.0.1 (2025-08-30)

Included figures in the vignette.


# topolow 2.1.0 (2025-10-02)

## New Features

### Subsampling for Computational Efficiency

* **Major Feature**: Added `opt_subsample` parameter to key optimization functions, enabling efficient parameter optimization on large datasets while maintaining final embedding quality.

* **New Functions**:
  - `check_matrix_connectivity()`: Validates that a dissimilarity matrix forms a connected graph
  - `subsample_dissimilarity_matrix()`: Creates random subsamples with automatic connectivity validation and adaptive size adjustment
  - `sanity_check_subsample()`: Validates subsample suitability for cross-validation

* **Enhanced Functions**:
  - `initial_parameter_optimization()`: Now accepts `opt_subsample` parameter
  - `run_adaptive_sampling()`: Now accepts `opt_subsample` parameter
  - `adaptive_MC_sampling()`: Now accepts `opt_subsample` parameter (internal)
  - `Euclidify()`: Now accepts `opt_subsample` parameter

### How Subsampling Works

When `opt_subsample` is specified:
1. Each parameter evaluation uses a random subsample of the specified size
2. Connectivity is automatically validated; disconnected subsamples are rejected
3. If connectivity fails, sample size needs to be increased
4. Different parameter evaluations use different subsamples for robustness
5. **Final embedding always uses the full dataset**

### Performance Benefits

- Speeds up parameter optimization by ~10-50x on large datasets (>500 points)
- Reduces memory usage proportional to subsample size
- Parameters found on subsamples generalize well to full data

### Breaking Changes

None. The `opt_subsample` parameter is optional (default: NULL = use full data).

### Other changes
- Function `prune_sparse_matrix()` is added to prune sparse dissimilarity matrices to a well-connected subset
- Package gridExtra is a required import now.

### Recommendations

- Datasets < 500 points: Use full data (`opt_subsample = NULL`)
- Datasets > 500 points: Recommended `opt_subsample = 200-500`
- Always ensure `opt_subsample >= folds` for reliable cross-validation

## Bug Fixes

- Conversion of matrices to numeric in "R/adaptive_sampling.R" are now properly handled by extract_numeric_values() function of the Topolow package.

## Improvements

* Enhanced connectivity checking using igraph
* Better error messages for disconnected data
* Adaptive strategies for handling sparse data
* Comprehensive logging of subsampling operations
* New diagnostic plots including MCMC exploration and parameter fit traces
