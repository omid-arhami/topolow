# topolow 0.3.2

*  Initial public release on GitHub (pre-CRAN). The first CRAN release was 1.0.0.
*  Introduces the Topolow algorithm, a physics-inspired method for antigenic cartography.
*  Provides robust mapping and complete positioning of all antigens, even with highly sparse datasets (>95% missing values).
*  Implements automatic, likelihood-based estimation to determine the optimal dimensionality of the antigenic map.
*  Includes functionality to calculate "antigenic velocity" vectors to quantify the rate and direction of antigenic drift.
*  Features tools for handling and processing cross-reactivity and binding affinity assay data, including those with thresholded values.
*  Demonstrates improved prediction accuracy and run-to-run stability compared to traditional MDS methods.

# topolow 1.0.0 (2025-07-01)

* All exported methods now include `\value` documentation describing the output's class, structure, and meaning.
* Examples for unexported functions have been omitted, and `\dontrun{}` wrappers have been removed. Slower examples are now wrapped in `\donttest{}` as appropriate.
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


# topolow 2.1.0

## Deprecations

* `create_diagnostic_plots()` is deprecated in favor of `plot_mcmc_diagnostics()`,
  so that every plotting function shares the `plot_*` prefix. The old function
  still works and warns; it will be removed in version 3.0.0.
  - **Migration**: `create_diagnostic_plots(...)` --> `plot_mcmc_diagnostics(...)`
  - **Migration**: argument `res` --> `dpi`

## New Features

### Subsampling for Computational Efficiency

* **Major Feature**: Added the optional `opt_subsample` parameter to key optimization functions, enabling efficient parameter optimization on large datasets while maintaining final embedding quality. Parameter optimization still works reliably with subsampling, because likelihoods of samples of the same size are comparable, allowing us to choose the optimal parameter values.

* **New Functions**:
  - `check_matrix_connectivity()`: Validates that a dissimilarity matrix forms a connected graph
  - `subsample_dissimilarity_matrix()`: Creates random subsamples with automatic connectivity validation and adaptive size adjustment
  - `sanity_check_subsample()`: Validates subsample suitability for cross-validation
  - `prune_sparse_matrix()`: Prunes sparse dissimilarity matrices to a well-connected subset
  - `plot_mcmc_diagnostics()`: Replacement for `create_diagnostic_plots()` (see Deprecations)
  - `plot_ll_improvement()`, `plot_performance_trace()`: Diagnostics for the parameter search
  - `plot_euclidify_diagnostics()`, `create_diagnostic_report()`: Diagnostics for `Euclidify()` output

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

The `opt_subsample` parameter is optional (default: NULL = use full data).

### Performance Benefits

- Substantially reduces the cost of parameter optimization on large datasets
  (>500 points), since each evaluation embeds a subsample rather than the full
  matrix
- Reduces memory usage proportional to subsample size
- Parameters found on subsamples generalize well to full data

### Recommendations

- Datasets < 500 points: Use full data (`opt_subsample = NULL`)
- Datasets > 500 points: Recommended `opt_subsample = 200-500`
- Always ensure `opt_subsample >= folds` for reliable cross-validation

### Epoch-Based Parameter Search

* `initial_parameter_optimization()` gained an `epochs` argument (default: 1).
  With more than one epoch the search runs as a simple evolutionary strategy:

  1. **Initialization**: start from the user-provided parameter ranges.
  2. **Epoch loop**, for each epoch:
     a. Generate `num_samples` points by Latin hypercube sampling within the
        current ranges.
     b. If `opt_subsample` is set, each evaluation uses its own random subsample.
     c. Evaluate the parameter sets by cross-validation, in parallel batches.
     d. **Range update** (after every epoch but the last): sort by NLL, keep the
        top 50%, and set the next epoch's range to
        `0.75 * min(survivors)` through `1.25 * max(survivors)`. This lets the
        search drift toward and zoom in on promising regions.
  3. **Finalization**: log-transform the results of the final epoch for direct
     use with adaptive sampling.

### Other changes

- Package `gridExtra` is a required import now.

## Bug Fixes

- `initial_parameter_optimization()` (and therefore `Euclidify()`) could fail on
  Windows with `invalid connection` whenever the parameter search ran in more
  than one batch. A cluster shutdown handler was registered per batch, and all
  of them referred to the last cluster, so it was stopped repeatedly while the
  earlier clusters leaked. Each batch's cluster is now stopped as soon as that
  batch finishes. The failure was reliably triggered by `max_cores = 2`, which
  is the configuration used on CRAN check machines.
- Conversion of matrices to numeric in `R/adaptive_sampling.R` is now handled
  by the package's `extract_numeric_values()` function.
- `run_adaptive_sampling()` no longer sleeps for a fixed 1.5 seconds on every
  call. Two `Sys.sleep()` calls guarded file visibility on networked
  filesystems (NFS, Lustre); they now poll for the expected state instead, so
  the wait is effectively free on a local filesystem while networked ones still
  get up to 5 seconds. This removed about a second from every `Euclidify()`
  run and has no effect on results.

## Documentation

* Clarified that the `*_range` arguments of `Euclidify()` are a **starting
  range for the search, not a hard bound**. Initial sampling draws inside the
  range, but adaptive refinement samples from a kernel density estimate of the
  best parameter sets found so far, which can propose values outside it. So
  `ndim_range = c(3, 6)` does not guarantee a 3- to 6-dimensional result; only
  the sanity limits (for example `1 <= ndim <= 50`) are enforced. This
  behaviour is unchanged from earlier versions and is now documented, in a new
  section of `?Euclidify`. Use `euclidean_embedding()` with an explicit `ndim`
  when a hard cap is required.

## Improvements

* Enhanced connectivity checking using igraph
* Better error messages for disconnected data
* Adaptive strategies for handling sparse data
* Comprehensive logging of subsampling operations
* New diagnostic plots including MCMC exploration and parameter fit traces

## C++ Backend for Core Optimization

* **This is the first release of topolow containing compiled code.** The core
  optimization loop of `euclidean_embedding()` has been rewritten in C++ using
  Rcpp, which speeds up the embedding substantially, most visibly on larger
  matrices. The remaining R-level loops in that function were replaced with
  vectorized operations.

* **New parameter for `euclidean_embedding()`**:
  - `convergence_check_freq`: How often, in iterations, to test for convergence
    (default: 3). Larger values reduce the bookkeeping overhead; smaller values
    stop more promptly.

* **Implementation details**:
  - **COO format**: edge data is held as a coordinate list, which avoids the
    zero-dropping behaviour of sparse matrix classes
  - **Shuffled sweep order**: the sweep order is permuted each iteration with a
    C++ Mersenne Twister (`std::mt19937`), which is what lets the configuration
    escape local optima
  - **Immediate updates**: Gauss-Seidel style in-place position updates, as in
    the original R implementation
  - **Single-pass error calculation**: the MAE over active constraints is
    accumulated in one pass over the edge list at each convergence check
  - **Cache-friendly layout**: edge data is stored in contiguous arrays
  - **Pre-computed factors**: degree-based normalization factors are computed
    once, before the iteration begins
  - **Direct memory access**: the inner loop writes to the raw column-major
    position buffer

* **Return value enhancement**: the `convergence` field of the returned
  `topolow` object now includes:
  - `achieved`: whether convergence was reached
  - `error`: final MAE on active constraints
  - `final_k`: final spring constant after cooling

* **Dependencies**: added `Rcpp` to `LinkingTo` (compile-time only; no new
  runtime dependency). The backend is Rcpp-only and links no BLAS or LAPACK
  routines, so no `src/Makevars` is required on any platform.

### Reproducibility and numerical equivalence

* The algorithm, its force model, cooling schedule and convergence criterion are
  unchanged from 2.0.1. Results from this version are **statistically
  equivalent** to 2.0.1 but **not bit-identical**, for two reasons:
  - The sweep order is randomized by a C++ generator seeded from
    `std::random_device`, independent of R's RNG. `set.seed()` fixes the
    starting configuration (drawn with `stats::runif()`) but not the sweep
    order, so repeated runs differ slightly. The R implementation in 2.0.1
    shuffled the sweep order in the same way.
  - The convergence MAE (`result$convergence$error`) is now accumulated in
    package code rather than by a BLAS routine, so it can differ in the last
    representable digit from earlier development builds. This is far below the
    run-to-run variation produced by the shuffle, and it makes the value
    consistent across platforms, which it previously was not. Embedded
    positions are unaffected.
* See the new `Reproducibility` section of `?euclidean_embedding` for what is
  and is not seeded, and for how to obtain a single canonical map.

### Backward Compatibility

- All existing code using `euclidean_embedding()` will work without modification
- Default parameter values preserve the original algorithm's behavior
- Output structure remains compatible with downstream functions (`Euclidify()`, parameter optimization, etc.)


# topolow 2.2.0

## Breaking Changes

* `preserve_order` now defaults to `TRUE` throughout the package. The affected
  functions are `euclidean_embedding()`, `Euclidify()`,
  `initial_parameter_optimization()`, `run_adaptive_sampling()`,
  `adaptive_MC_sampling()`, `likelihood_function()` and
  `subsample_dissimilarity_matrix()`. `bootstrap_embedding()` already defaulted
  to `TRUE` and is unchanged.
* In practice, `euclidean_embedding()` no longer reorders the dissimilarity
  matrix internally by default, so the rows of `positions` and `est_distances`
  come back in the order they were supplied in, and
  `subsample_dissimilarity_matrix()` returns its `selected_indices` sorted
  rather than in draw order. The spectral reordering was a convergence-speed
  optimisation rather than part of the objective, but it does change which
  random initialisation each point receives, so coordinates, `mae` and `iter`
  from a seeded run will not reproduce those of the same call under 2.1.0.
  Matrices whose input order carries no meaning may also need more iterations
  to converge than before.
  - **Migration**: pass `preserve_order = FALSE` explicitly to restore the
    previous behaviour.

## New Features

### Bootstrapping over sera

Parameter optimization and map building can now resample one group of points
while holding another fixed. The motivating case is antigenic cartography:
keep every antigen and draw a random subset of the sera, so that neither the
fitted parameters nor the final map depend on one particular serum panel.
Only the sera are drawn because no pair of antigens is ever measured against
each other: the antigen-antigen block of the matrix is empty, and antigens
constrain one another only through the sera they were both titrated against.
Varying the subset of references while preserving all the antigens therefore
covers the majority of the uncertainty in positions due to the sample, though
not the variability contributed by the choice of antigen panel. See
`?bootstrap_embedding` for the full argument.

* **New function `bootstrap_embedding()`**: runs `euclidean_embedding()` over
  many replicate submatrices and aggregates the estimated distances cell by
  cell into a mean (`est_distances`), a standard deviation (`sd_distances`) and
  a replicate count (`n_obs`). It also returns a consensus embedding of the
  mean matrix, a per-block calibration table against the full-data map, and the
  per-replicate coordinates rigidly aligned to that map.
* **New function `identify_antigens_sera()`**: splits a dissimilarity matrix
  into antigens and sera, by default on the `"V/"` and `"S/"` prefixes written
  by `process_antigenic_data()`, with explicit vectors or custom regular
  expressions as overrides.
* **New function `bootstrap_blobs()`**: turns the aligned replicate coordinates
  into per-point polygons. Regions cannot be derived from `sd_distances`,
  because a distance matrix is invariant to rotation and so carries no
  directional information; the coordinates are needed for that.
* **`plot_temporal_mapping()` and `plot_cluster_mapping()`** gained
  `uncertainty`, `uncertainty_conf_level` and `uncertainty_points`. Passing a
  `bootstrap_embedding()` result draws a shaded region beneath each point.
  `uncertainty_points` chooses which points get one and defaults to
  `"antigens"`, because a region on every point at once makes a dense map hard
  to read; pass `"all"` or `"antisera"` for the others.
* **New S3 methods** `print()` and `summary()` for `topolow_bootstrap` objects.

### Subsampling a restricted set of points

* `subsample_dissimilarity_matrix()` gained `always_include` and `sample_pool`,
  which split a subsample into a fixed part and a randomly drawn part.
  `sample_size` keeps its meaning: the total number of points returned. With
  both arguments at their `NULL` default the function draws exactly the indices
  it drew before, from the same random number stream, and returns an identical
  object.
* `initial_parameter_optimization()`, `run_adaptive_sampling()`,
  `adaptive_MC_sampling()` and `Euclidify()` gained `opt_subsample_sera` and
  `opt_sera`. `opt_subsample_sera` is a count of sera, or a fraction of them
  when strictly between 0 and 1; every antigen is kept. It is mutually
  exclusive with `opt_subsample`.
* The subset size is a **fixed count**, resolved once per call and used for
  every evaluation and every replicate, so likelihoods stay comparable. The
  parameter-search output records it in two new columns, `opt_subsample` (total
  points) and `opt_n_sera` (how many of those were sera), so the same count can
  be read back and handed to `bootstrap_embedding(n_sera = )`. A fraction is
  only a convenience: it is converted to a count immediately and the count is
  what gets reported. The one thing that can shrink an individual replicate is
  `on_disconnect = "prune"`, and `replicates$n_points` shows when that happened.
  - **Migration**: none. All new arguments default to `NULL`.

## Bug Fixes

* `Euclidify()` did not pass `opt_subsample` to `run_adaptive_sampling()`, so
  the adaptive refinement stage silently embedded the full matrix even when
  subsampling was requested. The argument is now forwarded. Existing calls
  that set `opt_subsample` will find the adaptive stage faster, and its
  likelihood estimates correspondingly noisier, than before.
* `subsample_dissimilarity_matrix()` raised
  `invalid format '%d'; use format %s for character objects` whenever
  `sample_size >= nrow(dissimilarity_matrix)` and `verbose = TRUE`, because a
  literal was passed as an argument to `sprintf()` rather than as part of the
  format string.

## Improvements

* The package now calls `parallelly::availableCores()` where it previously
  called `future::availableCores()`, and `Imports` lists `parallelly` instead of
  `future`. This is a one-for-one swap rather than a new dependency:
  `availableCores()` was the only thing the package used from `future`, and
  `future` re-exported it from `parallelly`, so the function being called is
  unchanged (`identical(future::availableCores, parallelly::availableCores)` is
  `TRUE`). Core detection, including the scheduler-aware behaviour on HPC
  systems, is therefore unaffected. `future` 1.82.0 deprecated the re-export, so
  this also removes a deprecation warning from every call that sizes a cluster.
  Parallel execution itself has always used the `parallel` package and is
  untouched.
  - **Migration**: none, unless your own scripts relied on `topolow` pulling in
    `future` as a dependency. If you call `future::plan()` yourself, install
    `future` explicitly.
* `reduce_dimensions()` now attaches the fitted `center`, `rotation` and
  `scale_factor` of its principal-component reduction to the returned data
  frame, so that further points can be projected into the same plane. Callers
  that ignore the attributes are unaffected.
* The package now defines an `.onUnload` hook that calls
  `library.dynam.unload()`, so `unloadNamespace("topolow")` releases the
  compiled shared library instead of leaving it mapped in the session. On
  Windows this is what allows the package to be reinstalled over itself from an
  R session that has already loaded it; previously the open `topolow.dll` kept
  a lock on the file and the install aborted with
  `cannot remove earlier installation, is it in use?` and rolled back to the
  previously installed version. The compiled routines themselves are unchanged.
