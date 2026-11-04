# Copyright (c) 2026 Omid Arhami omid.arhami@uga.edu
# R/bootstrap.R

#' Bootstrap (subsampling) functions for the topolow package


# =============================================================================
# Internal helpers
# =============================================================================

#' Rigid Procrustes Alignment of Two Coordinate Sets
#'
#' @description
#' Aligns `X` onto `Y` by the orthogonal Procrustes (Kabsch) solution computed
#' on their shared row names: rotation, reflection and translation only. No
#' scaling is fitted, because topolow's embedding scale is anchored by the
#' measured dissimilarities and fitting a dilation would absorb the very
#' design effect the bootstrap is meant to measure.
#'
#' @param X Numeric matrix of coordinates to align, with row names.
#' @param Y Numeric matrix of reference coordinates, with row names.
#'
#' @return A numeric matrix the same shape as `X`, with the rigid transform
#'   applied to every row of `X`.
#' @keywords internal
procrustes_align <- function(X, Y) {
  shared <- intersect(rownames(X), rownames(Y))
  if (length(shared) < 2) {
    stop("X and Y must share at least two named points")
  }
  A <- X[shared, , drop = FALSE]
  B <- Y[shared, , drop = FALSE]
  cA <- colMeans(A)
  cB <- colMeans(B)
  sv <- svd(crossprod(sweep(A, 2, cA), sweep(B, 2, cB)))
  rot <- sv$u %*% t(sv$v)
  out <- sweep(sweep(X, 2, cA) %*% rot, 2, cB, "+")
  rownames(out) <- rownames(X)
  out
}


#' Connected Components of a Dissimilarity Matrix
#'
#' @param m Square dissimilarity matrix.
#' @return The `igraph` components object for the observation graph.
#' @keywords internal
matrix_components <- function(m) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for connectivity checking.")
  }
  adj <- !is.na(m)
  diag(adj) <- FALSE
  g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected")
  igraph::components(g)
}


#' Draw One Bootstrap Replicate Submatrix
#'
#' @param m Square dissimilarity matrix.
#' @param always_idx Integer vector of row indices kept in every replicate.
#' @param pool_idx Integer vector of row indices the random draw comes from.
#' @param n_draw Integer. How many points to draw from `pool_idx`.
#' @param preserve_order Logical. Sort the selected indices before subsetting.
#' @param on_disconnect Character. One of "prune", "retry", "error".
#' @param max_attempts Integer. Re-draws allowed when `on_disconnect = "retry"`.
#'
#' @return A list with `matrix`, `pruned` (character vector of dropped names)
#'   and `attempts`, or NULL when a connected draw could not be obtained.
#' @keywords internal
draw_replicate_matrix <- function(m, always_idx, pool_idx, n_draw,
                                  preserve_order, on_disconnect, max_attempts) {
  n_try <- if (identical(on_disconnect, "retry")) max_attempts else 1L

  for (attempt in seq_len(n_try)) {
    idx <- c(always_idx, pool_idx[sample.int(length(pool_idx), size = n_draw)])
    if (preserve_order) idx <- sort(idx)
    sub <- m[idx, idx, drop = FALSE]

    comp <- matrix_components(sub)
    if (comp$no == 1) {
      return(list(matrix = sub, pruned = character(0), attempts = attempt))
    }

    if (identical(on_disconnect, "error")) {
      stop(sprintf("replicate submatrix is disconnected (%d components)", comp$no))
    }
    if (identical(on_disconnect, "prune")) {
      keep <- comp$membership == which.max(comp$csize)
      if (sum(keep) < 3) {
        return(NULL)
      }
      return(list(matrix = sub[keep, keep, drop = FALSE],
                  pruned = rownames(sub)[!keep],
                  attempts = attempt))
    }
  }
  NULL
}


#' Per-block Comparison of Two Distance Matrices
#'
#' @param boot Numeric matrix of bootstrap mean distances.
#' @param full Numeric matrix of full-data distances, same dimnames as `boot`.
#' @param antigens Character vector of antigen names.
#' @param sera Character vector of serum names.
#'
#' @return A data frame with one row per block, or NULL when no block has
#'   enough comparable cells.
#' @keywords internal
calibration_table <- function(boot, full, antigens, sera) {
  blocks <- list(
    list(name = "antigen-antigen", r = antigens, c = antigens, tri = TRUE),
    list(name = "antigen-serum",   r = antigens, c = sera,     tri = FALSE),
    list(name = "serum-serum",     r = sera,     c = sera,     tri = TRUE)
  )
  rows <- lapply(blocks, function(b) {
    if (length(b$r) < 2 || length(b$c) < 2) return(NULL)
    xb <- boot[b$r, b$c, drop = FALSE]
    xf <- full[b$r, b$c, drop = FALSE]
    if (b$tri) {
      keep <- upper.tri(xb)
      y <- xb[keep]; x <- xf[keep]
    } else {
      y <- as.vector(xb); x <- as.vector(xf)
    }
    ok <- is.finite(x) & is.finite(y) & x > 0
    x <- x[ok]; y <- y[ok]
    if (length(x) < 3) return(NULL)
    fit <- stats::lm(y ~ x)
    data.frame(
      block = b$name,
      n_cells = length(x),
      slope0 = sum(x * y) / sum(x^2),
      intercept = unname(stats::coef(fit)[1]),
      slope = unname(stats::coef(fit)[2]),
      r2 = summary(fit)$r.squared,
      ratio_of_medians = stats::median(y) / stats::median(x),
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}


# =============================================================================
# Main function
# =============================================================================

#' Bootstrap a Euclidean Embedding over a Resampled Subset of Points
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Runs [euclidean_embedding()] repeatedly. Each replicate embeds a submatrix
#' that keeps a fixed set of points and draws the rest at random, and the
#' estimated distances are then aggregated cell by cell into a mean, a standard
#' deviation and a replicate count.
#'
#' In antigenic cartography the fixed set is every antigen and the drawn set is
#' a random subset of the sera, which quantifies how much of a map's geometry
#' depends on which sera happened to be in the panel.
#'
#' @details
#' **The scheme is subsampling, not resampling with replacement.** Drawing sera
#' with replacement would place duplicated rows in the submatrix, giving a pair
#' of points a measured dissimilarity of zero, which the force model treats as a
#' rest-length-zero spring, and would distort the degree normalisation of every
#' antigen measured against that serum. Drawing *m* of *n* points without
#' replacement is *m*-out-of-*n* subsampling, equivalently a delete-*d*
#' jackknife with *d = n - m*.
#'
#' **The subset size is a fixed count.** `frac_sera` is a convenience that is
#' converted to a count once, before any replicate runs; `n_sera` sets that
#' count directly and overrides it. Every replicate then draws exactly that many
#' points, which is what keeps the replicates comparable with each other and
#' with the parameter search. To use the same size the search used, read
#' `opt_n_sera` from its output and pass it as `n_sera`. The only thing that can
#' shrink an individual replicate is `on_disconnect = "prune"`, and
#' `replicates$n_points` records it when it happens.
#'
#' **Aggregation.** A cell of the output is the mean over the replicates in
#' which both of its endpoints were present. Antigen-antigen cells appear in
#' every replicate, antigen-serum cells in about `frac_sera` of them, and
#' serum-serum cells in about `frac_sera^2` of them, so a useful rule is
#' `n_replicates >= 30 / frac_sera^2`. Cells seen fewer than `min_replicates`
#' times are `NA` in both `est_distances` and `sd_distances`.
#'
#' **Why distances and not coordinates.** Pairwise distances are invariant to
#' rotation, reflection and translation, so they can be averaged directly. The
#' embedding scale is anchored because spring rest lengths are the measured
#' dissimilarities in their own units, so no rescaling is applied before
#' averaging. Coordinates are *not* frame-free: the per-replicate positions
#' returned in `replicate_positions` are rigidly aligned to the full-data map
#' with [procrustes_align()] before being stored.
#'
#' **Blocks are not comparable.** Antigen-antigen and serum-serum cells are
#' never observed in any replicate, so their spread measures the stability of an
#' extrapolation. Antigen-serum cells have a spring behind them whenever that
#' serum was drawn. Summarise the three blocks separately rather than putting
#' them on one scale.
#'
#' @section What the standard deviations do and do not mean:
#' The reported standard deviation mixes two sources. One is the variability of
#' interest, caused by drawing a different subset of points. The other is
#' topolow's own run-to-run variability: the optimizer reshuffles its sweep
#' order with a generator independent of R's RNG, so repeated runs on an
#' identical matrix do not give identical coordinates (see the Reproducibility
#' section of [euclidean_embedding()]). Running with `frac_sera = 1` isolates
#' that second component, since then no point is ever dropped.
#'
#' The raw standard deviation is also not the sampling standard deviation of the
#' full-data estimate. For that interpretation multiply by
#' `settings$jackknife_factor`, which is `sqrt(f / (1 - f))` for a retained
#' fraction `f` (Shao and Wu, 1989); it equals exactly 1 at `frac_sera = 0.5`.
#' The correction is reported, never applied.
#'
#' Finally, nothing here represents uncertainty in the antigen panel, or
#' measurement error in the assay values themselves.
#'
#' @section Reproducibility:
#' `random_seed` fixes which points each replicate draws, and therefore `n_obs`
#' and `point_counts`. It does not fix the embeddings, so `est_distances` and
#' `sd_distances` differ slightly between runs even under `set.seed()`.
#'
#' @param dissimilarity_matrix Matrix. A square, symmetric dissimilarity matrix.
#'        May contain NA for missing measurements and character strings with
#'        `<` or `>` prefixes for thresholded measurements. Row names must be
#'        unique; when absent, points are named `P1`, `P2`, ... internally.
#' @param ndim Integer. Number of dimensions for the embedding space.
#' @param k0 Numeric. Initial spring constant controlling spring forces.
#' @param cooling_rate Numeric. Rate of spring constant decay per iteration.
#' @param c_repulsion Numeric. Repulsion constant controlling repulsive forces.
#' @param mapping_max_iter Integer. Maximum map optimization iterations per
#'        replicate. Default: 1000.
#' @param relative_epsilon Numeric. Convergence threshold. Default: 1e-4.
#' @param convergence_counter Integer. Consecutive iterations below the
#'        threshold before declaring convergence. Default: 5.
#' @param n_replicates Integer. Number of bootstrap replicates. Default: 100.
#' @param frac_sera Numeric in (0, 1]. Fraction of the drawable points (the
#'        sera, by default) retained in each replicate. Ignored when `n_sera` is
#'        supplied. Default: 0.8.
#' @param n_sera Integer or NULL. Exact number of drawable points retained in
#'        each replicate. Overrides `frac_sera`. Default: NULL.
#' @param sera Character or numeric vector, or NULL. The sera, as row names or
#'        row indices. NULL detects them with [identify_antigens_sera()].
#'        Ignored when `sample_pool` is supplied.
#' @param always_include Character or numeric vector, or NULL. Points kept in
#'        every replicate. NULL means every point that is not a serum.
#' @param sample_pool Character or numeric vector, or NULL. Points the random
#'        draw is taken from. NULL means the sera.
#' @param on_disconnect Character. What to do when a replicate's submatrix is
#'        disconnected: `"prune"` (default) keeps its largest connected
#'        component and records the dropped points, `"retry"` re-draws up to
#'        `max_attempts` times, `"error"` fails the replicate.
#' @param max_attempts Integer. Re-draws allowed when `on_disconnect = "retry"`.
#'        Default: 5.
#' @param min_replicates Integer. Minimum number of replicates a cell must
#'        appear in before its mean and standard deviation are reported.
#'        Default: 2.
#' @param consensus Logical. Whether to embed the mean distance matrix to obtain
#'        consensus coordinates. Default: TRUE.
#' @param keep_replicates Logical. Whether to retain the aligned per-replicate
#'        coordinates, which [bootstrap_blobs()] and the uncertainty overlays of
#'        [plot_temporal_mapping()] need. Default: TRUE.
#' @param max_cores Integer or NULL. Maximum cores. Default: NULL, meaning
#'        `parallelly::availableCores() - 1`.
#' @param random_seed Integer or NULL. Seeds the point draws. Default: NULL.
#' @param write_files Logical. Whether to write the aggregates and the replicate
#'        table as CSV. Default: FALSE.
#' @param output_dir Character. Required when `write_files` is TRUE.
#' @param verbose Logical. Whether to print progress messages. Default: FALSE.
#' @param preserve_order Logical. Passed to [euclidean_embedding()]. Default:
#'        TRUE, which keeps each replicate's points in input order.
#'
#' @return A `list` of class `topolow_bootstrap` containing:
#' \itemize{
#'   \item `est_distances`: Matrix of per-cell mean estimated distances.
#'   \item `sd_distances`: Matrix of per-cell standard deviations.
#'   \item `n_obs`: Integer matrix of replicates contributing to each cell.
#'   \item `positions`: Consensus coordinates, or NULL.
#'   \item `consensus_map`: The `topolow` object behind `positions`, or NULL.
#'   \item `replicate_positions`: List of aligned per-replicate coordinate
#'     matrices with NA rows for absent points, or NULL.
#'   \item `full_map`: The single full-data `topolow` embedding.
#'   \item `calibration`: Data frame comparing bootstrap means with the
#'     full-data distances, per block.
#'   \item `point_counts`: Named integer vector of replicates each point
#'     survived.
#'   \item `replicates`: Data frame with one row per replicate.
#'   \item `mae`: Mean absolute error of `est_distances` against the observed
#'     dissimilarities.
#'   \item `groups`: The fixed and drawable name vectors actually used.
#'   \item `parameters`: List of the embedding parameters.
#'   \item `settings`: List of the bootstrap settings, including
#'     `jackknife_factor`.
#' }
#'
#' @references
#' Shao, J. and Wu, C. F. J. (1989). A general theory for jackknife variance
#' estimation. *The Annals of Statistics*, 17(3), 1176-1197.
#'
#' @examples
#' # Four antigens and three sera; two sera are drawn per replicate.
#' nms <- c("V/a1", "V/a2", "V/a3", "V/a4",
#'          "S/s1", "S/s2", "S/s3", "S/s4", "S/s5")
#' m <- as.matrix(dist(matrix(rnorm(18), nrow = 9)))
#' dimnames(m) <- list(nms, nms)
#'
#' res <- bootstrap_embedding(
#'   dissimilarity_matrix = m, ndim = 2,
#'   k0 = 1.0, cooling_rate = 0.01, c_repulsion = 0.01,
#'   mapping_max_iter = 30, n_replicates = 6, n_sera = 4,
#'   consensus = FALSE, max_cores = 1
#' )
#' res$est_distances[1:3, 1:3]
#' res$n_obs[1:3, 1:3]
#'
#' @seealso
#' [euclidean_embedding()] for a single embedding,
#' [identify_antigens_sera()] for naming the two groups,
#' [bootstrap_blobs()] for turning the replicates into uncertainty regions.
#'
#' @importFrom stats lm median coef setNames
#' @importFrom parallelly availableCores
#' @importFrom parallel makeCluster stopCluster clusterExport clusterEvalQ parLapply mclapply
#' @importFrom utils write.csv
#' @export
bootstrap_embedding <- function(dissimilarity_matrix,
                                ndim,
                                k0,
                                cooling_rate,
                                c_repulsion,
                                mapping_max_iter = 1000,
                                relative_epsilon = 1e-4,
                                convergence_counter = 5,
                                n_replicates = 100,
                                frac_sera = 0.8,
                                n_sera = NULL,
                                sera = NULL,
                                always_include = NULL,
                                sample_pool = NULL,
                                on_disconnect = c("prune", "retry", "error"),
                                max_attempts = 5,
                                min_replicates = 2,
                                consensus = TRUE,
                                keep_replicates = TRUE,
                                max_cores = NULL,
                                random_seed = NULL,
                                write_files = FALSE,
                                output_dir,
                                verbose = FALSE,
                                preserve_order = TRUE) {

  on_disconnect <- match.arg(on_disconnect)

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------
  if (!is.matrix(dissimilarity_matrix)) {
    stop("dissimilarity_matrix must be a matrix")
  }
  if (nrow(dissimilarity_matrix) != ncol(dissimilarity_matrix)) {
    stop("dissimilarity_matrix must be square")
  }
  if (nrow(dissimilarity_matrix) < 3) {
    stop("dissimilarity_matrix must have at least 3 points")
  }
  if (!is.numeric(ndim) || length(ndim) != 1 || ndim < 1 || ndim != round(ndim)) {
    stop("ndim must be a positive integer")
  }
  if (!is.numeric(k0) || k0 <= 0) {
    stop("k0 must be a positive number")
  }
  if (!is.numeric(cooling_rate) || cooling_rate <= 0 || cooling_rate >= 1) {
    stop("cooling_rate must be between 0 and 1")
  }
  if (!is.numeric(c_repulsion) || c_repulsion <= 0) {
    stop("c_repulsion must be a positive number")
  }
  if (!is.numeric(n_replicates) || length(n_replicates) != 1 ||
      n_replicates < 2 || n_replicates != round(n_replicates)) {
    stop("n_replicates must be an integer greater than or equal to 2")
  }
  if (!is.numeric(min_replicates) || length(min_replicates) != 1 ||
      min_replicates < 2 || min_replicates != round(min_replicates)) {
    stop("min_replicates must be an integer greater than or equal to 2")
  }
  if (!is.numeric(frac_sera) || length(frac_sera) != 1 ||
      frac_sera <= 0 || frac_sera > 1) {
    stop("frac_sera must be a number greater than 0 and at most 1")
  }
  if (write_files && missing(output_dir)) {
    stop("An 'output_dir' must be provided when 'write_files' is TRUE.", call. = FALSE)
  }

  n_replicates <- as.integer(n_replicates)
  min_replicates <- as.integer(min_replicates)

  # Names are the key the aggregation is built on
  if (is.null(rownames(dissimilarity_matrix))) {
    auto_names <- paste0("P", seq_len(nrow(dissimilarity_matrix)))
    rownames(dissimilarity_matrix) <- auto_names
    colnames(dissimilarity_matrix) <- auto_names
    if (verbose) {
      cat("dissimilarity_matrix had no row names; using P1..Pn\n")
    }
  }
  if (anyDuplicated(rownames(dissimilarity_matrix))) {
    stop("dissimilarity_matrix must have unique row names")
  }

  all_names <- rownames(dissimilarity_matrix)
  n <- length(all_names)

  # ---------------------------------------------------------------------------
  # Resolve the fixed and drawable groups
  # ---------------------------------------------------------------------------
  if (is.null(always_include) || is.null(sample_pool)) {
    groups <- identify_antigens_sera(dissimilarity_matrix, sera = sera)
    if (is.null(sample_pool)) sample_pool <- groups$sera
    if (is.null(always_include)) always_include <- groups$antigens
  }
  always_idx <- resolve_point_indices(always_include, dissimilarity_matrix,
                                      "always_include")
  pool_idx <- resolve_point_indices(sample_pool, dissimilarity_matrix,
                                    "sample_pool")
  pool_idx <- setdiff(pool_idx, always_idx)
  if (length(pool_idx) < 2) {
    stop("sample_pool must contain at least 2 points that are not in always_include")
  }

  n_draw <- if (!is.null(n_sera)) {
    if (!is.numeric(n_sera) || length(n_sera) != 1 || n_sera < 1) {
      stop("n_sera must be NULL or a positive integer")
    }
    min(as.integer(n_sera), length(pool_idx))
  } else {
    max(1L, min(as.integer(round(frac_sera * length(pool_idx))), length(pool_idx)))
  }
  retained_fraction <- n_draw / length(pool_idx)
  if (n_draw == length(pool_idx)) {
    warning("every drawable point is retained in each replicate; replicates differ only by optimizer randomness")
  }

  if (verbose) {
    cat(sprintf("Bootstrap: %d replicates, %d fixed points, %d of %d drawn (f = %.2f)\n",
                n_replicates, length(always_idx), n_draw, length(pool_idx),
                retained_fraction))
    cat(sprintf("  Rule of thumb for stable serum-serum cells: n_replicates >= %d\n",
                ceiling(30 / retained_fraction^2)))
  }

  # ---------------------------------------------------------------------------
  # Full-data embedding: calibration reference and Procrustes target
  # ---------------------------------------------------------------------------
  if (verbose) cat("Embedding the full matrix (reference)...\n")
  full_map <- euclidean_embedding(
    dissimilarity_matrix = dissimilarity_matrix,
    ndim = ndim, mapping_max_iter = mapping_max_iter,
    k0 = k0, cooling_rate = cooling_rate, c_repulsion = c_repulsion,
    relative_epsilon = relative_epsilon,
    convergence_counter = convergence_counter,
    verbose = FALSE, preserve_order = preserve_order
  )
  ref_positions <- full_map$positions

  # ---------------------------------------------------------------------------
  # Replicate worker
  # ---------------------------------------------------------------------------
  run_replicate <- function(b) {
    seed_b <- if (!is.null(random_seed)) as.integer(random_seed) + b else NA_integer_
    if (!is.na(seed_b)) set.seed(seed_b)

    drawn <- tryCatch(
      draw_replicate_matrix(dissimilarity_matrix, always_idx, pool_idx, n_draw,
                            preserve_order, on_disconnect, max_attempts),
      error = function(e) e
    )
    if (inherits(drawn, "error")) {
      return(list(status = "subsample_failed", seed = seed_b,
                  message = substr(conditionMessage(drawn), 1, 100)))
    }
    if (is.null(drawn)) {
      return(list(status = "subsample_failed", seed = seed_b,
                  message = "no connected submatrix obtained"))
    }

    fit <- tryCatch(
      euclidean_embedding(
        dissimilarity_matrix = drawn$matrix,
        ndim = ndim, mapping_max_iter = mapping_max_iter,
        k0 = k0, cooling_rate = cooling_rate, c_repulsion = c_repulsion,
        relative_epsilon = relative_epsilon,
        convergence_counter = convergence_counter,
        verbose = FALSE, preserve_order = preserve_order
      ),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      return(list(status = "embedding_failed", seed = seed_b,
                  message = substr(conditionMessage(fit), 1, 100)))
    }

    list(
      status = "ok",
      seed = seed_b,
      dist = fit$est_distances,
      pos = fit$positions,
      n_points = nrow(drawn$matrix),
      n_pruned = length(drawn$pruned),
      mae = fit$mae,
      iter = fit$iter,
      converged = isTRUE(fit$convergence$achieved),
      message = NA_character_
    )
  }

  # ---------------------------------------------------------------------------
  # Core sizing and cluster setup (one cluster for the whole run)
  # ---------------------------------------------------------------------------
  available_cores <- parallelly::availableCores()
  cores_to_use <- if (is.null(max_cores)) {
    max(1L, as.integer(available_cores) - 1L)
  } else {
    max(1L, min(as.integer(max_cores), as.integer(available_cores)))
  }
  cores_to_use <- min(cores_to_use, n_replicates)

  cl <- NULL
  if (cores_to_use > 1 && .Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(cores_to_use)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, library(topolow))
    parallel::clusterExport(
      cl,
      c("dissimilarity_matrix", "always_idx", "pool_idx", "n_draw", "ndim",
        "k0", "cooling_rate", "c_repulsion", "mapping_max_iter",
        "relative_epsilon", "convergence_counter", "preserve_order",
        "on_disconnect", "max_attempts", "random_seed"),
      envir = environment()
    )
  }

  # ---------------------------------------------------------------------------
  # Accumulators
  # ---------------------------------------------------------------------------
  dn <- list(all_names, all_names)
  sum_mat <- matrix(0, n, n, dimnames = dn)
  sumsq_mat <- matrix(0, n, n, dimnames = dn)
  cnt_mat <- matrix(0L, n, n, dimnames = dn)
  point_counts <- stats::setNames(integer(n), all_names)
  clouds <- if (keep_replicates) vector("list", n_replicates) else NULL
  rep_rows <- vector("list", n_replicates)

  batch_size <- min(max(cores_to_use * 4L, 1L), n_replicates)
  batches <- split(seq_len(n_replicates),
                   ceiling(seq_len(n_replicates) / batch_size))

  for (batch in batches) {
    out <- if (!is.null(cl)) {
      parallel::parLapply(cl, batch, run_replicate)
    } else if (cores_to_use > 1) {
      parallel::mclapply(batch, run_replicate, mc.cores = cores_to_use)
    } else {
      lapply(batch, run_replicate)
    }

    for (j in seq_along(batch)) {
      b <- batch[j]
      r <- out[[j]]
      if (is.null(r) || !identical(r$status, "ok")) {
        rep_rows[[b]] <- data.frame(
          replicate = b,
          seed = if (is.null(r)) NA_integer_ else r$seed,
          n_points = NA_integer_, n_sera_drawn = n_draw, n_pruned = NA_integer_,
          mae = NA_real_, iter = NA_integer_, converged = NA,
          status = if (is.null(r)) "worker_failed" else r$status,
          message = if (is.null(r)) NA_character_ else r$message,
          stringsAsFactors = FALSE
        )
        next
      }

      d <- r$dist
      idx <- match(rownames(d), all_names)
      if (anyNA(idx)) {
        stop("replicate returned point names absent from the input matrix")
      }
      sum_mat[idx, idx] <- sum_mat[idx, idx] + d
      sumsq_mat[idx, idx] <- sumsq_mat[idx, idx] + d^2
      cnt_mat[idx, idx] <- cnt_mat[idx, idx] + 1L
      point_counts[idx] <- point_counts[idx] + 1L

      if (keep_replicates) {
        cloud <- matrix(NA_real_, n, ndim, dimnames = list(all_names, NULL))
        cloud[idx, ] <- procrustes_align(r$pos, ref_positions)
        clouds[[b]] <- cloud
      }

      rep_rows[[b]] <- data.frame(
        replicate = b, seed = r$seed, n_points = r$n_points,
        n_sera_drawn = n_draw, n_pruned = r$n_pruned,
        mae = r$mae, iter = r$iter, converged = r$converged,
        status = "ok", message = NA_character_,
        stringsAsFactors = FALSE
      )
    }
    rm(out)
    gc(verbose = FALSE)
    if (verbose) {
      cat(sprintf("  %d/%d replicates done\n", max(batch), n_replicates))
    }
  }

  replicates <- do.call(rbind, rep_rows)
  n_successful <- sum(replicates$status == "ok")

  if (n_successful == 0) {
    tab <- table(replicates$status)
    stop(sprintf("all %d bootstrap replicates failed (%s)",
                 n_replicates,
                 paste(sprintf("%s: %d", names(tab), as.integer(tab)),
                       collapse = "; ")))
  }
  if (n_successful < min_replicates) {
    stop(sprintf("only %d of %d replicates succeeded, which is fewer than min_replicates (%d)",
                 n_successful, n_replicates, min_replicates))
  }
  if (n_successful < n_replicates) {
    warning(sprintf("%d of %d bootstrap replicates failed; aggregates use %d replicates",
                    n_replicates - n_successful, n_replicates, n_successful))
  }
  if (keep_replicates) {
    clouds <- clouds[!vapply(clouds, is.null, logical(1))]
  }

  # ---------------------------------------------------------------------------
  # Finalise the aggregates
  # ---------------------------------------------------------------------------
  reported <- cnt_mat >= min_replicates
  safe_cnt <- pmax(cnt_mat, 1L)
  mean_mat <- sum_mat / safe_cnt
  var_mat <- (sumsq_mat - sum_mat^2 / safe_cnt) / pmax(cnt_mat - 1L, 1L)

  est_distances <- ifelse(reported, mean_mat, NA_real_)
  sd_distances <- ifelse(reported, sqrt(pmax(var_mat, 0)), NA_real_)
  dimnames(est_distances) <- dn
  dimnames(sd_distances) <- dn

  groups_used <- list(
    antigens = all_names[sort(always_idx)],
    sera = all_names[sort(pool_idx)]
  )

  if (any(cnt_mat == 0)) {
    warning("some point pairs were never sampled together; est_distances has NA entries")
  }
  if (any(point_counts == 0)) {
    warning(sprintf("%d points appeared in no replicate", sum(point_counts == 0)))
  }
  # Pruning, not the random draw, is what biases which replicates a point
  # contributes to. A drawn point appearing in about `retained_fraction` of the
  # replicates is the design working as intended, so only warn once something
  # was actually pruned, and then only for points that fall well short of what
  # the design implies for them.
  if (any(replicates$n_pruned > 0, na.rm = TRUE)) {
    expected <- ifelse(all_names %in% groups_used$antigens, 1, retained_fraction)
    short <- point_counts / n_successful < 0.8 * expected
    if (any(short)) {
      warning(sprintf("%d points were pruned from enough replicates to fall short of the sampling design; their standard deviations are conditional and understate uncertainty",
                      sum(short)))
    }
  }

  # MAE against the observed dissimilarities, matching euclidean_embedding()
  obs <- suppressWarnings(as.numeric(dissimilarity_matrix))
  fitted_vals <- as.vector(est_distances)
  ok <- !is.na(obs) & !is.na(fitted_vals)
  mae <- if (any(ok)) mean(abs(obs[ok] - fitted_vals[ok])) else NA_real_

  calibration <- calibration_table(est_distances,
                                   full_map$est_distances[all_names, all_names],
                                   groups_used$antigens, groups_used$sera)

  # ---------------------------------------------------------------------------
  # Consensus map
  # ---------------------------------------------------------------------------
  consensus_map <- NULL
  if (consensus) {
    if (verbose) cat("Embedding the mean distance matrix (consensus)...\n")
    consensus_map <- euclidean_embedding(
      dissimilarity_matrix = est_distances,
      ndim = ndim, mapping_max_iter = mapping_max_iter,
      k0 = k0, cooling_rate = cooling_rate, c_repulsion = c_repulsion,
      relative_epsilon = relative_epsilon,
      convergence_counter = convergence_counter,
      verbose = FALSE, preserve_order = preserve_order
    )
    consensus_map$positions <- procrustes_align(consensus_map$positions,
                                                ref_positions)
  }

  jackknife_factor <- if (retained_fraction < 1) {
    sqrt(retained_fraction / (1 - retained_fraction))
  } else {
    NA_real_
  }

  result <- structure(
    list(
      est_distances = est_distances,
      sd_distances = sd_distances,
      n_obs = cnt_mat,
      positions = if (is.null(consensus_map)) NULL else consensus_map$positions,
      consensus_map = consensus_map,
      replicate_positions = clouds,
      full_map = full_map,
      calibration = calibration,
      point_counts = point_counts,
      replicates = replicates,
      mae = mae,
      groups = groups_used,
      parameters = list(
        ndim = ndim, k0 = k0, cooling_rate = cooling_rate,
        c_repulsion = c_repulsion, mapping_max_iter = mapping_max_iter,
        relative_epsilon = relative_epsilon,
        convergence_counter = convergence_counter
      ),
      settings = list(
        n_replicates = n_replicates,
        n_successful = n_successful,
        n_drawn = n_draw,
        n_pool = length(pool_idx),
        retained_fraction = retained_fraction,
        jackknife_factor = jackknife_factor,
        min_replicates = min_replicates,
        on_disconnect = on_disconnect,
        random_seed = random_seed,
        preserve_order = preserve_order
      )
    ),
    class = "topolow_bootstrap"
  )

  # ---------------------------------------------------------------------------
  # Optional output files
  # ---------------------------------------------------------------------------
  if (write_files) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }
    utils::write.csv(est_distances,
                     file.path(output_dir, "bootstrap_mean_distances.csv"))
    utils::write.csv(sd_distances,
                     file.path(output_dir, "bootstrap_sd_distances.csv"))
    utils::write.csv(cnt_mat,
                     file.path(output_dir, "bootstrap_counts.csv"))
    utils::write.csv(replicates,
                     file.path(output_dir, "bootstrap_replicates.csv"),
                     row.names = FALSE)
    if (!is.null(consensus_map)) {
      utils::write.csv(consensus_map$positions,
                       file.path(output_dir, "bootstrap_consensus_positions.csv"))
    }
    if (verbose) cat("Bootstrap output written to:", output_dir, "\n")
  }

  result
}


# =============================================================================
# Uncertainty regions
# =============================================================================

#' Project Coordinates Through a Fitted Reduction
#'
#' @param X Numeric matrix of coordinates in the embedding space.
#' @param ctr Numeric vector. Centre of the fitted reduction.
#' @param rot Numeric matrix. Rotation (loadings) of the fitted reduction.
#' @param sf Numeric. Scalar applied after rotation.
#' @param k Integer. Number of components to keep.
#' @return A numeric matrix with `k` columns.
#' @keywords internal
project_coords <- function(X, ctr, rot, sf, k) {
  sweep(X, 2, ctr) %*% rot[, seq_len(k), drop = FALSE] * sf
}


#' Uncertainty Region for One Point
#'
#' @param pts Numeric matrix of projected replicate positions, two columns.
#' @param conf_level Numeric. Proportion of the spread the region should cover.
#' @param method Character. `"ellipse"` or `"density"`.
#' @param n_points Integer. Vertices used to draw an ellipse.
#' @return A list of two-column matrices of polygon vertices, possibly empty.
#' @keywords internal
point_blob <- function(pts, conf_level, method, n_points) {
  pts <- pts[stats::complete.cases(pts), , drop = FALSE]
  if (nrow(pts) < 3) {
    return(list())
  }

  if (identical(method, "density")) {
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("Package 'MASS' is required for the density method. Please install it.")
    }
    kd <- MASS::kde2d(pts[, 1], pts[, 2], n = 50)
    z <- sort(as.vector(kd$z), decreasing = TRUE)
    cz <- cumsum(z) / sum(z)
    lev <- z[which(cz >= conf_level)[1]]
    lines <- grDevices::contourLines(kd$x, kd$y, kd$z, levels = lev)
    return(lapply(lines, function(l) cbind(l$x, l$y)))
  }

  ctr <- colMeans(pts)
  S <- stats::cov(pts)
  e <- eigen(S, symmetric = TRUE)
  rad <- sqrt(stats::qchisq(conf_level, df = 2))
  ax <- sqrt(pmax(e$values, 0)) * rad
  th <- seq(0, 2 * pi, length.out = n_points + 1)[-(n_points + 1)]
  poly <- cbind(cos(th), sin(th)) %*% diag(ax, 2, 2) %*% t(e$vectors)
  list(sweep(poly, 2, ctr, "+"))
}


#' Build Uncertainty Blobs from Aligned Replicate Clouds
#'
#' @param clouds List of aligned replicate coordinate matrices.
#' @param ctr Numeric vector. Centre of the fitted reduction.
#' @param rot Numeric matrix. Rotation of the fitted reduction.
#' @param sf Numeric. Scalar applied after rotation.
#' @param point_names Character vector of the points to summarise.
#' @param conf_level Numeric. Proportion of the spread covered.
#' @param method Character. `"ellipse"` or `"density"`.
#' @param n_points Integer. Vertices used to draw an ellipse.
#' @return A data frame with `name`, `blob_id`, `order`, `dim1` and `dim2`.
#' @keywords internal
blobs_from_clouds <- function(clouds, ctr, rot, sf, point_names,
                              conf_level, method, n_points) {
  projected <- lapply(clouds, function(C) {
    out <- matrix(NA_real_, nrow(C), 2, dimnames = list(rownames(C), NULL))
    ok <- stats::complete.cases(C)
    if (any(ok)) {
      out[ok, ] <- project_coords(C[ok, , drop = FALSE], ctr, rot, sf, 2)
    }
    out
  })

  rows <- lapply(point_names, function(nm) {
    pts <- do.call(rbind, lapply(projected, function(P) {
      if (nm %in% rownames(P)) P[nm, , drop = FALSE] else NULL
    }))
    if (is.null(pts)) return(NULL)
    polys <- point_blob(pts, conf_level, method, n_points)
    if (length(polys) == 0) return(NULL)
    do.call(rbind, lapply(seq_along(polys), function(k) {
      p <- polys[[k]]
      data.frame(name = nm,
                 blob_id = paste0(nm, "__", k),
                 order = seq_len(nrow(p)),
                 dim1 = p[, 1],
                 dim2 = p[, 2],
                 stringsAsFactors = FALSE)
    }))
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) {
    return(data.frame(name = character(0), blob_id = character(0),
                      order = integer(0), dim1 = numeric(0), dim2 = numeric(0),
                      stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}


#' Uncertainty Regions from a Bootstrap Embedding
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Turns the aligned per-replicate coordinates of a [bootstrap_embedding()]
#' result into per-point polygons that can be drawn on a map. This is the
#' topolow equivalent of a bootstrap blob: the replicate positions of a point
#' are projected into the plotting plane and summarised by a region covering
#' `conf_level` of their spread.
#'
#' @details
#' The regions cannot be derived from `sd_distances`. A distance matrix is
#' invariant to rotation, reflection and translation, so it carries no
#' directional information and can only yield an isotropic radius. Elongated
#' regions require the coordinates themselves, which is why
#' [bootstrap_embedding()] stores them rigidly aligned to the full-data map.
#'
#' Because topolow embeddings generally have more than two dimensions, the
#' aligned clouds are projected with the same principal-component transform that
#' [reduce_dimensions()] applies to the plotted coordinates. The resulting
#' regions therefore understate uncertainty in the directions that projection
#' discards.
#'
#' @param bootstrap_result A `topolow_bootstrap` object produced by
#'        [bootstrap_embedding()] with `keep_replicates = TRUE`.
#' @param conf_level Numeric. Proportion of the replicate spread each region
#'        should cover. Default: 0.68.
#' @param method Character. `"ellipse"` (default) fits a covariance ellipse and
#'        needs no extra package; `"density"` draws a kernel-density contour and
#'        requires the suggested package `MASS`.
#' @param dim_config A `dim_reduction_config` object controlling the projection.
#'        Only `method = "pca"` is supported.
#'        Default: [new_dim_reduction_config()].
#' @param n_points Integer. Number of vertices used to draw an ellipse.
#'        Default: 64.
#'
#' @return A `data.frame` with one row per polygon vertex and the columns
#'   `name`, `antigen`, `antiserum`, `blob_id`, `order`, `dim1` and `dim2`.
#'   `blob_id` identifies a single closed polygon and is what a plotting layer
#'   should group on.
#'
#' @examples
#' \donttest{
#' nms <- c("V/a1", "V/a2", "V/a3", "V/a4",
#'          "S/s1", "S/s2", "S/s3", "S/s4", "S/s5")
#' m <- as.matrix(dist(matrix(rnorm(27), nrow = 9)))
#' dimnames(m) <- list(nms, nms)
#'
#' res <- bootstrap_embedding(
#'   dissimilarity_matrix = m, ndim = 2,
#'   k0 = 1.0, cooling_rate = 0.01, c_repulsion = 0.01,
#'   mapping_max_iter = 30, n_replicates = 10, n_sera = 4, max_cores = 1
#' )
#' blobs <- bootstrap_blobs(res)
#' head(blobs)
#' }
#'
#' @seealso
#' [bootstrap_embedding()] for the replicates,
#' [plot_temporal_mapping()] and [plot_cluster_mapping()], which draw these
#' regions directly when given the bootstrap object.
#'
#' @importFrom stats cov qchisq complete.cases
#' @importFrom grDevices contourLines
#' @export
bootstrap_blobs <- function(bootstrap_result,
                            conf_level = 0.68,
                            method = c("ellipse", "density"),
                            dim_config = new_dim_reduction_config(),
                            n_points = 64) {
  method <- match.arg(method)

  if (!inherits(bootstrap_result, "topolow_bootstrap")) {
    stop("bootstrap_result must be a topolow_bootstrap object")
  }
  if (is.null(bootstrap_result$replicate_positions)) {
    stop("bootstrap_result must be created with keep_replicates = TRUE")
  }
  if (!identical(dim_config$method, "pca")) {
    stop("dim_config$method must be 'pca': other reductions cannot project new points")
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1 ||
      conf_level <= 0 || conf_level >= 1) {
    stop("conf_level must be a number between 0 and 1")
  }

  pos <- bootstrap_result$positions
  if (is.null(pos)) pos <- bootstrap_result$full_map$positions

  df <- as.data.frame(pos)
  names(df) <- paste0("V", seq_len(ncol(pos)))
  df$name <- rownames(pos)
  reduced <- reduce_dimensions(df, dim_config)

  out <- blobs_from_clouds(
    clouds = bootstrap_result$replicate_positions,
    ctr = attr(reduced, "center"),
    rot = attr(reduced, "rotation"),
    sf = attr(reduced, "scale_factor"),
    point_names = rownames(pos),
    conf_level = conf_level,
    method = method,
    n_points = n_points
  )

  out$antigen <- out$name %in% bootstrap_result$groups$antigens
  out$antiserum <- out$name %in% bootstrap_result$groups$sera
  out[, c("name", "antigen", "antiserum", "blob_id", "order", "dim1", "dim2")]
}


# =============================================================================
# S3 methods
# =============================================================================

#' Print method for topolow_bootstrap objects
#'
#' @param x A `topolow_bootstrap` object returned by [bootstrap_embedding()].
#' @param ... Additional arguments (not used).
#' @return The object, invisibly. Called for the summary it prints.
#' @examples
#' nms <- c("V/a1", "V/a2", "V/a3", "V/a4",
#'          "S/s1", "S/s2", "S/s3", "S/s4", "S/s5")
#' m <- as.matrix(dist(matrix(rnorm(18), nrow = 9)))
#' dimnames(m) <- list(nms, nms)
#' res <- bootstrap_embedding(m, ndim = 2, k0 = 1, cooling_rate = 0.01,
#'                            c_repulsion = 0.01, mapping_max_iter = 30,
#'                            n_replicates = 6, n_sera = 4,
#'                            consensus = FALSE, max_cores = 1)
#' print(res)
#' @export
print.topolow_bootstrap <- function(x, ...) {
  s <- x$settings
  cat("topolow bootstrap result:\n")
  cat(sprintf("Replicates: %d of %d succeeded\n", s$n_successful, s$n_replicates))
  cat(sprintf("Points: %d (%d fixed, %d drawable)\n",
              nrow(x$est_distances), length(x$groups$antigens),
              length(x$groups$sera)))
  cat(sprintf("Drawn per replicate: %d of %d (f = %.2f)\n",
              s$n_drawn, s$n_pool, s$retained_fraction))
  cat(sprintf("Dimensions: %d\n", x$parameters$ndim))
  cat(sprintf("Mean replicate MAE: %.4f\n",
              mean(x$replicates$mae, na.rm = TRUE)))
  cat(sprintf("MAE of the mean distance matrix: %.4f\n", x$mae))
  n_missing <- sum(x$n_obs < s$min_replicates)
  cat(sprintf("Cells below min_replicates (%d): %d of %d\n",
              s$min_replicates, n_missing, length(x$n_obs)))
  invisible(x)
}


#' Summary method for topolow_bootstrap objects
#'
#' @description
#' Prints the information shown by [print.topolow_bootstrap()] plus the
#' embedding parameters, the spread of the standard deviations per block, the
#' calibration against the full-data map, and the replicate convergence rate.
#'
#' @param object A `topolow_bootstrap` object returned by
#'        [bootstrap_embedding()].
#' @param ... Additional arguments (not used).
#' @return The object, invisibly. Called for the summary it prints.
#' @examples
#' nms <- c("V/a1", "V/a2", "V/a3", "V/a4",
#'          "S/s1", "S/s2", "S/s3", "S/s4", "S/s5")
#' m <- as.matrix(dist(matrix(rnorm(18), nrow = 9)))
#' dimnames(m) <- list(nms, nms)
#' res <- bootstrap_embedding(m, ndim = 2, k0 = 1, cooling_rate = 0.01,
#'                            c_repulsion = 0.01, mapping_max_iter = 30,
#'                            n_replicates = 6, n_sera = 4,
#'                            consensus = FALSE, max_cores = 1)
#' summary(res)
#' @importFrom stats quantile prcomp
#' @export
summary.topolow_bootstrap <- function(object, ...) {
  print(object)
  p <- object$parameters
  s <- object$settings

  cat("\nParameters:\n")
  cat(sprintf("  k0: %.4f\n", p$k0))
  cat(sprintf("  cooling_rate: %.4f\n", p$cooling_rate))
  cat(sprintf("  c_repulsion: %.4f\n", p$c_repulsion))
  cat(sprintf("  on_disconnect: %s\n", s$on_disconnect))
  if (!is.na(s$jackknife_factor)) {
    cat(sprintf("  jackknife factor for full-data SD: %.3f (not applied)\n",
                s$jackknife_factor))
  }

  ag <- object$groups$antigens
  sr <- object$groups$sera
  blocks <- list(
    "antigen-antigen" = list(r = ag, c = ag, tri = TRUE),
    "antigen-serum"   = list(r = ag, c = sr, tri = FALSE),
    "serum-serum"     = list(r = sr, c = sr, tri = TRUE)
  )
  cat("\nStandard deviations by block (median [IQR], mean coverage):\n")
  for (nm in names(blocks)) {
    b <- blocks[[nm]]
    if (length(b$r) < 2 || length(b$c) < 2) next
    sdb <- object$sd_distances[b$r, b$c, drop = FALSE]
    nob <- object$n_obs[b$r, b$c, drop = FALSE]
    if (b$tri) {
      sdb <- sdb[upper.tri(sdb)]
      nob <- nob[upper.tri(nob)]
    }
    sdb <- sdb[is.finite(sdb)]
    if (length(sdb) == 0) next
    q <- stats::quantile(sdb, c(0.25, 0.5, 0.75), names = FALSE)
    cat(sprintf("  %-16s %.4f [%.4f, %.4f]   n_obs mean %.1f\n",
                nm, q[2], q[1], q[3], mean(nob)))
  }

  if (!is.null(object$calibration)) {
    cat("\nCalibration against the full-data map:\n")
    print(object$calibration, row.names = FALSE)
  }

  conv <- mean(object$replicates$converged, na.rm = TRUE)
  cat(sprintf("\nReplicates reaching convergence: %.1f%%\n", 100 * conv))
  failed <- object$replicates$status != "ok"
  if (any(failed)) {
    cat("Failures by status:\n")
    print(table(object$replicates$status[failed]))
  }

  if (!is.null(object$positions) && ncol(object$positions) > 2) {
    pc <- stats::prcomp(object$positions, center = TRUE, scale. = FALSE)
    frac <- sum(pc$sdev[1:2]^2) / sum(pc$sdev^2)
    cat(sprintf("Variance retained by the two plotted components: %.1f%%\n",
                100 * frac))
  }

  invisible(object)
}


#' Prepare Uncertainty Blobs for a Mapping Plot
#'
#' @description
#' Internal helper shared by [plot_temporal_mapping()] and
#' [plot_cluster_mapping()]. Accepts either a `topolow_bootstrap` object or a
#' data frame already produced by [bootstrap_blobs()], and returns blob
#' polygons expressed in the reduced coordinates of the plot.
#'
#' @param uncertainty A `topolow_bootstrap` object or a blob data frame.
#' @param reduced_df The data frame returned by [reduce_dimensions()], carrying
#'   the `center`, `rotation` and `scale_factor` attributes.
#' @param conf_level Numeric. Proportion of the replicate spread covered.
#' @param n_points Integer. Vertices used to draw an ellipse.
#' @param points Character. Which points get a region: `"antigens"`,
#'   `"antisera"` or `"all"`.
#'
#' @return A data frame with `name`, `blob_id`, `order`, `dim1` and `dim2`.
#' @keywords internal
prepare_uncertainty_blobs <- function(uncertainty, reduced_df,
                                      conf_level = 0.68, n_points = 64,
                                      points = c("antigens", "antisera", "all")) {
  points <- match.arg(points)
  needed <- c("name", "blob_id", "order", "dim1", "dim2")

  # Which points should carry a region. Antisera are excluded by default:
  # drawing every point's region at once makes a dense map unreadable.
  keep_names <- reduced_df$name
  if (identical(points, "antigens")) {
    keep_names <- reduced_df$name[as.logical(reduced_df$antigen)]
  } else if (identical(points, "antisera")) {
    keep_names <- reduced_df$name[as.logical(reduced_df$antiserum)]
  }
  keep_names <- keep_names[!is.na(keep_names)]

  if (is.data.frame(uncertainty)) {
    if (!all(needed %in% names(uncertainty))) {
      stop("uncertainty data frame must have columns: ",
           paste(needed, collapse = ", "))
    }
    return(uncertainty[uncertainty$name %in% keep_names, , drop = FALSE])
  }

  if (!inherits(uncertainty, "topolow_bootstrap")) {
    stop("uncertainty must be a topolow_bootstrap object or a blob data frame")
  }
  if (is.null(uncertainty$replicate_positions)) {
    stop("uncertainty must come from bootstrap_embedding() with keep_replicates = TRUE")
  }

  ctr <- attr(reduced_df, "center")
  rot <- attr(reduced_df, "rotation")
  sf <- attr(reduced_df, "scale_factor")
  if (is.null(ctr) || is.null(rot)) {
    stop("uncertainty regions require dim_config$method to be 'pca'")
  }
  if (ncol(uncertainty$replicate_positions[[1]]) != length(ctr)) {
    stop("uncertainty was computed in a different number of dimensions than df_coords")
  }
  if (!"name" %in% names(reduced_df)) {
    stop("df_coords must have a name column to match uncertainty regions")
  }

  shared <- intersect(keep_names,
                      rownames(uncertainty$replicate_positions[[1]]))
  if (length(shared) == 0) {
    if (length(intersect(reduced_df$name,
                         rownames(uncertainty$replicate_positions[[1]]))) == 0) {
      stop("no point names are shared between df_coords and uncertainty")
    }
    # Names match, but nothing of the requested type is left to draw.
    return(data.frame(name = character(0), blob_id = character(0),
                      order = integer(0), dim1 = numeric(0), dim2 = numeric(0),
                      stringsAsFactors = FALSE))
  }

  blobs_from_clouds(uncertainty$replicate_positions, ctr, rot, sf,
                    shared, conf_level, "ellipse", n_points)
}
