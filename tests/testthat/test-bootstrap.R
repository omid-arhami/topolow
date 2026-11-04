# tests/testthat/test-bootstrap.R

# Copyright (c) 2026 Omid Arhami omid.arhami@uga.edu

# Helper: a small antigen/serum dissimilarity matrix in the V/ S/ convention
make_vs_matrix <- function(n_v = 4, n_s = 3, ndim = 2, seed = 1) {
  set.seed(seed)
  m <- as.matrix(dist(matrix(rnorm((n_v + n_s) * ndim), ncol = ndim)))
  nm <- c(paste0("V/ag", seq_len(n_v)), paste0("S/sr", seq_len(n_s)))
  dimnames(m) <- list(nm, nm)
  m
}

run_small_bootstrap <- function(m, ...) {
  args <- list(
    dissimilarity_matrix = m, ndim = 2, k0 = 1, cooling_rate = 0.01,
    c_repulsion = 0.01, mapping_max_iter = 30, n_replicates = 5,
    n_sera = 2, max_cores = 1, random_seed = 42
  )
  args <- utils::modifyList(args, list(...))
  suppressWarnings(do.call(bootstrap_embedding, args))
}

# =============================================================================
# TEST GROUP 1: Procrustes alignment
# =============================================================================

test_that("procrustes_align recovers a known rotation", {
  set.seed(5)
  Y <- matrix(rnorm(20), ncol = 2)
  rownames(Y) <- paste0("p", 1:10)
  theta <- pi / 3
  R <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2, 2)
  X <- sweep(Y %*% R, 2, c(3, -2), "+")
  rownames(X) <- rownames(Y)

  aligned <- procrustes_align(X, Y)
  expect_equal(unname(aligned), unname(Y), tolerance = 1e-10)
})

test_that("procrustes_align preserves all pairwise distances", {
  set.seed(6)
  Y <- matrix(rnorm(24), ncol = 3)
  X <- matrix(rnorm(24), ncol = 3)
  rownames(X) <- rownames(Y) <- paste0("p", 1:8)

  aligned <- procrustes_align(X, Y)
  expect_equal(as.matrix(dist(aligned)), as.matrix(dist(X)), tolerance = 1e-10)
})

test_that("procrustes_align needs shared points", {
  X <- matrix(rnorm(4), ncol = 2, dimnames = list(c("a", "b"), NULL))
  Y <- matrix(rnorm(4), ncol = 2, dimnames = list(c("c", "d"), NULL))
  expect_error(procrustes_align(X, Y), "share at least two named points")
})

# =============================================================================
# TEST GROUP 2: Structure of the returned object
# =============================================================================

test_that("bootstrap_embedding returns the documented structure", {
  m <- make_vs_matrix()
  res <- run_small_bootstrap(m)

  expect_s3_class(res, "topolow_bootstrap")
  for (nm in c("est_distances", "sd_distances", "n_obs", "positions",
               "consensus_map", "replicate_positions", "full_map",
               "calibration", "point_counts", "replicates", "mae",
               "groups", "parameters", "settings")) {
    expect_true(nm %in% names(res), info = nm)
  }
  expect_equal(dim(res$est_distances), c(7L, 7L))
  expect_equal(rownames(res$est_distances), rownames(m))
  expect_equal(colnames(res$est_distances), colnames(m))
  expect_s3_class(res$replicates, "data.frame")
  expect_equal(nrow(res$replicates), 5L)
})

test_that("bootstrap aggregates are internally consistent", {
  m <- make_vs_matrix()
  res <- run_small_bootstrap(m)

  expect_equal(res$est_distances, t(res$est_distances))
  expect_equal(res$sd_distances, t(res$sd_distances))
  expect_true(all(diag(res$est_distances) == 0))
  expect_true(is.integer(res$n_obs))
  expect_true(all(res$n_obs <= res$settings$n_successful))
  expect_true(all(res$sd_distances >= 0, na.rm = TRUE))
})

test_that("min_replicates gates the mean and the standard deviation alike", {
  m <- make_vs_matrix()
  res <- run_small_bootstrap(m)

  below <- res$n_obs < res$settings$min_replicates
  expect_true(all(is.na(res$est_distances[below])))
  expect_true(all(is.na(res$sd_distances[below])))
  expect_true(all(!is.na(res$est_distances[!below])))
  expect_true(all(!is.na(res$sd_distances[!below])))
  # the two matrices are missing in exactly the same places
  expect_identical(is.na(res$est_distances), is.na(res$sd_distances))
})

test_that("min_replicates default is 2", {
  m <- make_vs_matrix()
  res <- run_small_bootstrap(m)
  expect_equal(res$settings$min_replicates, 2L)
})

test_that("preserve_order defaults to TRUE", {
  m <- make_vs_matrix()
  res <- run_small_bootstrap(m)
  expect_true(res$settings$preserve_order)
})

# =============================================================================
# TEST GROUP 3: The sampling design
# =============================================================================

test_that("every always_include point appears in every replicate", {
  m <- make_vs_matrix(n_v = 5, n_s = 4)
  res <- run_small_bootstrap(m, n_sera = 2)

  antigens <- res$groups$antigens
  expect_equal(length(antigens), 5L)
  expect_true(all(res$point_counts[antigens] == res$settings$n_successful))
  expect_true(all(res$point_counts[res$groups$sera] <= res$settings$n_successful))
})

test_that("drawing every serum gives full coverage", {
  m <- make_vs_matrix(n_v = 4, n_s = 3)
  res <- suppressWarnings(run_small_bootstrap(m, frac_sera = 1, n_sera = NULL))

  expect_true(all(res$n_obs == res$settings$n_successful))
  expect_false(any(is.na(res$est_distances)))
  expect_true(is.na(res$settings$jackknife_factor))
})

test_that("the jackknife factor is sqrt(f / (1 - f))", {
  m <- make_vs_matrix(n_v = 4, n_s = 4)
  res <- run_small_bootstrap(m, n_sera = 2)
  expect_equal(res$settings$retained_fraction, 0.5)
  expect_equal(res$settings$jackknife_factor, 1, tolerance = 1e-12)
})

test_that("point draws are reproducible under random_seed", {
  m <- make_vs_matrix(n_v = 5, n_s = 4)
  r1 <- run_small_bootstrap(m, random_seed = 123)
  r2 <- run_small_bootstrap(m, random_seed = 123)

  # The draws are seeded; the embeddings are not, so only the design is checked.
  expect_identical(r1$point_counts, r2$point_counts)
  expect_identical(r1$n_obs, r2$n_obs)
})

test_that("generic always_include and sample_pool override the prefixes", {
  m <- make_vs_matrix(n_v = 4, n_s = 3)
  nm <- rownames(m)
  res <- run_small_bootstrap(m,
                             always_include = nm[1:3],
                             sample_pool = nm[4:7],
                             n_sera = 2)
  expect_equal(res$groups$antigens, nm[1:3])
  expect_true(all(res$point_counts[nm[1:3]] == res$settings$n_successful))
})

# =============================================================================
# TEST GROUP 4: Consensus map and replicate clouds
# =============================================================================

test_that("the consensus map covers every point", {
  m <- make_vs_matrix()
  res <- run_small_bootstrap(m, consensus = TRUE)

  expect_s3_class(res$consensus_map, "topolow")
  expect_equal(nrow(res$positions), 7L)
  expect_equal(ncol(res$positions), 2L)
  expect_setequal(rownames(res$positions), rownames(m))
  expect_identical(res$positions, res$consensus_map$positions)
})

test_that("consensus = FALSE leaves the positions empty", {
  m <- make_vs_matrix()
  res <- run_small_bootstrap(m, consensus = FALSE)
  expect_null(res$positions)
  expect_null(res$consensus_map)
})

test_that("replicate clouds are full size with NA rows for absent points", {
  m <- make_vs_matrix(n_v = 4, n_s = 3)
  res <- run_small_bootstrap(m, n_sera = 2)

  expect_type(res$replicate_positions, "list")
  expect_equal(length(res$replicate_positions), res$settings$n_successful)
  cloud <- res$replicate_positions[[1]]
  expect_equal(dim(cloud), c(7L, 2L))
  expect_equal(rownames(cloud), rownames(m))
  # all four antigens present, one of the three sera missing
  expect_false(any(is.na(cloud[res$groups$antigens, ])))
  expect_equal(sum(stats::complete.cases(cloud)), 6L)
})

test_that("keep_replicates = FALSE drops the clouds", {
  m <- make_vs_matrix()
  res <- run_small_bootstrap(m, keep_replicates = FALSE)
  expect_null(res$replicate_positions)
  expect_error(bootstrap_blobs(res), "keep_replicates = TRUE")
})

# =============================================================================
# TEST GROUP 5: Blobs
# =============================================================================

test_that("bootstrap_blobs returns one closed polygon per point", {
  m <- make_vs_matrix(n_v = 4, n_s = 3)
  res <- run_small_bootstrap(m, n_replicates = 10, mapping_max_iter = 40)
  blobs <- bootstrap_blobs(res, n_points = 16)

  expect_s3_class(blobs, "data.frame")
  expect_named(blobs, c("name", "antigen", "antiserum", "blob_id",
                        "order", "dim1", "dim2"))
  expect_true(all(table(blobs$blob_id) == 16))
  expect_true(all(blobs$name %in% rownames(m)))
  expect_true(all(is.finite(blobs$dim1)))
  expect_true(all(is.finite(blobs$dim2)))
  expect_true(any(blobs$antigen))
  expect_true(any(blobs$antiserum))
})

test_that("a larger conf_level gives a larger region", {
  m <- make_vs_matrix(n_v = 4, n_s = 3)
  res <- run_small_bootstrap(m, n_replicates = 12, mapping_max_iter = 40)
  small <- bootstrap_blobs(res, conf_level = 0.5)
  large <- bootstrap_blobs(res, conf_level = 0.95)

  span <- function(b) diff(range(b$dim1)) + diff(range(b$dim2))
  expect_gt(span(large), span(small))
})

test_that("bootstrap_blobs validates its arguments", {
  m <- make_vs_matrix()
  res <- run_small_bootstrap(m)
  expect_error(bootstrap_blobs(list()), "topolow_bootstrap object")
  expect_error(bootstrap_blobs(res, conf_level = 1.5), "between 0 and 1")
  expect_error(
    bootstrap_blobs(res, dim_config = new_dim_reduction_config(method = "umap")),
    "pca"
  )
})

# =============================================================================
# TEST GROUP 6: Validation and failure handling
# =============================================================================

test_that("bootstrap_embedding rejects invalid arguments", {
  m <- make_vs_matrix()

  expect_error(
    bootstrap_embedding("not a matrix", ndim = 2, k0 = 1,
                        cooling_rate = 0.01, c_repulsion = 0.01),
    "dissimilarity_matrix must be a matrix"
  )
  expect_error(
    bootstrap_embedding(m, ndim = 2, k0 = 1, cooling_rate = 0.01,
                        c_repulsion = 0.01, n_replicates = 1),
    "n_replicates must be an integer greater than or equal to 2"
  )
  expect_error(
    bootstrap_embedding(m, ndim = 2, k0 = 1, cooling_rate = 0.01,
                        c_repulsion = 0.01, min_replicates = 1),
    "min_replicates must be an integer greater than or equal to 2"
  )
  expect_error(
    bootstrap_embedding(m, ndim = 2, k0 = 1, cooling_rate = 0.01,
                        c_repulsion = 0.01, frac_sera = 0),
    "frac_sera must be a number greater than 0"
  )
  expect_error(
    bootstrap_embedding(m, ndim = 2, k0 = -1, cooling_rate = 0.01,
                        c_repulsion = 0.01),
    "k0 must be a positive number"
  )

  dup <- m
  rownames(dup)[2] <- rownames(dup)[1]
  colnames(dup) <- rownames(dup)
  expect_error(
    bootstrap_embedding(dup, ndim = 2, k0 = 1, cooling_rate = 0.01,
                        c_repulsion = 0.01),
    "unique row names"
  )
})

test_that("write_files without output_dir is an error, and nothing is written otherwise", {
  m <- make_vs_matrix()
  expect_error(
    bootstrap_embedding(m, ndim = 2, k0 = 1, cooling_rate = 0.01,
                        c_repulsion = 0.01, n_replicates = 2,
                        write_files = TRUE),
    "'output_dir' must be provided"
  )

  out_dir <- tempfile()
  dir.create(out_dir, recursive = TRUE)
  run_small_bootstrap(m)
  expect_equal(length(list.files(out_dir)), 0)

  run_small_bootstrap(m, write_files = TRUE, output_dir = out_dir)
  for (f in c("bootstrap_mean_distances.csv", "bootstrap_sd_distances.csv",
              "bootstrap_counts.csv", "bootstrap_replicates.csv",
              "bootstrap_consensus_positions.csv")) {
    expect_true(file.exists(file.path(out_dir, f)), info = f)
  }
  unlink(out_dir, recursive = TRUE)
})

test_that("an antigen that can never connect fails every replicate", {
  m <- make_vs_matrix(n_v = 4, n_s = 3)
  # Isolate one antigen: no observed dissimilarity to anything
  m["V/ag1", ] <- NA
  m[, "V/ag1"] <- NA
  m["V/ag1", "V/ag1"] <- 0

  expect_error(
    suppressWarnings(bootstrap_embedding(
      m, ndim = 2, k0 = 1, cooling_rate = 0.01, c_repulsion = 0.01,
      mapping_max_iter = 20, n_replicates = 3, n_sera = 2,
      on_disconnect = "error", max_cores = 1, consensus = FALSE
    )),
    "replicates failed"
  )
})

test_that("on_disconnect = prune drops the isolated point and keeps going", {
  m <- make_vs_matrix(n_v = 4, n_s = 3)
  m["V/ag1", ] <- NA
  m[, "V/ag1"] <- NA
  m["V/ag1", "V/ag1"] <- 0

  res <- suppressWarnings(bootstrap_embedding(
    m, ndim = 2, k0 = 1, cooling_rate = 0.01, c_repulsion = 0.01,
    mapping_max_iter = 20, n_replicates = 4, n_sera = 2,
    on_disconnect = "prune", max_cores = 1, consensus = FALSE
  ))

  expect_equal(res$settings$n_successful, 4L)
  expect_equal(unname(res$point_counts["V/ag1"]), 0L)
  expect_true(all(res$replicates$n_pruned >= 1))
})

test_that("an unmatched point count warns rather than silently dropping", {
  m <- make_vs_matrix(n_v = 4, n_s = 3)
  m["V/ag1", ] <- NA
  m[, "V/ag1"] <- NA
  m["V/ag1", "V/ag1"] <- 0

  # Several warnings fire together here, so collect them all and check one.
  seen <- character()
  withCallingHandlers(
    bootstrap_embedding(
      m, ndim = 2, k0 = 1, cooling_rate = 0.01, c_repulsion = 0.01,
      mapping_max_iter = 20, n_replicates = 4, n_sera = 2,
      on_disconnect = "prune", max_cores = 1, consensus = FALSE
    ),
    warning = function(cond) {
      seen <<- c(seen, conditionMessage(cond))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("appeared in no replicate", seen)))
})

# =============================================================================
# TEST GROUP 7: Parallel execution
# =============================================================================

test_that("bootstrap_embedding runs with more than one core", {
  skip_on_cran()
  skip_if_not_installed("parallel")

  m <- make_vs_matrix(n_v = 5, n_s = 4)
  res <- run_small_bootstrap(m, n_replicates = 4, max_cores = 2)

  expect_s3_class(res, "topolow_bootstrap")
  expect_equal(nrow(res$replicates), 4L)
  expect_gt(res$settings$n_successful, 0L)
})

test_that("every replicate draws the same number of points", {
  m <- make_vs_matrix(n_v = 5, n_s = 4)
  res <- run_small_bootstrap(m, n_replicates = 8, n_sera = 3)

  expect_length(unique(res$replicates$n_sera_drawn), 1)
  expect_equal(unique(res$replicates$n_sera_drawn), 3)
  expect_length(unique(res$replicates$n_points), 1)
  expect_equal(unique(res$replicates$n_points), 8L)
  expect_true(all(res$replicates$n_pruned == 0))
})

test_that("frac_sera resolves to a fixed count that the object reports", {
  m <- make_vs_matrix(n_v = 5, n_s = 10)
  res <- run_small_bootstrap(m, n_replicates = 4, frac_sera = 0.8,
                             n_sera = NULL)

  expect_equal(res$settings$n_drawn, 8L)
  expect_equal(res$settings$n_pool, 10L)
  expect_equal(res$settings$retained_fraction, 0.8)
  expect_equal(unique(res$replicates$n_points), 13L)
})
