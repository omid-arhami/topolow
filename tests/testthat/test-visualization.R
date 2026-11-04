# tests/testthat/test-visualization.R

# Copyright (c) 2025 Omid Arhami omid.arhami@uga.edu

test_that("scatterplot_fitted_vs_true creates valid plots", {
  true_dist <- matrix(runif(100, 1, 10), 10, 10)
  pred_dist <- true_dist + rnorm(100)

  plots <- scatterplot_fitted_vs_true(
    dissimilarity_matrix = true_dist,
    p_dissimilarity_mat = pred_dist,
    save_plot = FALSE
  )

  expect_true(ggplot2::is_ggplot(plots$scatter_plot))
  expect_true(ggplot2::is_ggplot(plots$residuals_plot))
})

test_that("plot_network_structure creates a valid plot", {
  adj_mat <- matrix(runif(25), 5, 5)
  diag(adj_mat) <- 0
  net_analysis <- analyze_network_structure(adj_mat)

  p <- plot_network_structure(net_analysis)
  expect_true(ggplot2::is_ggplot(p))
})

test_that("plot_mcmc_diagnostics works with temporary files", {
  # Create dummy chain files in a temporary directory
  temp_dir <- tempdir()
  chain_files <- character(2)
  par_names <- c("log_N", "log_k0", "log_cooling_rate", "log_c_repulsion")
  sample_data <- data.frame(
    log_N = rnorm(100), log_k0 = rnorm(100),
    log_cooling_rate = rnorm(100), log_c_repulsion = rnorm(100)
  )
  for (i in 1:2) {
    chain_files[i] <- file.path(temp_dir, paste0("chain", i, ".csv"))
    write.csv(sample_data, chain_files[i], row.names = FALSE)
  }

  p <- plot_mcmc_diagnostics(chain_files, mutual_size = 50, save_plot = FALSE)
  # grid.arrange returns a gtable object
  expect_s3_class(p, "gtable")

  unlink(chain_files)
})

test_that("plot_temporal_mapping creates valid plot", {
  # Create test data
  test_df <- data.frame(
    V1 = rnorm(10),
    V2 = rnorm(10),
    antigen = rep(c(TRUE, FALSE), 5),
    antiserum = rep(c(FALSE, TRUE), 5),
    name = paste0(ifelse(rep(c(TRUE, FALSE), 5), "V/", "S/"), "strain", 1:10),
    year = 2000:2009
  )
  
  plot <- plot_temporal_mapping(test_df, ndim=2)
  expect_true(ggplot2::is_ggplot(plot))
  
  # Test configuration objects
  aesthetic_config <- new_aesthetic_config()
  layout_config <- new_layout_config()
  plot <- plot_temporal_mapping(
    test_df, ndim=2,
    aesthetic_config = aesthetic_config,
    layout_config = layout_config
  )
  expect_true(ggplot2::is_ggplot(plot))
})

test_that("plot_cluster_mapping creates valid plot", {
  test_df <- data.frame(
    V1 = rnorm(10),
    V2 = rnorm(10),
    antigen = rep(c(TRUE, FALSE), 5),
    antiserum = rep(c(FALSE, TRUE), 5),
    cluster = rep(c("A", "B"), 5),
    name = paste0(ifelse(rep(c(TRUE, FALSE), 5), "V/", "S/"), "strain", 1:10)
  )
  
  plot <- plot_cluster_mapping(test_df, ndim=2, )
  expect_true(ggplot2::is_ggplot(plot))
})

test_that("plot configuration objects work correctly", {
  config <- new_aesthetic_config(
    point_size = 2,
    point_alpha = 0.5
  )
  expect_s3_class(config, "aesthetic_config")
  
  config <- new_layout_config(
    width = 8,
    height = 6
  )
  expect_s3_class(config, "layout_config")
})

# =============================================================================
# Uncertainty overlays and the fitted reduction transform
# =============================================================================

make_boot_fixture <- function(n_v = 5, n_s = 4, ndim = 3, seed = 21) {
  set.seed(seed)
  m <- as.matrix(dist(matrix(rnorm((n_v + n_s) * ndim), ncol = ndim)))
  nm <- c(paste0("V/ag", seq_len(n_v)), paste0("S/sr", seq_len(n_s)))
  dimnames(m) <- list(nm, nm)

  boot <- suppressWarnings(bootstrap_embedding(
    m, ndim = ndim, k0 = 1, cooling_rate = 0.01, c_repulsion = 0.01,
    mapping_max_iter = 40, n_replicates = 10, frac_sera = 0.75,
    max_cores = 1, random_seed = 3
  ))

  df <- as.data.frame(boot$positions)
  names(df) <- paste0("V", seq_len(ndim))
  df$name <- rownames(boot$positions)
  df$antigen <- startsWith(df$name, "V/")
  df$antiserum <- startsWith(df$name, "S/")
  df$year <- seq(1970, by = 2, length.out = nrow(df))
  df$cluster <- factor(rep(c("c1", "c2"), length.out = nrow(df)))

  list(boot = boot, df = df, ndim = ndim)
}

test_that("reduce_dimensions exposes a transform that reproduces its output", {
  skip_if_not_installed("ggplot2")
  set.seed(31)
  X <- matrix(rnorm(40), ncol = 4)
  rownames(X) <- paste0("p", 1:10)
  df <- as.data.frame(X)
  names(df) <- paste0("V", 1:4)
  df$name <- rownames(X)

  red <- reduce_dimensions(df, new_dim_reduction_config())

  ctr <- attr(red, "center")
  rot <- attr(red, "rotation")
  sf <- attr(red, "scale_factor")
  expect_false(is.null(ctr))
  expect_false(is.null(rot))
  expect_true(is.numeric(sf))

  proj <- sweep(X, 2, ctr) %*% rot[, 1:2] * sf
  expect_equal(unname(proj), unname(as.matrix(red[, c("dim1", "dim2")])),
               tolerance = 1e-10)
})

test_that("plot_temporal_mapping draws an uncertainty layer beneath the points", {
  skip_if_not_installed("ggplot2")
  fx <- make_boot_fixture()

  p0 <- plot_temporal_mapping(fx$df, ndim = fx$ndim)
  p1 <- plot_temporal_mapping(fx$df, ndim = fx$ndim, uncertainty = fx$boot)

  expect_s3_class(p1, "ggplot")
  expect_equal(length(p1$layers), length(p0$layers) + 1)
  expect_s3_class(p1$layers[[1]]$geom, "GeomPolygon")
})

test_that("plot_cluster_mapping draws an uncertainty layer beneath the points", {
  skip_if_not_installed("ggplot2")
  fx <- make_boot_fixture()

  p0 <- plot_cluster_mapping(fx$df, ndim = fx$ndim)
  p1 <- plot_cluster_mapping(fx$df, ndim = fx$ndim, uncertainty = fx$boot)

  expect_s3_class(p1, "ggplot")
  expect_equal(length(p1$layers), length(p0$layers) + 1)
  expect_s3_class(p1$layers[[1]]$geom, "GeomPolygon")
})

test_that("a precomputed blob data frame can be passed as uncertainty", {
  skip_if_not_installed("ggplot2")
  fx <- make_boot_fixture()
  blobs <- bootstrap_blobs(fx$boot)

  p <- plot_temporal_mapping(fx$df, ndim = fx$ndim, uncertainty = blobs)
  expect_s3_class(p$layers[[1]]$geom, "GeomPolygon")
})

test_that("uncertainty overlays reject mismatched input", {
  skip_if_not_installed("ggplot2")
  fx <- make_boot_fixture()

  expect_error(
    plot_temporal_mapping(fx$df, ndim = fx$ndim, uncertainty = list(a = 1)),
    "topolow_bootstrap object or a blob data frame"
  )
  expect_error(
    plot_temporal_mapping(fx$df, ndim = fx$ndim,
                          uncertainty = data.frame(x = 1)),
    "uncertainty data frame must have columns"
  )

  no_clouds <- fx$boot
  no_clouds$replicate_positions <- NULL
  expect_error(
    plot_temporal_mapping(fx$df, ndim = fx$ndim, uncertainty = no_clouds),
    "keep_replicates = TRUE"
  )
})

test_that("a plot with uncertainty renders to a file", {
  skip_if_not_installed("ggplot2")
  fx <- make_boot_fixture()
  p <- plot_temporal_mapping(fx$df, ndim = fx$ndim, uncertainty = fx$boot)

  f <- tempfile(fileext = ".png")
  ggplot2::ggsave(f, p, width = 6, height = 5, dpi = 72)
  expect_true(file.exists(f))
  expect_gt(file.info(f)$size, 1000)
  unlink(f)
})

test_that("uncertainty regions cover antigens only by default", {
  skip_if_not_installed("ggplot2")
  fx <- make_boot_fixture()

  p <- plot_temporal_mapping(fx$df, ndim = fx$ndim, uncertainty = fx$boot)
  drawn <- unique(p$layers[[1]]$data$name)

  expect_setequal(drawn, fx$boot$groups$antigens)
  expect_false(any(startsWith(drawn, "S/")))
})

test_that("uncertainty_points selects which points get a region", {
  skip_if_not_installed("ggplot2")
  fx <- make_boot_fixture()

  n_drawn <- function(p) length(unique(p$layers[[1]]$data$name))

  p_ag <- plot_temporal_mapping(fx$df, ndim = fx$ndim, uncertainty = fx$boot,
                                uncertainty_points = "antigens")
  p_sr <- plot_temporal_mapping(fx$df, ndim = fx$ndim, uncertainty = fx$boot,
                                uncertainty_points = "antisera")
  p_all <- plot_temporal_mapping(fx$df, ndim = fx$ndim, uncertainty = fx$boot,
                                 uncertainty_points = "all")

  expect_equal(n_drawn(p_ag), length(fx$boot$groups$antigens))
  expect_equal(n_drawn(p_sr), length(fx$boot$groups$sera))
  expect_equal(n_drawn(p_all), n_drawn(p_ag) + n_drawn(p_sr))
  expect_true(all(startsWith(unique(p_sr$layers[[1]]$data$name), "S/")))
})

test_that("plot_cluster_mapping honours uncertainty_points too", {
  skip_if_not_installed("ggplot2")
  fx <- make_boot_fixture()

  p_ag <- plot_cluster_mapping(fx$df, ndim = fx$ndim, uncertainty = fx$boot)
  p_all <- plot_cluster_mapping(fx$df, ndim = fx$ndim, uncertainty = fx$boot,
                                uncertainty_points = "all")

  expect_equal(length(unique(p_ag$layers[[1]]$data$name)),
               length(fx$boot$groups$antigens))
  expect_gt(length(unique(p_all$layers[[1]]$data$name)),
            length(unique(p_ag$layers[[1]]$data$name)))
})

test_that("a precomputed blob data frame is filtered the same way", {
  skip_if_not_installed("ggplot2")
  fx <- make_boot_fixture()
  blobs <- bootstrap_blobs(fx$boot)

  p <- plot_temporal_mapping(fx$df, ndim = fx$ndim, uncertainty = blobs)
  expect_setequal(unique(p$layers[[1]]$data$name), fx$boot$groups$antigens)
})

test_that("uncertainty_points rejects an unknown choice", {
  skip_if_not_installed("ggplot2")
  fx <- make_boot_fixture()
  expect_error(
    plot_temporal_mapping(fx$df, ndim = fx$ndim, uncertainty = fx$boot,
                          uncertainty_points = "viruses"),
    "should be one of"
  )
})
