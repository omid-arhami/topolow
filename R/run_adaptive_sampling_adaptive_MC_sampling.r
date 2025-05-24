run_adaptive_sampling <- function(initial_samples_file,
                                  scenario_name,
                                  distance_matrix,
                                  num_parallel_jobs = 5,
                                  max_cores = NULL,
                                  num_samples = 10,
                                  mapping_max_iter = 1000, 
                                  relative_epsilon = 1e-4,
                                  folds = 20, 
                                  time = "8:00:00",
                                  memory = "10G",
                                  output_dir = NULL,
                                  use_slurm = FALSE,
                                  cider = FALSE,
                                  verbose = FALSE) {
  # Parameter names and validation
  par_names <- c("log_N", "log_k0", "log_cooling_rate", "log_c_repulsion")
  if (!is.numeric(num_samples) || num_samples < 1 || num_samples != round(num_samples)) {
    stop("num_samples must be a positive integer")
  }
  iterations <- ceiling(num_samples / num_parallel_jobs)
  integer_params <- list(
    mapping_max_iter = mapping_max_iter,
    folds = folds,
    num_parallel_jobs = num_parallel_jobs,
    iterations = iterations
  )
  for (p in names(integer_params)) {
    v <- integer_params[[p]]
    if (!is.numeric(v) || v < 1 || v != round(v)) stop(sprintf("%s must be a positive integer", p))
  }
  if (!is.numeric(relative_epsilon) || relative_epsilon <= 0) stop("relative_epsilon must be positive")
  if (!is.character(scenario_name) || length(scenario_name) != 1) stop("scenario_name must be a single string")
  if (use_slurm && !has_slurm()) stop("SLURM requested but not available")

  # Determine cores
  available_cores <- parallel::detectCores()
  if (is.null(max_cores)) {
    max_cores <- max(1, available_cores - 1)
  } else {
    if (!is.numeric(max_cores) || max_cores < 1 || max_cores != round(max_cores)) {
      stop("max_cores must be a positive integer")
    }
    max_cores <- min(max_cores, available_cores)
  }
  if (verbose) cat(sprintf("Using %d jobs with up to %d cores\n", num_parallel_jobs, max_cores))

  # Setup directories
  if (is.null(output_dir)) output_dir <- getwd()
  adaptive_dir <- file.path(output_dir, "adaptive_sampling_jobs")
  param_dir    <- file.path(output_dir, "model_parameters")
  for (d in c(adaptive_dir, param_dir)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

  # Prepare initial and results files
  results_file <- file.path(param_dir, paste0(scenario_name, "_model_parameters.csv"))
  if (!file.exists(initial_samples_file)) stop("initial_samples_file not found: ", initial_samples_file)
  init <- read.csv(initial_samples_file, stringsAsFactors = FALSE)
  req <- c(par_names, "NLL")
  if (!all(req %in% names(init))) stop("Missing columns in initial samples: ", paste(setdiff(req, names(init)), collapse = ", "))
  init <- init[complete.cases(init[, req]) & apply(init[, req], 1, function(x) all(is.finite(x))), ]
  if (nrow(init) == 0) stop("No valid initial samples after filtering")
  if (!file.exists(results_file)) file.copy(initial_samples_file, results_file)

  # Per-job temp helper
  make_temp <- function(i) file.path(adaptive_dir, sprintf("job_%02d_%s.csv", i, scenario_name))

    if (use_slurm) {
    # 1) dump the distance matrix
    dist_file <- file.path(adaptive_dir, paste0(scenario_name, "_distance_matrix.rds"))
    saveRDS(distance_matrix, dist_file)

    # 2) schedule each sampler to its own temp CSV
    job_ids <- character(num_parallel_jobs)
    for (i in seq_len(num_parallel_jobs)) {
      temp_f <- make_temp(i)
      file.copy(initial_samples_file, temp_f, overwrite=TRUE)
      args <- c(
        temp_f,
        dist_file,
        as.character(mapping_max_iter),
        as.character(relative_epsilon),
        "1",               # inner batch_size
        scenario_name,
        as.character(i), 
        as.character(iterations),
        output_dir,
        as.character(folds)
      )
      sbatch <- create_slurm_script(
        job_name    = paste0("job", i, "_", scenario_name),
        script_path = system.file("scripts/run_adaptive_sampling_jobs.R", package="topolow"),
        args        = args,
        num_cores   = 1,
        time        = time,
        memory      = memory,
        output_file= file.path(adaptive_dir, paste0(i, ".out")),
        error_file = file.path(adaptive_dir, paste0(i, ".err"))
      )
      job_ids[i] <- submit_job(sbatch, use_slurm=TRUE, cider=cider)
      if (verbose) cat("Submitted sampler job", job_ids[i], "→", temp_f, "\n")
    }

    # 3) build the gather script (with an SBATCH header!)
    dep      <- paste(job_ids, collapse=":")
    results_file <- file.path(param_dir, paste0(scenario_name, "_model_parameters.csv"))
    gather_R <- file.path(adaptive_dir, paste0("gather_", scenario_name, ".R"))
    gather_code <- c(
      "#!/usr/bin/env Rscript",
      paste0("#SBATCH --dependency=afterok:", dep),
      "",
      "args <- commandArgs(trailingOnly=TRUE)",
      "# arg1 = results_file, arg2 = adaptive_dir, arg3 = scenario_name, arg4 = output_dir",
      "results_file <- args[1]",
      "adaptive_dir <- args[2]",
      "scenario_name <- args[3]",
      "output_dir <- args[4]",
      "",
      "init <- read.csv(results_file, stringsAsFactors=FALSE)",
      "n0   <- nrow(init)",
      "files <- list.files(adaptive_dir, pattern=paste0(\"job_.*_\", scenario_name, \"\\\\.csv\"), full.names=TRUE)",
      "new_list <- lapply(files, function(f) {",
      "  df <- read.csv(f, stringsAsFactors=FALSE)",
      "  if (nrow(df) > n0) df[(n0+1):nrow(df), ] else NULL",
      "})",
      "all <- do.call(rbind, c(list(init), new_list))",
      "all <- all[, names(init), drop=FALSE]",
      "",
      "out_dir <- file.path(output_dir, 'model_parameters')",
      "if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)",
      "final <- file.path(out_dir, paste0(scenario_name, '_model_parameters.csv'))",
      "write.csv(all, final, row.names=FALSE)",
      "file.remove(files)"
    )
    writeLines(gather_code, con=gather_R)
    Sys.chmod(gather_R, "0755")

    # 4) submit it (now the script itself carries the dependency)
    gather_job <- create_slurm_script(
      job_name    = paste0("gather_", scenario_name),
      script_path = gather_R,
      args        = c(results_file, adaptive_dir, scenario_name, output_dir),
      num_cores   = 1,
      time        = "00:10:00",
      memory      = "1G",
      output_file = file.path(adaptive_dir, "gather.out"),
      error_file  = file.path(adaptive_dir, "gather.err")
    )
    submit_job(gather_job, use_slurm=TRUE, cider=cider)
    if (verbose) cat("Gather job submitted with afterok:", dep, "\n")
    return(invisible(NULL))
  }


  # Local parallel execution
  temps <- vapply(seq_len(num_parallel_jobs), make_temp, FUN.VALUE = "")
  for (i in seq_along(temps)) file.copy(initial_samples_file, temps[i], overwrite = TRUE)
  if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(min(num_parallel_jobs, max_cores))
    parallel::clusterExport(cl, c("adaptive_MC_sampling", "distance_matrix", "mapping_max_iter",
                                 "relative_epsilon", "folds", "scenario_name", "output_dir"),
                              envir = environment())
    parallel::clusterEvalQ(cl, library(topolow))
    parallel::parLapply(cl, seq_along(temps), function(i) {
      adaptive_MC_sampling(samples_file    = temps[i],
                           distance_matrix = distance_matrix,
                           iterations      = iterations,
                           batch_size      = 1,
                           mapping_max_iter= mapping_max_iter,
                           relative_epsilon= relative_epsilon,
                           folds           = folds,
                           num_cores       = 1,
                           scenario_name   = scenario_name,
                           output_dir      = output_dir,
                           verbose         = FALSE)
    })
    parallel::stopCluster(cl)
  } else {
    parallel::mclapply(seq_along(temps), function(i) {
      adaptive_MC_sampling(samples_file    = temps[i],
                           distance_matrix = distance_matrix,
                           iterations      = iterations,
                           batch_size      = 1,
                           mapping_max_iter= mapping_max_iter,
                           relative_epsilon= relative_epsilon,
                           folds           = folds,
                           num_cores       = 1,
                           scenario_name   = scenario_name,
                           output_dir      = output_dir,
                           verbose         = FALSE)
    }, mc.cores = min(num_parallel_jobs, max_cores))
  }

  # Gather local results
  init2 <- read.csv(initial_samples_file, stringsAsFactors = FALSE)
  n0   <- nrow(init2)
  new_list <- lapply(temps, function(f) {
    df <- read.csv(f, stringsAsFactors = FALSE)
    if (nrow(df) > n0) df[(n0+1):nrow(df), ] else NULL
  })
  all <- do.call(rbind, c(list(init2), new_list))
  all <- all[, names(init2), drop=FALSE](rbind, c(list(init2), new_list))
  write.csv(all, results_file, row.names = FALSE)
  file.remove(temps)
  if (verbose) cat("Local jobs complete; results in", results_file, "\n")
  invisible(NULL)
}

adaptive_MC_sampling <- function(samples_file,
                                 distance_matrix,
                                 iterations = 1,
                                 batch_size = 1,
                                 mapping_max_iter,
                                 relative_epsilon,
                                 folds = 20,
                                 num_cores = 1,
                                 scenario_name,
                                 output_dir = NULL,
                                 verbose = FALSE) {
  # (unchanged; safe append with filelock)
}
adaptive_MC_sampling <- function(samples_file, 
                                 distance_matrix,
                                 iterations = 1, 
                                 batch_size = 1, # Now fixed to 1 by design
                                 mapping_max_iter, 
                                 relative_epsilon,
                                 folds = 20, 
                                 num_cores = 1,
                                 scenario_name, 
                                 output_dir = NULL,
                                 verbose = FALSE) {
  # Require filelock for safe concurrent writes
  if (!requireNamespace("filelock", quietly = TRUE)) {
    stop("Package 'filelock' is required for safe writes. Please install it.")
  }

  # Handle parallel processing setup
  use_parallelism <- num_cores > 1
  if (use_parallelism) {
    if (verbose) cat("Setting up parallel processing\n")
    if (.Platform$OS.type == "windows") {
      if (verbose) cat("Using parallel cluster for Windows\n")
      cl <- parallel::makeCluster(num_cores)
      on.exit(parallel::stopCluster(cl))
      parallel::clusterExport(cl, c("distance_matrix", "mapping_max_iter", 
                                   "relative_epsilon", "folds"),
                             envir = environment())
      parallel::clusterEvalQ(cl, { library(topolow) })
    }
  }

  # Handle output directory
  if (is.null(output_dir)) output_dir <- getwd()

  # Create results directory (not used for writes when using samples_file)
  param_dir <- file.path(output_dir, "model_parameters")
  if (!dir.exists(param_dir)) dir.create(param_dir, recursive = TRUE, showWarnings = FALSE)

  par_names <- c("log_N", "log_k0", "log_cooling_rate", "log_c_repulsion")
  key_cols <- c(par_names, "Holdout_MAE", "NLL")

  for (iter in seq_len(iterations)) {
    if (verbose) cat(sprintf("\nStarting iteration %d of %d\n", iter, iterations))

    # Read current samples (always from samples_file)
    current_samples <- read.csv(samples_file, stringsAsFactors = FALSE)
    for (col in key_cols[key_cols %in% names(current_samples)]) {
      current_samples[[col]] <- as.numeric(as.character(current_samples[[col]]))
    }
    current_samples <- current_samples[apply(current_samples, 1, 
                                           function(row) all(is.finite(row))), ]
    current_samples <- na.omit(current_samples)

    if (nrow(current_samples) == 0) {
      warning("No valid samples remaining after filtering")
      break
    }

    # Burn-in and convergence checks
    if (nrow(current_samples) > 10) {
      burn_in <- min(round(nrow(current_samples) * 0.3), nrow(current_samples) - 5)
      current_samples <- current_samples[-seq_len(burn_in), ]
    }
    if (nrow(current_samples) > 500) {
      conv_check <- check_gaussian_convergence(
        data = current_samples[, par_names],
        window_size = 500,
        tolerance = 0.002
      )
      if (conv_check$converged) {
        if (verbose) cat("Convergence achieved at iteration", iter, "\n")
        break
      }
    }

    # Generate new samples via KDE
    new_samples <- generate_kde_samples(samples = current_samples, n = batch_size)
    for (col in par_names) new_samples[[col]] <- as.numeric(new_samples[[col]])

    # Define likelihood evaluation
    evaluate_sample <- function(i) {
      N <- round(exp(new_samples[i, "log_N"]))
      k0 <- exp(new_samples[i, "log_k0"])
      cooling_rate <- exp(new_samples[i, "log_cooling_rate"])
      c_repulsion <- exp(new_samples[i, "log_c_repulsion"])
      inner_cores <- if (use_parallelism) 1 else min(folds, num_cores)
      tryCatch({
        likelihood_function(
          distance_matrix = distance_matrix,
          mapping_max_iter = mapping_max_iter,
          relative_epsilon = relative_epsilon,
          N = N,
          k0 = k0,
          cooling_rate = cooling_rate,
          c_repulsion = c_repulsion,
          folds = folds,
          num_cores = inner_cores
        )
      }, error = function(e) {
        if (verbose) cat("Error in likelihood calculation:", e$message, "\n")
        NA
      })
    }

    # Evaluate new samples
    if (use_parallelism) {
      if (.Platform$OS.type == "windows") {
        parallel::clusterExport(cl, c("evaluate_sample", "new_samples"), envir = environment())
        new_likelihoods <- parallel::parLapply(cl, seq_len(nrow(new_samples)), evaluate_sample)
      } else {
        new_likelihoods <- parallel::mclapply(seq_len(nrow(new_samples)), evaluate_sample, mc.cores = num_cores)
      }
    } else {
      new_likelihoods <- lapply(seq_len(nrow(new_samples)), evaluate_sample)
    }

    # Filter valid results
    valid_results <- !sapply(new_likelihoods, is.null) &
                     !sapply(new_likelihoods, function(x) all(is.na(unlist(x))))

    if (any(valid_results)) {
      likelihoods_mat <- do.call(rbind, lapply(new_likelihoods[valid_results], unlist))
      new_likelihoods_df <- as.data.frame(likelihoods_mat)
      colnames(new_likelihoods_df) <- c("Holdout_MAE", "NLL")

      valid_samples <- new_samples[valid_results, , drop = FALSE]
      valid_samples$Holdout_MAE <- as.numeric(new_likelihoods_df$Holdout_MAE)
      valid_samples$NLL <- as.numeric(new_likelihoods_df$NLL)

      for (col in key_cols[key_cols %in% names(valid_samples)]) {
        valid_samples[[col]] <- as.numeric(as.character(valid_samples[[col]]))
      }
      valid_samples <- valid_samples[apply(valid_samples, 1, function(row) all(is.finite(row))), ]

      if (nrow(valid_samples) > 0) {
        result_file <- samples_file
        lock_path <- paste0(result_file, ".lock")
        lock <- filelock::lock(lock_path)

        if (!file.exists(result_file)) {
          write.table(valid_samples, result_file,
                      sep = ",", row.names = FALSE, col.names = TRUE,
                      append = FALSE, quote = TRUE)
        } else {
          write.table(valid_samples, result_file,
                      sep = ",", row.names = FALSE, col.names = FALSE,
                      append = TRUE, quote = TRUE)
        }
        filelock::unlock(lock)

        if (verbose) cat(sprintf("Safely appended %d new valid samples to %s\n", 
                                 nrow(valid_samples), result_file))
      } else if (verbose) {
        cat("No valid samples in this iteration\n")
      }
    } else if (verbose) {
      cat("All likelihood evaluations failed in this iteration\n")
    }
  }

  # Return final samples
  result_file <- samples_file
  if (file.exists(result_file)) {
    final_samples <- read.csv(result_file, stringsAsFactors = FALSE)
    return(final_samples)
  } else {
    warning("No results file created")
    return(NULL)
  }
}
