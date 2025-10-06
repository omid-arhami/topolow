
process_param_set <- function(i) {
tryCatch({
    # Get parameters for this evaluation
    N <- lhs_params$N[i]
    k0 <- lhs_params$k0[i]
    c_repulsion <- lhs_params$c_repulsion[i]
    cooling_rate <- lhs_params$cooling_rate[i]
    
    if (verbose) {
    cat(sprintf("\nEvaluating sample %d/%d: N=%d, k0=%.3f, c_rep=%.4f, cool=%.4f\n",
                i, num_samples, N, k0, c_repulsion, cooling_rate))
    }
    
    # ======================================================================
    # SUBSAMPLING (if enabled)
    # ======================================================================
    local_matrix <- full_dissimilarity_matrix
    
    if (use_subsampling) {
    if (verbose) {
        cat(sprintf("  Subsampling to %d points...\n", opt_subsample))
    }
    
    subsample_result <- tryCatch({
        subsample_dissimilarity_matrix(
        dissimilarity_matrix = full_dissimilarity_matrix,
        sample_size = opt_subsample,
        max_attempts = 5,
        verbose = verbose
        )
    }, error = function(e) {
        # CHANGE: Track subsample failure with details
        failure_diagnostics$subsample_failures <<- failure_diagnostics$subsample_failures + 1
        failure_diagnostics$failure_messages <<- c(
        failure_diagnostics$failure_messages,
        sprintf("Sample %d: Subsample failed - %s", i, e$message)
        )
        if (verbose) {
        cat("  X Subsampling failed:", e$message, "\n")
        }
        return(NULL)
    })
    
    if (is.null(subsample_result)) {
        return(NULL)
    }

    if (!subsample_result$is_connected) {
        # Track connectivity failure with details
        failure_diagnostics$subsample_failures <<- failure_diagnostics$subsample_failures + 1
        failure_diagnostics$failure_messages <<- c(
        failure_diagnostics$failure_messages,
        sprintf("Sample %d: Disconnected subsample (%d components, %.1f%% complete)",
                i, subsample_result$n_components, 
                subsample_result$completeness * 100)
        )
        if (verbose) {
        cat("  X Failed to obtain connected subsample. Skipping this parameter set.\n")
        }
        return(NULL)
    }
    
    local_matrix <- subsample_result$subsampled_matrix
    
    # Sanity check the subsampled data
    sanity_result <- sanity_check_subsample(
        local_matrix,
        folds = folds,
        verbose = verbose
    )
    
    if (!sanity_result$all_checks_passed && verbose) {
        cat("  [!] Subsample sanity checks raised warnings (proceeding anyway)\n")
    }
    }
    
    # ======================================================================
    # CROSS-VALIDATION EVALUATION
    # ======================================================================
    result <- tryCatch({
    likelihood_function(
        dissimilarity_matrix = local_matrix,
        mapping_max_iter = mapping_max_iter,
        relative_epsilon = relative_epsilon,
        N = N,
        k0 = k0,
        cooling_rate = cooling_rate,
        c_repulsion = c_repulsion,
        folds = folds,
        num_cores = 1
    )
    }, error = function(e) {
    # Track CV/embedding failure with details
    if (grepl("fold|sparse", e$message, ignore.case = TRUE)) {
        failure_diagnostics$cv_fold_failures <<- failure_diagnostics$cv_fold_failures + 1
    } else {
        failure_diagnostics$embedding_failures <<- failure_diagnostics$embedding_failures + 1
    }
    failure_diagnostics$failure_messages <<- c(
        failure_diagnostics$failure_messages,
        sprintf("Sample %d: Evaluation failed - %s", i, 
                substr(e$message, 1, 100))
    )
    if (verbose) {
        cat("  X Evaluation error:", e$message, "\n")
    }
    return(list(Holdout_MAE = NA, NLL = NA))
    })
    
    # ======================================================================
    # CHECK MAE REASONABLENESS
    # ======================================================================
    # Check if result is valid
    if (is.null(result) || is.na(result$Holdout_MAE) || is.na(result$NLL)) {
    if (is.null(result)) {
        failure_diagnostics$other_failures <<- failure_diagnostics$other_failures + 1
        failure_diagnostics$failure_messages <<- c(
        failure_diagnostics$failure_messages,
        sprintf("Sample %d: Returned NULL result", i)
        )
    }
    return(NULL)
    }

    if (!is.na(result$Holdout_MAE) && !is.null(result$Holdout_MAE)) {
    # Get median of observed dissimilarities for comparison
    # (handles threshold indicators):
    observed_dissim <- extract_numeric_values(local_matrix)
    median_dissim <- median(observed_dissim[!is.na(observed_dissim)])
    
    if (!is.na(median_dissim) && result$Holdout_MAE > (2 * median_dissim)) {
        if (verbose) {
        cat(sprintf("  [!] High MAE (%.3f) vs median dissimilarity (%.3f). Poor fit suspected.\n",
                    result$Holdout_MAE, median_dissim))
        }
    }
    }
    
    # ======================================================================
    # COMPILE RESULTS
    # ======================================================================
    result_row <- data.frame(
    N = N,
    k0 = k0,
    cooling_rate = cooling_rate,
    c_repulsion = c_repulsion,
    Holdout_MAE = result$Holdout_MAE,
    NLL = result$NLL
    )
    
    # Add subsampling metadata if used
    if (use_subsampling) {
    result_row$opt_subsample <- opt_subsample
    result_row$original_n_points <- matrix_size
    }
    
    if (verbose) {
    cat(sprintf("  [OK] MAE: %.4f, NLL: %.2f\n",
                result$Holdout_MAE, result$NLL))
    }
    
    return(result_row)
    
}, error = function(e) {
    if (verbose) {
    cat(sprintf("  X Error evaluating sample %d: %s\n", i, e$message))
    }
    return(NULL)
})
}