


# Then AFTER parallel execution, aggregate the failures:

# Process batch in parallel (Windows example)
if (.Platform$OS.type == "windows") {
  cl <- parallel::makeCluster(min(cores_to_use, length(batch_indices)))
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  parallel::clusterEvalQ(cl, library(topolow))
  parallel::clusterExport(cl, c(...), envir = environment())
  
  batch_results <- parallel::parLapply(cl, batch_indices, process_param_set)
} else {
  batch_results <- parallel::mclapply(batch_indices, process_param_set,
                                      mc.cores = min(cores_to_use, length(batch_indices)))
}

# AGGREGATE FAILURE DIAGNOSTICS (this is the key fix)
failure_diagnostics <- list(
  subsample_failures = 0,
  cv_fold_failures = 0,
  embedding_failures = 0,
  other_failures = 0,
  failure_messages = c()
)

for (res in batch_results) {
  if (!is.null(res$failure_type)) {
    # Count failure by type
    if (res$failure_type == "subsample") {
      failure_diagnostics$subsample_failures <- failure_diagnostics$subsample_failures + 1
    } else if (res$failure_type == "cv_fold") {
      failure_diagnostics$cv_fold_failures <- failure_diagnostics$cv_fold_failures + 1
    } else if (res$failure_type == "embedding") {
      failure_diagnostics$embedding_failures <- failure_diagnostics$embedding_failures + 1
    } else {
      failure_diagnostics$other_failures <- failure_diagnostics$other_failures + 1
    }
    
    # Store message
    failure_diagnostics$failure_messages <- c(
      failure_diagnostics$failure_messages,
      sprintf("Sample %d: %s", res$sample_id, res$failure_message)
    )
  }
}

# Now extract valid results
valid_results <- Filter(function(x) !is.null(x$result), batch_results)
valid_results <- lapply(valid_results, function(x) x$result)