# Copyright (c) 2024 Omid Arhami omid.arhami@uga.edu
# Licensed under MIT License
# R/jobs.R

#' Job Submission Utilities
#' 
#' @description
#' Functions for submitting and managing computational jobs using SLURM
#' or running them locally. Provides consistent interface regardless of 
#' execution environment.
#'
#' @keywords internal
"_PACKAGE"

#' Create SLURM Job Script
#'
#' @description
#' Creates a SLURM batch script with specified parameters and resource requests.
#'
#' @param r_module Character. R module to load (default: "R/4.3.2-foss-2022b")
#' @param job_name Name of the job
#' @param script_path Path to R script to execute
#' @param args Vector of command line arguments
#' @param num_cores Number of CPU cores to request
#' @param output_file Path for job output file
#' @param error_file Path for job error file
#' @param time Time limit (default: "8:00:00")
#' @param memory Memory request (default: "14G")
#' @param partition SLURM partition (default: "rohani_p")
#' @return Path to created script file
#' @export
create_slurm_script <- function(job_name, script_path, args, num_cores,
                                output_file, error_file,
                                time = "8:00:00", memory = "14G",
                                partition = "rohani_p", 
                                r_module = "R/4.3.2-foss-2022b") {
  if (!has_slurm()) {
    warning("SLURM appears to not be available on this system")
  }

  if (!file.exists(script_path)) {
    stop("script_path does not exist: ", script_path)
  }

  # Make module loading configurable
  module_cmd <- if (!is.null(r_module)) {
    paste0("module load ", r_module, "\n")
  } else {
    ""  # No module loading if NULL
  }

  slurm_script <- paste0(
    "#!/bin/bash\n",
    "#SBATCH --partition=", partition, "\n",
    "#SBATCH --job-name=", job_name, "\n",
    "#SBATCH --output=", output_file, "\n", 
    "#SBATCH --error=", error_file, "\n",
    "#SBATCH --ntasks=1\n",
    "#SBATCH --cpus-per-task=", num_cores, "\n",
    "#SBATCH --time=", time, "\n",
    "#SBATCH --mem=", memory, "\n\n",
     module_cmd,  # Use configured module
    "Rscript ", script_path, " ", paste(args, collapse = " "), "\n"
  )

  script_file <- tempfile(fileext = ".sh")
  writeLines(slurm_script, script_file)

  return(script_file)
}


#' Submit Job to SLURM or Run Locally
#'
#' @description
#' Submits a job to SLURM if available, otherwise runs locally.
#' Provides consistent interface for both execution modes.
#'
#' @param script_file Path to script file 
#' @param use_slurm Logical; whether to use SLURM if available
#' @param cider Logical; whether to use cider_qos queue
#' @return Exit status code (invisible)
#' @export
submit_job <- function(script_file, use_slurm = TRUE, cider = FALSE) {
  if (!file.exists(script_file)) {
    stop("script_file does not exist: ", script_file)
  }
  
  if (use_slurm && has_slurm()) {
    cmd <- if(cider) {
      paste("sbatch -q cider_qos -p batch", script_file)
    } else {
      paste("sbatch", script_file)
    }
    status <- system(cmd)
  } else {
    status <- system(paste("Rscript", script_file))
  }
  
  return(invisible(status))
}


#' Check if SLURM is Available
#' 
#' @return Logical indicating if SLURM commands are available
#' @keywords internal
has_slurm <- function() {
  system("which sbatch", ignore.stdout = TRUE) == 0
}


#' Check Status of Submitted Job
#' @param job_id Character. SLURM job ID
#' @return Character job status or NA if not found
#' @export
check_job_status <- function(job_id) {
  if (!has_slurm()$available) {
    warning("SLURM not available")
    return(NA)
  }
  
  status <- system(
    paste("squeue -h -j", job_id, "-o %t"), 
    intern = TRUE,
    ignore.stderr = TRUE
  )
  
  if (length(status) == 0) {
    # Check if job completed
    sacct <- system(
      paste("sacct -j", job_id, "-o State -n"),
      intern = TRUE,
      ignore.stderr = TRUE
    )
    if (length(sacct) > 0) {
      return(trimws(sacct[1]))
    }
    return(NA)
  }
  
  return(status)
}