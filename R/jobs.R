# Copyright (c) 2024 Omid Arhami o.arhami@gmail.com
# Licensed under Pre-Publication Academic License https://github.com/omid-arhami/topolow/blob/main/LICENSE

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
#' @param r_module Character. R module to load (default: "R/4.4.1-foss-2022b")
#' @param job_name Name of the job
#' @param script_path Path to R script to execute
#' @param args Vector of command line arguments
#' @param num_cores Number of CPU cores to request
#' @param output_file Path for job output file
#' @param error_file Path for job error file
#' @param time Time limit (default: "8:00:00")
#' @param memory Memory request (default: "14G")
#' @param partition SLURM partition (default: "rohani_p")
#' @param working_dir Working directory (default: NULL)
#' @param extra_sbatch_args Additional SBATCH arguments (default: NULL)
#' @return Path to created script file
#' @export
create_slurm_script <- function(job_name,
                                script_path,
                                args,
                                num_cores,
                                output_file,
                                error_file,
                                time = "8:00:00",
                                memory = "4G",
                                partition = "rohani_p",
                                r_module = "R/4.4.1-foss-2022b",
                                working_dir = NULL,
                                extra_sbatch_args = NULL) {
    if (!has_slurm()) {
        warning("SLURM appears to not be available on this system")
    }

    if (!file.exists(script_path)) {
        stop("script_path does not exist: ", script_path)
    }

    # Build SBATCH header
    sbatch_header <- c(
        "#!/bin/bash",
        paste0("#SBATCH --partition=", partition),
        paste0("#SBATCH --job-name=", job_name),
        paste0("#SBATCH --output=", output_file),
        paste0("#SBATCH --error=", error_file),
        "#SBATCH --ntasks=1",
        paste0("#SBATCH --cpus-per-task=", num_cores),
        paste0("#SBATCH --time=", time),
        paste0("#SBATCH --mem=", memory)
    )

    if (!is.null(working_dir)) {
        sbatch_header <- c(sbatch_header,
                           paste0("#SBATCH --chdir=", working_dir))
    }
    if (!is.null(extra_sbatch_args)) {
        sbatch_header <- c(sbatch_header, extra_sbatch_args)
    }

    module_cmd <- if (!is.null(r_module)) paste0("module load ", r_module) else ""
    # Package installation checks
    package_cmds <- paste(
        "# Load required packages",
        "Rscript -e \"if(!require(reshape2)) install.packages('reshape2', repos='https://cloud.r-project.org')\"",
        "Rscript -e \"if(!require(data.table)) install.packages('data.table', repos='https://cloud.r-project.org')\"",
        "Rscript -e \"if(!require(dplyr)) install.packages('dplyr', repos='https://cloud.r-project.org')\"",
        sep = "
    "
    )
    # --- START OF FIX ---
    # Process arguments: quote them unless they are shell variables (start with $)
    processed_args <- vapply(args, function(arg) {
        if (is.character(arg) && startsWith(arg, "$")) {
            return(arg) # Don't quote shell variables
        } else {
            return(shQuote(arg)) # Quote all other arguments
        }
    }, FUN.VALUE = character(1))

    # Assemble the final command
    r_command <- paste("Rscript", shQuote(script_path), paste(processed_args, collapse = " "))
    # --- END OF FIX ---

    # Final script assembly
    slurm_script <- paste0(
        paste(sbatch_header, collapse = "\n"), "\n\n",
        module_cmd, "\n\n",
        r_command, "\n"
    )

    script_file <- tempfile(fileext = ".sh")
    writeLines(slurm_script, script_file)
    Sys.chmod(script_file, "0755")
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
    # Build sbatch arguments
    sbatch_args <- character()
    if (cider) {
      sbatch_args <- c(sbatch_args, "-q", "cider_qos", "-p", "batch")
    }
    # Ask sbatch to echo back the job ID
    sbatch_args <- c(sbatch_args, "--parsable", script_file)
    
    # system2 will return the stdout (the job ID) as a character vector
    job_id <- system2("sbatch", sbatch_args, stdout = TRUE)
    return(job_id)
    
  } else {
    # Fallback: run the script locally with Rscript
    status <- system2("Rscript", script_file)
    return(status)
  }
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
  if (!has_slurm()) {
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