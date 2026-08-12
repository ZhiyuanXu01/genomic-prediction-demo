# Reproduce the complete genomic-prediction analysis.
#
# Run from the project root with:
# Rscript run_all.R

project_file <- "genomic-prediction-demo.Rproj"

if (!file.exists(project_file)) {
  stop(
    paste0(
      "Run this script from the project root, ",
      "which must contain ",
      project_file,
      "."
    ),
    call. = FALSE
  )
}

required_packages <- c(
  "BGLR",
  "ggplot2",
  "here",
  "rrBLUP"
)

package_available <- vapply(
  required_packages,
  requireNamespace,
  logical(1),
  quietly = TRUE
)

missing_packages <- required_packages[
  !package_available
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing required packages: ",
      paste(
        missing_packages,
        collapse = ", "
      ),
      ". Run renv::restore() before executing ",
      "the analysis."
    ),
    call. = FALSE
  )
}

analysis_scripts <- file.path(
  "R",
  c(
    "01_load_and_qc.R",
    "02_relationship_matrices.R",
    "03_fit_blup.R",
    "04_cross_validation.R",
    "05_partial_genotyping.R",
    "06_summarise_results.R"
  )
)

missing_scripts <- analysis_scripts[
  !file.exists(analysis_scripts)
]

if (length(missing_scripts) > 0) {
  stop(
    paste0(
      "Missing analysis scripts: ",
      paste(
        missing_scripts,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}

dir.create(
  "results",
  showWarnings = FALSE
)

dir.create(
  "figures",
  showWarnings = FALSE
)

analysis_environment <- new.env(
  parent = globalenv()
)

timing_list <- vector(
  "list",
  length(analysis_scripts)
)

analysis_start_time <- Sys.time()

cat(
  "\n========================================\n",
  "Genomic Prediction Demo\n",
  "Complete reproducible analysis\n",
  "========================================\n",
  sep = ""
)

for (
  script_index in
  seq_along(analysis_scripts)
) {
  script_path <-
    analysis_scripts[[script_index]]
  
  script_start_time <- Sys.time()
  
  cat(
    "\n----------------------------------------\n",
    sprintf(
      "Running %s (%d/%d)\n",
      script_path,
      script_index,
      length(analysis_scripts)
    ),
    "----------------------------------------\n",
    sep = ""
  )
  
  tryCatch(
    source(
      script_path,
      local = analysis_environment,
      echo = FALSE
    ),
    error = function(error_condition) {
      stop(
        paste0(
          "Analysis failed in ",
          script_path,
          ": ",
          conditionMessage(
            error_condition
          )
        ),
        call. = FALSE
      )
    }
  )
  
  script_elapsed_seconds <- as.numeric(
    difftime(
      Sys.time(),
      script_start_time,
      units = "secs"
    )
  )
  
  timing_list[[script_index]] <- data.frame(
    script = script_path,
    elapsed_seconds =
      script_elapsed_seconds,
    status = "completed",
    stringsAsFactors = FALSE
  )
  
  cat(
    sprintf(
      "\nCompleted %s in %.1f seconds.\n",
      script_path,
      script_elapsed_seconds
    )
  )
}

run_timing <- do.call(
  rbind,
  timing_list
)

rownames(run_timing) <- NULL

total_elapsed_seconds <- as.numeric(
  difftime(
    Sys.time(),
    analysis_start_time,
    units = "secs"
  )
)

stopifnot(
  nrow(run_timing) ==
    length(analysis_scripts),
  all(
    run_timing$status ==
      "completed"
  ),
  file.exists(
    file.path(
      "results",
      "final_model_performance.csv"
    )
  ),
  file.exists(
    file.path(
      "results",
      "final_paired_improvements.csv"
    )
  ),
  file.exists(
    file.path(
      "figures",
      "repeated_cv_model_comparison.png"
    )
  ),
  file.exists(
    file.path(
      "figures",
      "partial_genotyping_model_comparison.png"
    )
  )
)

cat(
  "\n========================================\n",
  "Analysis completed successfully\n",
  "========================================\n",
  sep = ""
)

print(run_timing)

cat(
  sprintf(
    "\nTotal elapsed time: %.1f seconds\n",
    total_elapsed_seconds
  )
)

cat(
  paste0(
    "Final tables: results/",
    "final_model_performance.csv and results/",
    "final_paired_improvements.csv\n"
  )
)

