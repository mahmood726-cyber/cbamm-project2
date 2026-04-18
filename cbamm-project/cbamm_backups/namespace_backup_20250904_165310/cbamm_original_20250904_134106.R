
#' Comprehensive Bayesian and Advanced Meta-Analysis Methods
#'
#' @description 
#' Implements transportability weighting, HKSJ adjustment, robust variance
#' estimation, PET-PEESE, Bayesian methods, and novel conflict detection
#' in a unified meta-analysis framework.
#'
#' @param data A data.frame containing study-level data with required columns:
#'   \code{study_id}, \code{yi} (effect sizes), \code{se} (standard errors),
#'   \code{study_type} (RCT/observational)
#' @param config Configuration object from \code{\link{cbamm_config}}
#' @param target_population List with transportability targets (if using transport weighting)
#'
#' @return Object of class 'cbamm_results'
#' @export
cbamm <- function(data, config = cbamm_config(), target_population = NULL) {
  
  # Start timing
  start_time <- Sys.time()
  
  # Input validation
  validate_cbamm_data(data)
  
  if (config$output$verbose) {
    message("Starting CBAMM analysis...")
    message("Configuration: ", length(names(config$methods)[sapply(config$methods, isTRUE)]), " methods enabled")
  }
  
  # Initialize results structure
  results <- structure(
    list(
      meta_results = NULL,
      transport_results = NULL,
      bayesian_results = NULL,
      sensitivity_results = NULL,
      advisor_recommendations = NULL,
      conflict_detection = NULL,
      plots = list(),
      config = config,
      diagnostics = list(
        warnings = character(0),
        errors = character(0),
        computation_time = NULL,
        convergence_issues = character(0)
      ),
      data_summary = summarize_input_data(data)
    ),
    class = "cbamm_results"
  )
  
  # Core meta-analysis
  if (config$output$verbose) message("Running core meta-analysis...")
  results$meta_results <- safe_run_core_meta_analysis(data, config)
  
  # Transport weighting (if enabled)
  if (config$methods$transport) {
    if (config$output$verbose) message("Computing transport weights...")
    tryCatch({
      results$transport_results <- run_transport_analysis(data, target_population, config)
      # Update data with transport weights for subsequent analyses
      if (!is.null(results$transport_results) && "transport_weights" %in% names(results$transport_results)) {
        data$transport_weights <- results$transport_results$transport_weights
        data$analysis_weights <- results$transport_results$analysis_weights
      }
    }, error = function(e) {
      warning("Transport analysis failed: ", e$message)
      results$diagnostics$warnings <- c(results$diagnostics$warnings, paste("Transport:", e$message))
    })
  }
  
  # Add default weights if not present
  if (!"analysis_weights" %in% names(data)) {
    data$analysis_weights <- rep(1, nrow(data))
  }
  
  # Conflict detection and advisor (FIXED parameter passing)
  if (config$methods$conflict_detection) {
    if (config$output$verbose) message("Running conflict detection...")
    tryCatch({
      # Fixed: Pass individual parameters instead of config object
      results$conflict_detection <- detect_study_conflicts(
        data = data, 
        threshold = 0.15,  # Extract from config or use default
        k_candidates = 2:4  # Extract from config or use default
      )
    }, error = function(e) {
      warning("Conflict detection failed: ", e$message)
      results$diagnostics$warnings <- c(results$diagnostics$warnings, paste("Conflict:", e$message))
    })
  }
  
  # Bayesian analysis (if enabled)  
  if (config$methods$bayesian) {
    if (config$output$verbose) message("Running Bayesian analysis...")
    tryCatch({
      results$bayesian_results <- run_bayesian_analysis(data, config)
    }, error = function(e) {
      warning("Bayesian analysis failed: ", e$message)
      results$diagnostics$errors <- c(results$diagnostics$errors, 
                                    paste("Bayesian:", e$message))
    })
  }
  
  # Sensitivity analyses
  if (config$output$verbose) message("Running sensitivity analyses...")
  tryCatch({
    results$sensitivity_results <- run_sensitivity_analyses(data, config)
  }, error = function(e) {
    warning("Sensitivity analysis failed: ", e$message)
    results$diagnostics$warnings <- c(results$diagnostics$warnings, paste("Sensitivity:", e$message))
  })
  
  # Generate advisor recommendations
  if (config$output$verbose) message("Generating recommendations...")
  tryCatch({
    results$advisor_recommendations <- generate_advisor_recommendations(results, config)
  }, error = function(e) {
    warning("Advisor recommendations failed: ", e$message)
    results$diagnostics$warnings <- c(results$diagnostics$warnings, paste("Advisor:", e$message))
  })
  
  # Generate plots
  if (config$output$plots) {
    if (config$output$verbose) message("Creating visualizations...")
    tryCatch({
      results$plots <- generate_cbamm_plots(results, config)
    }, error = function(e) {
      warning("Plot generation failed: ", e$message)
      results$diagnostics$warnings <- c(results$diagnostics$warnings, paste("Plots:", e$message))
    })
  }
  
  # Final diagnostics
  results$diagnostics$computation_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  if (config$output$verbose) {
    message("CBAMM analysis completed in ", 
            round(results$diagnostics$computation_time, 2), " seconds")
    
    # Report any issues
    if (length(results$diagnostics$warnings) > 0) {
      message("Warnings encountered: ", length(results$diagnostics$warnings))
    }
    if (length(results$diagnostics$errors) > 0) {
      message("Errors encountered: ", length(results$diagnostics$errors))
    }
  }
  
  return(results)
}

# S3 Methods for cbamm_results

#' Print method for cbamm_results
#' @export
print.cbamm_results <- function(x, ...) {
  cat("CBAMM Meta-Analysis Results\n")
  cat("===========================\n\n")
  
  # Data summary
  cat("Studies analyzed:", x$data_summary$n_studies, "\n")
  if (!is.null(x$data_summary$study_types)) {
    cat("Study types:", paste(names(x$data_summary$study_types), 
                            x$data_summary$study_types, 
                            sep = "=", collapse = ", "), "\n")
  }
  
  # Main results
  if (!is.null(x$meta_results$transport)) {
    cat("\n--- Transport-Weighted Results ---\n")
    report_meta_result(x$meta_results$transport, "Primary Analysis")
  }
  
  if (!is.null(x$meta_results$grade)) {
    cat("\n--- GRADE-Weighted Results ---\n")
    report_meta_result(x$meta_results$grade, "GRADE Analysis")
  }
  
  # Advisor recommendations
  if (!is.null(x$advisor_recommendations) && length(x$advisor_recommendations$recommendations) > 0) {
    cat("\n--- CBAMM Advisor Recommendations ---\n")
    for (rec in x$advisor_recommendations$recommendations) {
      cat(" •", rec, "\n")
    }
  }
  
  # Diagnostics
  if (length(x$diagnostics$warnings) > 0) {
    cat("\nWarnings:", length(x$diagnostics$warnings), "\n")
  }
  
  cat("\nComputation time:", round(x$diagnostics$computation_time, 2), "seconds\n")
  
  invisible(x)
}

#' Summary method for cbamm_results
#' @export
summary.cbamm_results <- function(object, ...) {
  print(object)
  
  # Additional detailed output
  if (!is.null(object$conflict_detection)) {
    cat("\n--- Conflict Detection ---\n")
    cat("Clustering method: K-means with K =", object$conflict_detection$K, "\n")
    cat("Effect size difference:", round(object$conflict_detection$delta, 3), "\n")
    if (object$conflict_detection$threshold_met) {
      cat("⚠ Significant conflicts detected between study clusters\n")
    } else {
      cat("✓ No significant conflicts detected\n")
    }
  }
  
  if (!is.null(object$transport_results)) {
    cat("\n--- Transport Analysis ---\n")
    if (!is.null(object$transport_results$balance)) {
      cat("Age balance improvement:", 
          round(object$transport_results$balance$pre_age, 1), "→",
          round(object$transport_results$balance$post_age, 1), "\n")
    }
  }
  
  invisible(object)
}

#' Plot method for cbamm_results
#' @export
plot.cbamm_results <- function(x, which = "all", ...) {
  if (length(x$plots) == 0) {
    message("No plots available. Set config$output$plots = TRUE to generate plots.")
    return(invisible(NULL))
  }
  
  if (which == "all") {
    # Display all available plots
    for (plot_name in names(x$plots)) {
      if (!is.null(x$plots[[plot_name]])) {
        message("Displaying:", plot_name)
        print(x$plots[[plot_name]])
      }
    }
  } else {
    # Display specific plot
    if (which %in% names(x$plots) && !is.null(x$plots[[which]])) {
      print(x$plots[[which]])
    } else {
      message("Plot '", which, "' not available. Available plots: ", 
              paste(names(x$plots), collapse = ", "))
    }
  }
  
  invisible(x)
}



#' Validate CBAMM Input Data.
#' 
#' Validates that input data has required columns and proper format
#' 
#' @param data Input data frame
#' @keywords internal
validate_cbamm_data <- function(data) {
  
  # Check data frame
  if (!is.data.frame(data)) {
    stop("data must be a data.frame")
  }
  
  if (nrow(data) < 3) {
    stop("Minimum of 3 studies required for meta-analysis")
  }
  
  # Required columns
  required <- c("study_id", "yi", "se")
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "), 
         "\nRequired: study_id, yi (effect sizes), se (standard errors)")
  }
  
  # Check for missing values in critical columns
  critical_na <- sapply(data[required], function(x) any(is.na(x)))
  if (any(critical_na)) {
    stop("Missing values detected in: ", paste(names(critical_na)[critical_na], collapse = ", "))
  }
  
  # Check numeric columns
  if (!is.numeric(data$yi)) stop("yi (effect sizes) must be numeric")
  if (!is.numeric(data$se)) stop("se (standard errors) must be numeric")
  if (any(data$se <= 0)) stop("All standard errors must be positive")
  
  # Check for extreme values that might cause issues
  if (any(abs(data$yi) > 10)) {
    warning("Large effect sizes detected (|yi| > 10). Consider checking data.")
  }
  
  if (any(data$se > 5)) {
    warning("Large standard errors detected (se > 5). Consider checking data.")
  }
  
  # Study type (optional but recommended)
  if ("study_type" %in% names(data)) {
    valid_types <- c("RCT", "observational", "cohort", "case_control", "cross_sectional")
    invalid_types <- setdiff(unique(data$study_type), valid_types)
    if (length(invalid_types) > 0) {
      warning("Unrecognized study types: ", paste(invalid_types, collapse = ", "))
    }
  }
  
  invisible(TRUE)
}

#' Summarize Input Data.
#' 
#' Creates summary of input data for results object
#' 
#' @param data Input data frame
#' @return List with data summary
#' @keywords internal
summarize_input_data <- function(data) {
  list(
    n_studies = nrow(data),
    effect_size_range = range(data$yi),
    se_range = range(data$se),
    study_types = if ("study_type" %in% names(data)) table(data$study_type) else NULL,
    sample_sizes = if ("n" %in% names(data)) summary(data$n) else NULL,
    has_covariates = sum(!names(data) %in% c("study_id", "yi", "se", "vi", "n", "study_type")) > 0
  )
}
