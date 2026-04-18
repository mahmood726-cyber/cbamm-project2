
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

#' Comprehensive Bayesian and Advanced Meta-Analysis Methods (Performance Optimized)
#'
#' @description 
#' High-performance implementation with transportability weighting, HKSJ adjustment, 
#' robust variance estimation, PET-PEESE, Bayesian methods, and conflict detection.
#' Optimized for speed while maintaining analytical rigor.
#'
#' @param data A data.frame containing study-level data with required columns:
#'   \code{study_id}, \code{yi} (effect sizes), \code{se} (standard errors),
#'   \code{study_type} (RCT/observational)
#' @param config Configuration object from \code{\link{cbamm_config}}
#' @param target_population List with transportability targets (if using transport weighting)
#' @param performance_mode Character: "fast" (default), "comprehensive", or "balanced"
#'
#' @return Object of class "cbamm_results"
#' @export
cbamm <- function(data, config = cbamm_config(), target_population = NULL, 
                  performance_mode = "fast") {
  
  start_time <- Sys.time()
  
  # Input validation
  validate_cbamm_data(data)
  
  # Apply performance optimizations based on mode
  config <- optimize_config_for_performance(config, performance_mode, data)
  
  if (config$output$verbose) {
    message("Starting CBAMM analysis (", performance_mode, " mode)...")
    message("Configuration: ", length(names(config$methods)[sapply(config$methods, isTRUE)]), " methods enabled")
    message("Dataset size: ", nrow(data), " studies")
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
        convergence_issues = character(0),
        optimizations_applied = character(0),
        performance_mode = performance_mode
      ),
      data_summary = summarize_input_data(data)
    ),
    class = "cbamm_results"
  )
  
  # Track optimizations applied
  optimizations <- character(0)
  
  # Add default weights if not present
  if (!"analysis_weights" %in% names(data)) {
    data$analysis_weights <- rep(1, nrow(data))
  }
  
  # OPTIMIZED: Core meta-analysis
  if (config$output$verbose) message("Running optimized core meta-analysis...")
  
  tryCatch({
    fit_transport <- fast_robust_rma(
      yi = data$yi, 
      sei = data$se, 
      data = data,
      method = config$estimators[1],
      weights = data$analysis_weights,
      use_hksj = config$methods$hksj
    )
    results$meta_results <- list(transport = fit_transport)
    optimizations <- c(optimizations, "fast_meta_analysis")
    
  }, error = function(e) {
    if (config$output$verbose) message("Fast meta-analysis failed, using fallback")
    results$meta_results <- safe_run_core_meta_analysis(data, config)
  })
  
  # Transport weighting
  if (config$methods$transport) {
    if (config$output$verbose) message("Computing transport weights...")
    tryCatch({
      results$transport_results <- run_transport_analysis(data, target_population, config)
      if (!is.null(results$transport_results) && "transport_weights" %in% names(results$transport_results)) {
        data$transport_weights <- results$transport_results$transport_weights
        data$analysis_weights <- results$transport_results$analysis_weights
      }
    }, error = function(e) {
      warning("Transport analysis failed: ", e$message)
      results$diagnostics$warnings <- c(results$diagnostics$warnings, paste("Transport:", e$message))
    })
  }
  
  # Conditional Bayesian analysis
  if (config$methods$bayesian) {
    if (should_skip_bayesian(data, min_studies = 6)) {
      if (config$output$verbose) message("Skipping Bayesian analysis (dataset too small)")
      optimizations <- c(optimizations, "skipped_bayesian_small_n")
    } else {
      if (config$output$verbose) message("Running optimized Bayesian analysis...")
      tryCatch({
        results$bayesian_results <- run_bayesian_analysis(data, config)
        optimizations <- c(optimizations, "optimized_bayesian")
      }, error = function(e) {
        warning("Bayesian analysis failed: ", e$message)
        results$diagnostics$errors <- c(results$diagnostics$errors, paste("Bayesian:", e$message))
      })
    }
  }
  
  # Optimized sensitivity analysis
  if (config$output$verbose) message("Running optimized sensitivity analyses...")
  
  tryCatch({
    results$sensitivity_results <- run_optimized_sensitivity_analyses(data, config, performance_mode)
    optimizations <- c(optimizations, "optimized_sensitivity")
  }, error = function(e) {
    warning("Sensitivity analysis failed: ", e$message)
    results$diagnostics$warnings <- c(results$diagnostics$warnings, paste("Sensitivity:", e$message))
  })
  
  # Conflict detection
  if (config$methods$conflict_detection) {
    if (config$output$verbose) message("Running conflict detection...")
    tryCatch({
      results$conflict_detection <- detect_study_conflicts(
        data = data, 
        threshold = 0.15,
        k_candidates = 2:4
      )
    }, error = function(e) {
      warning("Conflict detection failed: ", e$message)
      results$diagnostics$warnings <- c(results$diagnostics$warnings, paste("Conflict:", e$message))
    })
  }
  
  # Generate advisor recommendations
  if (config$output$verbose) message("Generating recommendations...")
  tryCatch({
    results$advisor_recommendations <- generate_advisor_recommendations(results, config)
  }, error = function(e) {
    warning("Advisor recommendations failed: ", e$message)
    results$diagnostics$warnings <- c(results$diagnostics$warnings, paste("Advisor:", e$message))
  })
  
  # Conditional plot generation
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
  results$diagnostics$optimizations_applied <- optimizations
  
  if (config$output$verbose) {
    message("CBAMM analysis completed in ", round(results$diagnostics$computation_time, 2), " seconds")
    
    if (length(optimizations) > 0) {
      message("Performance optimizations: ", paste(optimizations, collapse = ", "))
    }
  }
  
  return(results)
}
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

#' Summarize Input Data
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
