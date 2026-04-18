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
#' @return Object of class 'cbamm_results' containing:
#' \describe{
#'   \item{meta_results}{Traditional and robust meta-analysis results}
#'   \item{transport_results}{Transportability-weighted results (if enabled)}
#'   \item{bayesian_results}{Bayesian meta-analysis results (if enabled)}
#'   \item{sensitivity_results}{Publication bias and sensitivity analyses}
#'   \item{advisor_recommendations}{CBAMM advisor recommendations}
#'   \item{plots}{Generated plots}
#'   \item{config}{Configuration used}
#'   \item{diagnostics}{Model diagnostics and warnings}
#' }
#'
#' @examples
#' \dontrun{
#' # Simulate example data
#' data <- simulate_cbamm_data(n_rct = 10, n_obs = 15)
#' 
#' # Run analysis with default settings
#' results <- cbamm(data)
#' 
#' # Custom configuration
#' config <- cbamm_config(
#'   methods = list(transport = TRUE, bayesian = FALSE),
#'   estimators = c("REML", "PM")
#' )
#' results <- cbamm(data, config = config)
#' 
#' # With target population for transportability
#' target_pop <- list(age_mean = 65, female_prop = 0.6)
#' results <- cbamm(data, target_population = target_pop)
#' }
#' 
#' @references
#' Hartung, J., & Knapp, G. (2001). A refined method for the meta‐analysis of controlled clinical trials with binary outcome. Statistics in Medicine, 20(24), 3875-3889.
#' 
#' Egger, M., Smith, G. D., Schneider, M., & Minder, C. (1997). Bias in meta-analysis detected by a simple, graphical test. BMJ, 315(7109), 629-634.
#' 
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
  results$meta_results <- run_core_meta_analysis(data, config)
  
  # Transport weighting (if enabled)
  if (config$methods$transport) {
    if (config$output$verbose) message("Computing transport weights...")
    results$transport_results <- run_transport_analysis(data, target_population, config)
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
  results$sensitivity_results <- run_sensitivity_analyses(data, config)
  
  # Conflict detection and advisor
  if (config$methods$conflict_detection) {
    if (config$output$verbose) message("Running conflict detection...")
    results$conflict_detection <- detect_study_conflicts(data)
  }
  
  # Generate advisor recommendations
  if (config$output$verbose) message("Generating recommendations...")
  results$advisor_recommendations <- generate_advisor_recommendations(results, config)
  
  # Generate plots
  if (config$output$plots) {
    if (config$output$verbose) message("Creating visualizations...")
    results$plots <- generate_cbamm_plots(results, config)
  }
  
  # Final diagnostics
  results$diagnostics$computation_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  if (config$output$verbose) {
    message("CBAMM analysis completed in ", 
            round(results$diagnostics$computation_time, 2), " seconds")
  }
  
  return(results)
}

#' Validate CBAMM Input Data
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
    valid_types <- c("RCT", "OBS", "MR", "observational", "cohort", "case_control", "cross_sectional")
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
