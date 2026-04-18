#' Package configuration settings..
#' 
#' Creates configuration object for CBAMM analysis
#'
#' @param methods List of methods to enable/disable
#' @param estimators Vector of variance estimators to try
#' @param output List of output options  
#' @param transportability List of transportability settings
#' @param bayesian List of Bayesian method settings
#' @param sensitivity List of sensitivity analysis settings
#'
#' @return List containing configuration options
#'
#' @examples
#' # Default configuration
#' config <- cbamm_config()
#' 
#' # Custom configuration
#' config <- cbamm_config(
#'   methods = list(transport = FALSE, hksj = TRUE),
#'   estimators = c("REML", "DL")
#' )
cbamm_config <- function(
  methods = list(
    transport = TRUE,
    hksj = TRUE,
    bayesian = TRUE,
    robust_variance = TRUE,
    pet_peese = TRUE,
    conflict_detection = TRUE,
    missing_studies = TRUE,
    grade_weighting = FALSE
  ),
  estimators = c("REML", "DL", "PM", "ML"),
  output = list(
    export = FALSE,
    interactive = FALSE,
    verbose = TRUE,
    plots = TRUE
  ),
  transportability = list(
    method = "entropy",  # or "covariate_balance"
    target_population = NULL,
    covariates = NULL
  ),
  bayesian = list(
    method = "brms",  # or "jags"
    chains = 4,
    iter = 2000,
    warmup = 1000,
    cores = 1,
    prior = "weakly_informative"
  ),
  sensitivity = list(
    egger_test = TRUE,
    begg_test = TRUE,
    trim_fill = TRUE,
    copas_selection = FALSE,
    publication_bias_methods = c("egger", "begg", "trim_fill")
  )
) {
  
  # Validate configuration
  validate_config(methods, estimators, output, transportability, bayesian, sensitivity)
  
  structure(
    list(
      methods = methods,
      estimators = estimators,
      output = output,
      transportability = transportability,
      bayesian = bayesian,
      sensitivity = sensitivity,
      timestamp = Sys.time()
    ),
    class = "cbamm_config"
  )
}

#' Validate Package configuration settings..
#' 
#' Internal function to validate configuration parameters
#' @keywords internal
validate_config <- function(methods, estimators, output, transportability, bayesian, sensitivity) {
  
  # Validate estimators
  valid_estimators <- c("REML", "DL", "PM", "ML", "EB", "SJ")
  invalid_est <- setdiff(estimators, valid_estimators)
  if (length(invalid_est) > 0) {
    stop("Invalid estimators: ", paste(invalid_est, collapse = ", "), 
         "\nValid options: ", paste(valid_estimators, collapse = ", "))
  }
  
  # Validate Bayesian method
  if (methods$bayesian) {
    if (!bayesian$method %in% c("brms", "jags")) {
      stop("bayesian$method must be 'brms' or 'jags'")
    }
    
    if (bayesian$chains < 1 || bayesian$iter < 100) {
      stop("Bayesian settings: chains >= 1, iter >= 100")
    }
  }
  
  # Validate transportability settings
  if (!is.null(methods$transport) && isTRUE(methods$transport)) {
    if (!transportability$method %in% c("entropy", "covariate_balance", "ipsw")) {
      stop("transportability$method must be 'entropy', 'covariate_balance', or 'ipsw'")
    }
  }
  
  invisible(TRUE)
}

#' Print method for cbamm_config
print.cbamm_config <- function(x, ...) {
  cat("CBAMM Configuration\n")
  cat("===================\n\n")
  
  cat("Enabled Methods:\n")
  enabled <- names(x$methods)[sapply(x$methods, isTRUE)]
  cat(" ", paste(enabled, collapse = ", "), "\n\n")
  
  cat("Variance Estimators:", paste(x$estimators, collapse = ", "), "\n")
  
  if (x$methods$bayesian) {
    cat("Bayesian Method:", x$bayesian$method, "\n")
  }
  
  if (x$methods$transport) {
    cat("Transport Method:", x$transportability$method, "\n")
  }
  
  cat("\nCreated:", format(x$timestamp), "\n")
  invisible(x)
}

#' Update Package configuration settings..
#' 
#' Update an existing configuration object
#' 
#' @param config Existing cbamm_config object
#' @param ... Named arguments to update
#' 
#' @return Updated configuration object
update_config <- function(config, ...) {
  if (!inherits(config, "cbamm_config")) {
    stop("config must be a cbamm_config object")
  }
  
  updates <- list(...)
  
  # Update nested lists properly
  for (name in names(updates)) {
    if (name %in% names(config) && is.list(config[[name]]) && is.list(updates[[name]])) {
      config[[name]] <- modifyList(config[[name]], updates[[name]])
    } else {
      config[[name]] <- updates[[name]]
    }
  }
  
  # Re-validate
  do.call(validate_config, config[c("methods", "estimators", "output", 
                                   "transportability", "bayesian", "sensitivity")])
  
  config$timestamp <- Sys.time()
  config
}


#' Helper to properly merge config methods
#' @keywords internal
merge_methods <- function(defaults, updates) {
  # Start with all default methods
  result <- defaults
  
  # Update with user-specified values, preserving FALSE values
  if (!is.null(updates)) {
    for (name in names(updates)) {
      result[[name]] <- updates[[name]]
    }
  }
  
  return(result)
}

