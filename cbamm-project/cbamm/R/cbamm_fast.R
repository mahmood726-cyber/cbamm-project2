
#' Fast Meta-Analysis
#' 
#' Quick meta-analysis without Bayesian components
#' 
#' @param data Data frame with yi and se columns
#' @return cbamm_fast object with results
#' @export
cbamm_fast <- function(data,
                       config = cbamm_config(),
                       target_population = NULL,
                       method = NULL,
                       ...) {
  normalized_data <- .normalize_cbamm_public_data(data)

  if (!is.null(method)) {
    config$estimators[1] <- method
  }

  result <- cbamm(
    data = normalized_data,
    config = config,
    target_population = target_population,
    performance_mode = "fast"
  )

  class(result) <- unique(c("cbamm_fast", class(result)))
  result
}

#' @export
print.cbamm_fast <- function(x, ...) {
  if (inherits(x, "cbamm_results")) {
    return(print.cbamm_results(x, ...))
  }

  cat("\nFast Meta-Analysis Results\n")
  cat("==========================\n")
  cat("Studies:", x$summary$n_studies, "\n")
  cat("Effect:", round(x$summary$estimate, 3), "\n")
  cat("95% CI: [", round(x$summary$ci_lower, 3), ", ", 
      round(x$summary$ci_upper, 3), "]\n", sep = "")
  cat("I² heterogeneity:", round(x$summary$I2, 1), "%\n")
  invisible(x)
}

