
#' Fast Meta-Analysis
#' 
#' Quick meta-analysis without Bayesian components
#' 
#' @param data Data frame with yi and se columns
#' @return cbamm_fast object with results
#' @export
cbamm_fast <- function(data) {
  if (!is.data.frame(data)) {
    stop("Data must be a data.frame")
  }
  
  if (all(c("yi", "se") %in% names(data))) {
    weights <- 1 / data$se^2
    effect <- sum(data$yi * weights) / sum(weights)
    se <- sqrt(1 / sum(weights))
    Q <- sum(weights * (data$yi - effect)^2)
    df <- nrow(data) - 1
    I2 <- max(0, 100 * (Q - df) / Q)
    
    result <- list(
      summary = list(
        estimate = effect,
        se = se,
        ci_lower = effect - 1.96 * se,
        ci_upper = effect + 1.96 * se,
        I2 = I2,
        n_studies = nrow(data)
      ),
      data = data,
      method = "FAST",
      time_elapsed = 0
    )
    
    class(result) <- "cbamm_fast"
    return(result)
  } else {
    stop("Data must contain 'yi' and 'se' columns")
  }
}

#' @export
print.cbamm_fast <- function(x, ...) {
  cat("\nFast Meta-Analysis Results\n")
  cat("==========================\n")
  cat("Studies:", x$summary$n_studies, "\n")
  cat("Effect:", round(x$summary$estimate, 3), "\n")
  cat("95% CI: [", round(x$summary$ci_lower, 3), ", ", 
      round(x$summary$ci_upper, 3), "]\n", sep = "")
  cat("I² heterogeneity:", round(x$summary$I2, 1), "%\n")
  invisible(x)
}

