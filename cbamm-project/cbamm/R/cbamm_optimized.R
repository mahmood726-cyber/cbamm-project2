
#' Optimized CBAMM Analysis
#' 
#' Fast meta-analysis with configurable modes
#' 
#' @param data Data frame with effect sizes
#' @param mode Analysis mode: "quick", "standard", or "full"
#' @param effect_col Name of effect size column
#' @param se_col Name of standard error column
#' @param study_col Name of study ID column
#' 
#' @return cbamm_optimized object
#' @export
cbamm_optimized <- function(data,
                           mode = "fast",
                           config = cbamm_config(),
                           target_population = NULL,
                           effect_col = "yi",
                           se_col = "se",
                           study_col = "study_id",
                           ...) {
  performance_mode <- .normalize_cbamm_mode(mode)
  normalized_data <- .normalize_cbamm_public_data(
    data = data,
    effect_col = effect_col,
    se_col = se_col,
    study_col = study_col
  )

  result <- cbamm(
    data = normalized_data,
    config = config,
    target_population = target_population,
    performance_mode = performance_mode
  )

  result$diagnostics$requested_mode <- mode
  class(result) <- unique(c("cbamm_optimized", class(result)))
  result
}

#' @export
print.cbamm_optimized <- function(x, ...) {
  if (inherits(x, "cbamm_results")) {
    return(print.cbamm_results(x, ...))
  }

  cat("\nCBAMM Results (", x$mode, " mode)\n", sep = "")
  cat("========================\n")
  cat("Studies:", x$summary$n_studies, "\n")
  cat("Effect:", round(x$summary$estimate, 3), "\n")
  cat("95% CI: [", round(x$summary$ci_lower, 3), ", ", 
      round(x$summary$ci_upper, 3), "]\n", sep = "")
  cat("I²:", round(x$summary$I2, 1), "%\n")
  cat("Time:", round(x$execution_time, 3), "seconds\n")
  invisible(x)
}

