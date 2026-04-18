
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
                           mode = "quick",
                           effect_col = "yi",
                           se_col = "se", 
                           study_col = "study_id") {
  
  start_time <- Sys.time()
  
  if (!is.data.frame(data)) {
    stop("Data must be a data.frame")
  }
  
  if (!all(c(effect_col, se_col) %in% names(data))) {
    stop(paste("Data must contain", effect_col, "and", se_col, "columns"))
  }
  
  results <- list(data = data, mode = mode, timestamp = Sys.time())
  
  cat(sprintf("Running CBAMM analysis (%s mode)...\n", mode))
  cat("  Core meta-analysis...")
  
  weights <- 1 / data[[se_col]]^2
  weighted_effect <- sum(data[[effect_col]] * weights) / sum(weights)
  pooled_se <- sqrt(1 / sum(weights))
  Q <- sum(weights * (data[[effect_col]] - weighted_effect)^2)
  df <- nrow(data) - 1
  tau2 <- max(0, (Q - df) / (sum(weights) - sum(weights^2) / sum(weights)))
  I2 <- max(0, 100 * (Q - df) / Q)
  
  results$summary <- list(
    estimate = weighted_effect,
    se = pooled_se,
    ci_lower = weighted_effect - 1.96 * pooled_se,
    ci_upper = weighted_effect + 1.96 * pooled_se,
    tau2 = tau2,
    I2 = I2,
    Q = Q,
    n_studies = nrow(data)
  )
  
  cat(" done\n")
  
  if (mode == "standard" || mode == "full") {
    cat("  Sensitivity analysis...")
    loo_results <- lapply(1:min(nrow(data), 20), function(i) {
      data_loo <- data[-i, ]
      weights_loo <- 1 / data_loo[[se_col]]^2
      effect_loo <- sum(data_loo[[effect_col]] * weights_loo) / sum(weights_loo)
      list(excluded = i, estimate = effect_loo)
    })
    results$sensitivity <- loo_results
    cat(" done\n")
  }
  
  if (mode == "full") {
    cat("  Note: For full Bayesian, use original cbamm()\n")
  }
  
  elapsed <- difftime(Sys.time(), start_time, units = "secs")
  results$execution_time <- as.numeric(elapsed)
  class(results) <- c("cbamm_optimized", "list")
  cat(sprintf("\nCompleted in %.3f seconds\n", elapsed))
  
  return(results)
}

#' @export
print.cbamm_optimized <- function(x, ...) {
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

