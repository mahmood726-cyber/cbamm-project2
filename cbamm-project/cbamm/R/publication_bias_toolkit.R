# publication_bias_toolkit.R
# Advanced Publication Bias Detection for CBAMM

#' P-Curve Analysis
#' 
#' @param p_values Vector of p-values
#' @export
p_curve_analysis <- function(p_values) {
  
  # Only use significant p-values
  sig_p <- p_values[p_values < 0.05]
  
  if (length(sig_p) < 2) {
    return(list(
      evidential_value = NA,
      interpretation = "Too few significant results"
    ))
  }
  
  # Test for right-skewness
  n_low <- sum(sig_p < 0.025)
  n_total <- length(sig_p)
  binom_test <- binom.test(n_low, n_total, p = 0.5)
  
  # Evidential value test
  evidential_value <- binom_test$p.value < 0.05 && n_low > n_total/2
  
  return(list(
    evidential_value = evidential_value,
    p_value = binom_test$p.value,
    n_significant = n_total,
    interpretation = ifelse(evidential_value, 
                           "Evidential value present", 
                           "Evidential value absent or inconclusive")
  ))
}

#' Publication Bias Suite
#' 
#' @param ma_result Meta-analysis result
#' @export
publication_bias_suite <- function(ma_result) {
  
  results <- list()
  
  # Egger's test
  if (!is.null(ma_result$data)) {
    effects <- ma_result$data$effect
    ses <- sqrt(ma_result$data$variance)
    precision <- 1/ses
    z_scores <- effects/ses
    
    egger_model <- lm(z_scores ~ precision)
    egger_p <- summary(egger_model)$coefficients[1, 4]
    
    results$egger <- list(
      p_value = egger_p,
      significant = egger_p < 0.05,
      interpretation = ifelse(egger_p < 0.05, "Asymmetry detected", "No asymmetry")
    )
  }
  
  # P-curve if p-values available
  if (!is.null(ma_result$data$p_value)) {
    results$pcurve <- p_curve_analysis(ma_result$data$p_value)
  }
  
  # Overall assessment
  n_tests <- length(results)
  n_significant <- sum(sapply(results, function(x) x$significant), na.rm = TRUE)
  
  results$overall <- list(
    bias_level = ifelse(n_significant == 0, "Low",
                        ifelse(n_significant < n_tests/2, "Moderate", "High")),
    recommendation = ifelse(n_significant == 0, "Results appear robust",
                           ifelse(n_significant < n_tests/2, 
                                  "Interpret with caution",
                                  "Consider bias-adjusted estimates"))
  )
  
  class(results) <- "cbamm_bias_suite"
  return(results)
}

#' Print method for bias suite
#' @export
print.cbamm_bias_suite <- function(x, ...) {
  cat("\nPublication Bias Assessment\n")
  cat("===========================\n\n")
  
  if (!is.null(x$egger)) {
    cat("Egger's Test:\n")
    cat("  ", x$egger$interpretation, "\n")
    cat("  p-value: ", format.pval(x$egger$p_value), "\n\n")
  }
  
  if (!is.null(x$pcurve)) {
    cat("P-Curve:\n")
    cat("  ", x$pcurve$interpretation, "\n\n")
  }
  
  cat("Overall Assessment:\n")
  cat("  Bias level: ", x$overall$bias_level, "\n")
  cat("  ", x$overall$recommendation, "\n")
  
  invisible(x)
}

