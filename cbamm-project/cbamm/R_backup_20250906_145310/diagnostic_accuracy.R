# diagnostic_accuracy.R
# Diagnostic Test Accuracy Meta-Analysis Functions for CBAMM Package

#' Diagnostic Test Accuracy Meta-Analysis
#' 
#' @param data Data frame containing diagnostic test data
#' @param TP Column name for true positives
#' @param FP Column name for false positives
#' @param FN Column name for false negatives
#' @param TN Column name for true negatives
#' @param study_id Column name for study identifier
#' @param method Analysis method
#' @param add Continuity correction value
#' @param ... Additional arguments
#' @return Object of class cbamm_dta
#' @export
diagnostic_accuracy <- function(data, 
                               TP = "TP", 
                               FP = "FP", 
                               FN = "FN", 
                               TN = "TN",
                               study_id = "study_id",
                               method = c("bivariate", "hsroc", "univariate"),
                               add = 0.5, 
                               ...) {
  
  method <- match.arg(method)
  
  # Extract columns
  tp <- data[[TP]]
  fp <- data[[FP]]
  fn <- data[[FN]]
  tn <- data[[TN]]
  
  if (!is.null(study_id) && study_id %in% names(data)) {
    study <- data[[study_id]]
  } else {
    study <- seq_len(nrow(data))
  }
  
  # Apply continuity correction if needed
  zero_cells <- (tp == 0) | (fp == 0) | (fn == 0) | (tn == 0)
  if (any(zero_cells) && add > 0) {
    tp[zero_cells] <- tp[zero_cells] + add
    fp[zero_cells] <- fp[zero_cells] + add
    fn[zero_cells] <- fn[zero_cells] + add
    tn[zero_cells] <- tn[zero_cells] + add
  }
  
  # Calculate sensitivity and specificity
  sensitivity <- tp / (tp + fn)
  specificity <- tn / (tn + fp)
  
  # Simple pooling (can be enhanced with bivariate models)
  sens_pooled <- mean(sensitivity, na.rm = TRUE)
  spec_pooled <- mean(specificity, na.rm = TRUE)
  
  # Calculate confidence intervals (simple percentile method)
  sens_ci <- quantile(sensitivity, c(0.025, 0.975), na.rm = TRUE)
  spec_ci <- quantile(specificity, c(0.025, 0.975), na.rm = TRUE)
  
  # Create result object
  result <- list(
    summary = list(
      sensitivity = sens_pooled,
      specificity = spec_pooled,
      ci_sensitivity = sens_ci,
      ci_specificity = spec_ci
    ),
    studies = data.frame(
      study_id = study,
      sensitivity = sensitivity,
      specificity = specificity,
      tp = tp,
      fp = fp,
      fn = fn,
      tn = tn
    ),
    method = method,
    data = data
  )
  
  class(result) <- c("cbamm_dta", class(result))
  return(result)
}

#' @export
summary.cbamm_dta <- function(object, ...) {
  cat("\nDiagnostic Test Accuracy Meta-Analysis\n")
  cat("=======================================\n")
  cat("Method:", object$method, "\n")
  cat("Studies:", length(unique(object$studies$study_id)), "\n\n")
  cat("Pooled Estimates:\n")
  cat("  Sensitivity:", round(object$summary$sensitivity, 3), "\n")
  cat("  Specificity:", round(object$summary$specificity, 3), "\n")
  cat("  95% CI Sensitivity: [", round(object$summary$ci_sensitivity[1], 3), 
      ", ", round(object$summary$ci_sensitivity[2], 3), "]\n", sep = "")
  cat("  95% CI Specificity: [", round(object$summary$ci_specificity[1], 3),
      ", ", round(object$summary$ci_specificity[2], 3), "]\n", sep = "")
  invisible(object)
}

#' Simulate Diagnostic Test Accuracy Data
#' 
#' @param n_studies Number of studies to simulate
#' @param sens_mean Mean sensitivity
#' @param spec_mean Mean specificity
#' @param n_diseased Number of diseased patients per study
#' @param n_healthy Number of healthy patients per study
#' @return Data frame with diagnostic test data
#' @export
simulate_dta_data <- function(n_studies = 10, 
                             sens_mean = 0.8, 
                             spec_mean = 0.9,
                             n_diseased = 100,
                             n_healthy = 100) {
  
  # Simulate true sensitivities and specificities with some variation
  sens_true <- pmin(pmax(rnorm(n_studies, sens_mean, 0.1), 0.1), 0.99)
  spec_true <- pmin(pmax(rnorm(n_studies, spec_mean, 0.05), 0.1), 0.99)
  
  # Generate 2x2 table data
  TP <- rbinom(n_studies, n_diseased, sens_true)
  FN <- n_diseased - TP
  TN <- rbinom(n_studies, n_healthy, spec_true)
  FP <- n_healthy - TN
  
  return(data.frame(
    study_id = 1:n_studies,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN
  ))
}

#' Print method for diagnostic accuracy results
#' @export
print.cbamm_dta <- function(x, ...) {
  summary(x)
}
