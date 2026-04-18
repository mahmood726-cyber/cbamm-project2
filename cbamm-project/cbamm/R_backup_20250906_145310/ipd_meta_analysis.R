# ipd_meta_analysis.R
# Individual Participant Data Meta-Analysis Functions for CBAMM Package

#' Individual Participant Data Meta-Analysis
#' 
#' @param ipd_data List of individual participant data frames
#' @param aggregate_data Optional aggregate data frame
#' @param outcome Outcome variable name
#' @param treatment Treatment variable name
#' @param covariates Covariate names
#' @param method Analysis method
#' @param model_type Model type
#' @param clustering_var Clustering variable
#' @param interaction_terms Interaction terms
#' @param centering Whether to center covariates
#' @param ... Additional arguments
#' @return Object of class cbamm_ipd
#' @export
ipd_meta_analysis <- function(ipd_data,
                             aggregate_data = NULL,
                             outcome = "outcome",
                             treatment = "treatment",
                             covariates = NULL,
                             method = c("one-stage", "two-stage", "combined"),
                             model_type = c("fixed", "random", "multilevel"),
                             clustering_var = "study_id",
                             interaction_terms = NULL,
                             centering = TRUE,
                             ...) {
  
  method <- match.arg(method)
  model_type <- match.arg(model_type)
  
  if (!is.list(ipd_data)) {
    stop("ipd_data must be a list of data frames")
  }
  
  # Combine all IPD data
  all_data <- do.call(rbind, ipd_data)
  
  # Add study IDs if not present
  if (!clustering_var %in% names(all_data)) {
    study_ids <- rep(seq_along(ipd_data), sapply(ipd_data, nrow))
    all_data[[clustering_var]] <- study_ids
  }
  
  # Simple analysis (can be enhanced)
  formula_str <- paste(outcome, "~", treatment)
  if (!is.null(covariates)) {
    formula_str <- paste(formula_str, "+", paste(covariates, collapse = " + "))
  }
  
  # Fit model
  model <- lm(as.formula(formula_str), data = all_data)
  
  # Extract results
  coef_summary <- summary(model)$coefficients
  
  # Find treatment effect
  trt_row <- which(rownames(coef_summary) == treatment)
  if (length(trt_row) == 0) {
    trt_row <- 2  # Assume second row if not found by name
  }
  
  result <- list(
    estimate = coef_summary[trt_row, 1],
    se = coef_summary[trt_row, 2],
    p_value = coef_summary[trt_row, 4],
    ci_lower = coef_summary[trt_row, 1] - 1.96 * coef_summary[trt_row, 2],
    ci_upper = coef_summary[trt_row, 1] + 1.96 * coef_summary[trt_row, 2],
    model = model,
    data = all_data,
    method = method,
    model_type = model_type
  )
  
  class(result) <- c("cbamm_ipd", class(result))
  return(result)
}

#' Simulate IPD Data
#' 
#' @param n_studies Number of studies
#' @param n_patients Number of patients per study
#' @param treatment_effect True treatment effect
#' @param between_study_sd Between-study standard deviation
#' @return List of data frames
#' @export
simulate_ipd_data <- function(n_studies = 5, 
                             n_patients = 100,
                             treatment_effect = 0.5, 
                             between_study_sd = 0.2) {
  
  ipd_list <- list()
  
  for (i in 1:n_studies) {
    # Study-specific effect
    study_effect <- rnorm(1, treatment_effect, between_study_sd)
    
    # Generate patient data
    n <- n_patients
    treatment <- rbinom(n, 1, 0.5)
    age <- rnorm(n, 50, 10)
    
    # Generate outcome
    outcome <- study_effect * treatment + 0.02 * age + rnorm(n, 0, 1)
    
    ipd_list[[i]] <- data.frame(
      study_id = i,
      patient_id = 1:n,
      treatment = treatment,
      outcome = outcome,
      age = age
    )
  }
  
  return(ipd_list)
}

#' Print method for IPD results
#' @export
print.cbamm_ipd <- function(x, ...) {
  cat("IPD Meta-Analysis Results\n")
  cat("-------------------------\n")
  cat("Method:", x$method, "\n")
  cat("Model type:", x$model_type, "\n")
  cat("Effect:", round(x$estimate, 3), "\n")
  cat("SE:", round(x$se, 3), "\n")
  cat("95% CI: [", round(x$ci_lower, 3), ", ", round(x$ci_upper, 3), "]\n", sep = "")
  cat("p-value:", format.pval(x$p_value), "\n")
  invisible(x)
}

#' Summary method for IPD results
#' @export
summary.cbamm_ipd <- function(object, ...) {
  print(object)
}
