# advanced_methods.R
# Advanced Meta-Analysis Methods for CBAMM Package

#' Calculate Prediction Intervals
#' 
#' @param effect Effect size estimate
#' @param se Standard error
#' @param tau2 Between-study variance (heterogeneity)
#' @param level Confidence level (default 0.95)
#' @export
calculate_prediction_interval <- function(effect, se, tau2 = 0, level = 0.95) {
  pred_se <- sqrt(se^2 + tau2)
  z_crit <- qnorm((1 + level) / 2)
  pi_lower <- effect - z_crit * pred_se
  pi_upper <- effect + z_crit * pred_se
  
  return(list(
    prediction_interval = c(lower = pi_lower, upper = pi_upper),
    level = level,
    prediction_se = pred_se,
    tau2 = tau2
  ))
}

#' Simulate IPD Data
#' @export
simulate_ipd_data <- function(n_studies = 5, n_patients = 100, effect = 0.5) {
  ipd_list <- list()
  for (i in 1:n_studies) {
    n <- n_patients
    treatment <- rbinom(n, 1, 0.5)
    outcome <- effect * treatment + rnorm(n)
    ipd_list[[i]] <- data.frame(
      study_id = i,
      patient_id = 1:n,
      treatment = treatment,
      outcome = outcome
    )
  }
  return(ipd_list)
}

#' Simple IPD Meta-Analysis
#' @export
ipd_meta_analysis <- function(ipd_data, outcome = "outcome", treatment = "treatment") {
  # Combine all data
  all_data <- do.call(rbind, ipd_data)
  
  # Simple fixed effects model
  model <- lm(as.formula(paste(outcome, "~", treatment)), data = all_data)
  
  # Extract results
  coef_summary <- summary(model)$coefficients
  trt_row <- which(rownames(coef_summary) == treatment)
  
  result <- list(
    estimate = coef_summary[trt_row, 1],
    se = coef_summary[trt_row, 2],
    p_value = coef_summary[trt_row, 4],
    ci_lower = coef_summary[trt_row, 1] - 1.96 * coef_summary[trt_row, 2],
    ci_upper = coef_summary[trt_row, 1] + 1.96 * coef_summary[trt_row, 2],
    model = model,
    data = all_data
  )
  
  class(result) <- "cbamm_ipd"
  return(result)
}

#' Simulate Diagnostic Test Data
#' @export
simulate_dta_data <- function(n_studies = 10, sens = 0.8, spec = 0.9) {
  TP <- rbinom(n_studies, 100, sens)
  FN <- 100 - TP
  TN <- rbinom(n_studies, 100, spec)
  FP <- 100 - TN
  
  return(data.frame(
    study_id = 1:n_studies,
    TP = TP, FP = FP, FN = FN, TN = TN
  ))
}

#' Simple Diagnostic Accuracy Analysis
#' @export
diagnostic_accuracy <- function(data, TP = "TP", FP = "FP", FN = "FN", TN = "TN") {
  tp <- data[[TP]]
  fp <- data[[FP]]
  fn <- data[[FN]]
  tn <- data[[TN]]
  
  # Calculate sensitivity and specificity
  sensitivity <- tp / (tp + fn)
  specificity <- tn / (tn + fp)
  
  result <- list(
    sensitivity = mean(sensitivity),
    specificity = mean(specificity),
    sens_ci = quantile(sensitivity, c(0.025, 0.975)),
    spec_ci = quantile(specificity, c(0.025, 0.975)),
    studies = data.frame(
      sensitivity = sensitivity,
      specificity = specificity
    )
  )
  
  class(result) <- "cbamm_dta"
  return(result)
}

#' Simulate Sequential Data
#' @export
simulate_sequential_data <- function(n_studies = 10, true_effect = 0.3) {
  effects <- rnorm(n_studies, true_effect, 0.1)
  ses <- runif(n_studies, 0.05, 0.2)
  
  return(data.frame(
    study_id = 1:n_studies,
    effect = effects,
    se = ses,
    variance = ses^2
  ))
}

#' Simple Living Review
#' @export
living_systematic_review <- function(data, effect = "effect", se = "se", alpha = 0.05) {
  effects <- data[[effect]]
  ses <- data[[se]]
  n_studies <- nrow(data)
  
  # Calculate cumulative effects
  cumulative_effects <- numeric(n_studies)
  cumulative_z <- numeric(n_studies)
  
  for (i in 1:n_studies) {
    weights <- 1 / ses[1:i]^2
    cumulative_effects[i] <- sum(effects[1:i] * weights) / sum(weights)
    cumulative_se <- sqrt(1 / sum(weights))
    cumulative_z[i] <- cumulative_effects[i] / cumulative_se
  }
  
  # Simple boundaries (O'Brien-Fleming approximation)
  boundaries <- qnorm(1 - alpha/2) / sqrt((1:n_studies) / n_studies)
  
  result <- list(
    cumulative_effects = cumulative_effects,
    cumulative_z = cumulative_z,
    boundaries = boundaries,
    stopped = any(abs(cumulative_z) > boundaries),
    data = data
  )
  
  class(result) <- "cbamm_living"
  return(result)
}

#' Cross-Design Synthesis
#' @export
cross_design_synthesis <- function(rct_data, obs_data) {
  # Pool RCT data
  rct_weights <- 1 / rct_data$variance
  rct_effect <- sum(rct_data$effect * rct_weights) / sum(rct_weights)
  rct_se <- sqrt(1 / sum(rct_weights))
  
  # Pool observational data
  obs_weights <- 1 / obs_data$variance
  obs_effect <- sum(obs_data$effect * obs_weights) / sum(obs_weights)
  obs_se <- sqrt(1 / sum(obs_weights))
  
  # Estimate bias
  bias <- obs_effect - rct_effect
  
  # Combined estimate (simple average weighted by precision)
  combined_effect <- (rct_effect / rct_se^2 + obs_effect / obs_se^2) / 
                    (1 / rct_se^2 + 1 / obs_se^2)
  combined_se <- sqrt(1 / (1 / rct_se^2 + 1 / obs_se^2))
  
  return(list(
    rct = list(effect = rct_effect, se = rct_se),
    obs = list(effect = obs_effect, se = obs_se),
    combined = list(effect = combined_effect, se = combined_se),
    bias = bias
  ))
}

#' Print method for IPD results
#' @export
print.cbamm_ipd <- function(x, ...) {
  cat("IPD Meta-Analysis Results
")
  cat("-------------------------
")
  cat("Effect:", round(x$estimate, 3), "
")
  cat("95% CI: [", round(x$ci_lower, 3), ", ", round(x$ci_upper, 3), "]
", sep = "")
  cat("p-value:", format.pval(x$p_value), "
")
  invisible(x)
}

#' Print method for DTA results
#' @export
print.cbamm_dta <- function(x, ...) {
  cat("Diagnostic Test Accuracy
")
  cat("------------------------
")
  cat("Sensitivity:", round(x$sensitivity, 3), "
")
  cat("Specificity:", round(x$specificity, 3), "
")
  invisible(x)
}

#' Print method for Living Review
#' @export
print.cbamm_living <- function(x, ...) {
  cat("Living Systematic Review
")
  cat("-----------------------
")
  cat("Studies:", nrow(x$data), "
")
  cat("Stopped:", x$stopped, "
")
  cat("Current Z:", round(tail(x$cumulative_z, 1), 3), "
")
  invisible(x)
}

