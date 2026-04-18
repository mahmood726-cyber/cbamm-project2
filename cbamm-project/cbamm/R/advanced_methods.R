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

#' Bias-Adjusted Cross-Design Synthesis
#'
#' Integrates pooling, publication bias adjustment, and heterogeneity
#' in a single cross-design synthesis.
#'
#' @param rct_data Data frame with RCT study-level data (effect, variance)
#' @param obs_data Data frame with observational study-level data
#' @param bias_offset Prior belief about observational bias (default 0)
#' @param bias_inflation Scaling factor for observational variance (default 1)
#' @param asymmetry_p_threshold P-value threshold used to flag funnel asymmetry
#' @param use_pet_peese Whether to apply PET-PEESE style correction when
#'   asymmetry is detected
#'
#' @return A list containing RCT, observational, and combined summaries with
#'   publication bias diagnostics and heterogeneity estimates
#' @export
cross_design_synthesis <- function(rct_data,
                                   obs_data,
                                   bias_offset = 0,
                                   bias_inflation = 1,
                                   asymmetry_p_threshold = 0.10,
                                   use_pet_peese = TRUE) {

  .validate_cross_design_data <- function(d, label) {
    if (is.null(d) || nrow(d) == 0) {
      return(data.frame(effect = numeric(0), variance = numeric(0)))
    }

    if (!all(c("effect", "variance") %in% names(d))) {
      stop(label, " data must contain 'effect' and 'variance' columns")
    }

    d <- d[is.finite(d$effect) & is.finite(d$variance) & d$variance > 0, , drop = FALSE]

    if (nrow(d) == 0) {
      return(data.frame(effect = numeric(0), variance = numeric(0)))
    }

    d
  }

  .random_effects_pool <- function(effect, variance) {
    k <- length(effect)
    if (k == 0) {
      return(list(
        effect = NA_real_,
        se = NA_real_,
        ci_lower = NA_real_,
        ci_upper = NA_real_,
        tau2 = NA_real_,
        I2 = NA_real_,
        Q = NA_real_,
        Q_p = NA_real_,
        k = 0L,
        weights = numeric(0)
      ))
    }

    w_fe <- 1 / variance
    mu_fe <- sum(w_fe * effect) / sum(w_fe)

    if (k == 1) {
      se_fe <- sqrt(variance[1])
      return(list(
        effect = effect[1],
        se = se_fe,
        ci_lower = effect[1] - 1.96 * se_fe,
        ci_upper = effect[1] + 1.96 * se_fe,
        tau2 = 0,
        I2 = 0,
        Q = 0,
        Q_p = 1,
        k = 1L,
        weights = 1
      ))
    }

    Q <- sum(w_fe * (effect - mu_fe)^2)
    df <- k - 1
    c_val <- sum(w_fe) - sum(w_fe^2) / sum(w_fe)
    tau2 <- if (is.finite(c_val) && c_val > 0) max(0, (Q - df) / c_val) else 0

    w_re <- 1 / (variance + tau2)
    mu_re <- sum(w_re * effect) / sum(w_re)
    se_re <- sqrt(1 / sum(w_re))
    i2 <- if (is.finite(Q) && Q > 0 && Q > df) 100 * (Q - df) / Q else 0

    list(
      effect = mu_re,
      se = se_re,
      ci_lower = mu_re - 1.96 * se_re,
      ci_upper = mu_re + 1.96 * se_re,
      tau2 = tau2,
      I2 = i2,
      Q = Q,
      Q_p = stats::pchisq(Q, df = df, lower.tail = FALSE),
      k = k,
      weights = w_re
    )
  }

  .publication_bias_diagnostics <- function(effect, variance, threshold) {
    k <- length(effect)
    out <- list(
      k = k,
      egger_intercept = NA_real_,
      egger_p = NA_real_,
      pet_intercept = NA_real_,
      pet_se = NA_real_,
      pet_p = NA_real_,
      peese_intercept = NA_real_,
      peese_se = NA_real_,
      chosen_adjustment = "none",
      adjusted_effect = NA_real_,
      adjusted_se = NA_real_,
      asymmetry_detected = FALSE,
      mean_p_value = NA_real_
    )

    if (k < 3) {
      return(out)
    }

    sei <- sqrt(variance)
    precision <- 1 / sei
    z_scores <- effect / sei
    p_values <- 2 * stats::pnorm(-abs(z_scores))
    out$mean_p_value <- mean(p_values, na.rm = TRUE)

    if (stats::sd(precision) == 0) {
      return(out)
    }

    w <- 1 / variance
    egger_fit <- try(stats::lm(z_scores ~ precision, weights = w), silent = TRUE)
    pet_fit <- try(stats::lm(effect ~ sei, weights = w), silent = TRUE)
    peese_fit <- try(stats::lm(effect ~ I(sei^2), weights = w), silent = TRUE)

    if (!inherits(egger_fit, "try-error")) {
      egger_coef <- summary(egger_fit)$coefficients
      out$egger_intercept <- unname(egger_coef[1, 1])
      out$egger_p <- unname(egger_coef[1, 4])
      out$asymmetry_detected <- is.finite(out$egger_p) && out$egger_p < threshold
    }

    if (!inherits(pet_fit, "try-error")) {
      pet_coef <- summary(pet_fit)$coefficients
      out$pet_intercept <- unname(pet_coef[1, 1])
      out$pet_se <- unname(pet_coef[1, 2])
      out$pet_p <- unname(pet_coef[1, 4])
    }

    if (!inherits(peese_fit, "try-error")) {
      peese_coef <- summary(peese_fit)$coefficients
      out$peese_intercept <- unname(peese_coef[1, 1])
      out$peese_se <- unname(peese_coef[1, 2])
    }

    if (isTRUE(out$asymmetry_detected) && is.finite(out$peese_intercept)) {
      out$adjusted_effect <- out$peese_intercept
      out$adjusted_se <- out$peese_se
      out$chosen_adjustment <- "PEESE"
    } else if (isTRUE(out$asymmetry_detected) && is.finite(out$pet_intercept)) {
      out$adjusted_effect <- out$pet_intercept
      out$adjusted_se <- out$pet_se
      out$chosen_adjustment <- "PET"
    }

    out
  }

  .summarize_design <- function(d,
                                label,
                                bias_offset = 0,
                                bias_inflation = 1,
                                threshold = 0.10,
                                use_pet_peese = TRUE) {
    d <- .validate_cross_design_data(d, label)

    if (nrow(d) == 0) {
      return(list(
        effect = NA_real_,
        raw_effect = NA_real_,
        bias_adjusted_effect = NA_real_,
        se = NA_real_,
        ci_lower = NA_real_,
        ci_upper = NA_real_,
        tau2 = NA_real_,
        I2 = NA_real_,
        Q = NA_real_,
        Q_p = NA_real_,
        k = 0L,
        publication_bias = .publication_bias_diagnostics(numeric(0), numeric(0), threshold),
        design = label
      ))
    }

    working <- d
    working$effect_adj <- working$effect - bias_offset
    working$variance_adj <- working$variance * bias_inflation

    pooled <- .random_effects_pool(working$effect_adj, working$variance_adj)
    pub_bias <- .publication_bias_diagnostics(working$effect_adj, working$variance_adj, threshold)

    chosen_effect <- pooled$effect
    chosen_se <- pooled$se

    adjustment_applied <- isTRUE(use_pet_peese) &&
        isTRUE(pub_bias$asymmetry_detected) &&
        is.finite(pub_bias$adjusted_effect)

    if (adjustment_applied) {
      chosen_effect <- pub_bias$adjusted_effect
      if (is.finite(pub_bias$adjusted_se)) {
        chosen_se <- max(chosen_se, pub_bias$adjusted_se, na.rm = TRUE)
      }
    }

    list(
      effect = chosen_effect,
      raw_effect = pooled$effect,
      bias_adjusted_effect = if (adjustment_applied) pub_bias$adjusted_effect else NA_real_,
      se = chosen_se,
      ci_lower = chosen_effect - 1.96 * chosen_se,
      ci_upper = chosen_effect + 1.96 * chosen_se,
      tau2 = pooled$tau2,
      I2 = pooled$I2,
      Q = pooled$Q,
      Q_p = pooled$Q_p,
      k = pooled$k,
      publication_bias = pub_bias,
      design = label
    )
  }

  rct_res <- .summarize_design(
    d = rct_data,
    label = "RCT",
    threshold = asymmetry_p_threshold,
    use_pet_peese = use_pet_peese
  )

  obs_res <- .summarize_design(
    d = obs_data,
    label = "OBS",
    bias_offset = bias_offset,
    bias_inflation = bias_inflation,
    threshold = asymmetry_p_threshold,
    use_pet_peese = use_pet_peese
  )

  combined_designs <- data.frame(
    design = c("RCT", "OBS"),
    effect = c(rct_res$effect, obs_res$effect),
    variance = c(rct_res$se^2, obs_res$se^2)
  )
  combined_designs <- combined_designs[
    is.finite(combined_designs$effect) &
      is.finite(combined_designs$variance) &
      combined_designs$variance > 0,
    ,
    drop = FALSE
  ]

  combined_res <- .random_effects_pool(combined_designs$effect, combined_designs$variance)
  design_weights <- if (length(combined_res$weights) == nrow(combined_designs)) {
    stats::setNames(as.numeric(combined_res$weights / sum(combined_res$weights)), combined_designs$design)
  } else {
    numeric(0)
  }

  all_effects <- numeric(0)
  all_variances <- numeric(0)
  if (!is.null(rct_data) && nrow(rct_data) > 0) {
    all_effects <- c(all_effects, rct_data$effect)
    all_variances <- c(all_variances, rct_data$variance)
  }
  if (!is.null(obs_data) && nrow(obs_data) > 0) {
    all_effects <- c(all_effects, obs_data$effect - bias_offset)
    all_variances <- c(all_variances, obs_data$variance * bias_inflation)
  }

  list(
    rct = rct_res,
    obs = obs_res,
    combined = list(
      effect = combined_res$effect,
      se = combined_res$se,
      ci_lower = combined_res$ci_lower,
      ci_upper = combined_res$ci_upper,
      tau2 = combined_res$tau2,
      I2 = combined_res$I2,
      Q = combined_res$Q,
      Q_p = combined_res$Q_p,
      k = combined_res$k,
      design_weights = design_weights
    ),
    heterogeneity = list(
      within_design = c(RCT = rct_res$tau2, OBS = obs_res$tau2),
      between_design = combined_res$tau2
    ),
    publication_bias = list(
      rct = rct_res$publication_bias,
      obs = obs_res$publication_bias,
      combined = .publication_bias_diagnostics(all_effects, all_variances, asymmetry_p_threshold)
    ),
    bias_parameters = list(
      offset = bias_offset,
      inflation = bias_inflation,
      asymmetry_p_threshold = asymmetry_p_threshold,
      use_pet_peese = use_pet_peese
    )
  )
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
