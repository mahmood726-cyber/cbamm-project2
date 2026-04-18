# ipd_meta_analysis.R
# Individual Participant Data Meta-Analysis Functions for CBAMM Package

#' Individual Participant Data Meta-Analysis
#' 
#' @param ipd_data List of individual participant data frames
#' @param aggregate_data Optional aggregate data frame
#' @param outcome Outcome variable name
#' @param treatment Treatment variable name
#' @param covariates Covariate names
#' @param method Analysis method: "one-stage" (default) or "two-stage"
#' @param model_type Model type: "fixed" or "random"
#' @param clustering_var Clustering variable (usually "study_id")
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
                             method = c("one-stage", "two-stage"),
                             model_type = c("random", "fixed"),
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
  
  # Ensure clustering variable is present
  if (!clustering_var %in% names(all_data)) {
    study_ids <- rep(seq_along(ipd_data), sapply(ipd_data, nrow))
    all_data[[clustering_var]] <- factor(study_ids)
  } else {
    all_data[[clustering_var]] <- factor(all_data[[clustering_var]])
  }
  
  if (method == "one-stage") {
    return(.run_one_stage_ipd(all_data, outcome, treatment, covariates, 
                             model_type, clustering_var, interaction_terms, ...))
  } else {
    return(.run_two_stage_ipd(ipd_data, outcome, treatment, covariates, 
                             model_type, ...))
  }
}

#' Run One-Stage IPD Meta-Analysis
#' @keywords internal
.run_one_stage_ipd <- function(all_data, outcome, treatment, covariates, 
                              model_type, clustering_var, interaction_terms, ...) {
  
  if (model_type == "random") {
    if (!check_package("lme4", quietly = TRUE)) {
      stop("lme4 package required for one-stage random-effects IPD")
    }
    
    # Formula: outcome ~ treatment + covariates + (treatment | study_id)
    # Random intercepts for studies and random slopes for treatment
    formula_str <- paste(outcome, "~", treatment)
    if (!is.null(covariates)) {
      formula_str <- paste(formula_str, "+", paste(covariates, collapse = " + "))
    }
    
    # Add random effects
    formula_str <- paste(formula_str, "+ (", treatment, "|", clustering_var, ")")
    
    # Add interactions with covariates if requested or present
    if (!is.null(covariates)) {
      interactions <- paste(paste0(treatment, ":", covariates), collapse = " + ")
      formula_str <- paste(formula_str, "+", interactions)
    }
    
    model <- lme4::lmer(as.formula(formula_str), data = all_data, ...)
    sum_mod <- summary(model)
    
    # Extract treatment effect
    coefs <- sum_mod$coefficients
    trt_idx <- which(rownames(coefs) == treatment)
    
    result <- list(
      estimate = coefs[trt_idx, "Estimate"],
      se = coefs[trt_idx, "Std. Error"],
      p_value = 2 * pnorm(abs(coefs[trt_idx, "t value"]), lower.tail = FALSE),
      model = model,
      sum_mod = sum_mod
    )
  } else {
    # Fixed effects (Stratified)
    formula_str <- paste(outcome, "~", treatment, "*", clustering_var)
    if (!is.null(covariates)) {
      formula_str <- paste(formula_str, "+", paste(covariates, collapse = " + "))
    }
    
    model <- lm(as.formula(formula_str), data = all_data)
    sum_mod <- summary(model)
    
    # This is a simplification; ideally use FE with study-specific intercepts
    result <- list(
      estimate = coef(model)[treatment],
      se = sqrt(vcov(model)[treatment, treatment]),
      p_value = sum_mod$coefficients[treatment, 4],
      model = model,
      sum_mod = sum_mod
    )
  }
  
  result$ci_lower <- result$estimate - 1.96 * result$se
  result$ci_upper <- result$estimate + 1.96 * result$se
  result$method <- "one-stage"
  result$model_type <- model_type
  
  class(result) <- c("cbamm_ipd", "list")
  return(result)
}

#' Run Two-Stage IPD Meta-Analysis
#' @keywords internal
.run_two_stage_ipd <- function(ipd_list, outcome, treatment, covariates, model_type, ...) {
  # 1. Estimate effect in each study
  stage1_results <- lapply(ipd_list, function(d) {
    fml <- paste(outcome, "~", treatment)
    if (!is.null(covariates)) fml <- paste(fml, "+", paste(covariates, collapse = " + "))
    
    fit <- lm(as.formula(fml), data = d)
    sum_fit <- summary(fit)$coefficients
    
    if (treatment %in% rownames(sum_fit)) {
      return(c(yi = sum_fit[treatment, 1], vi = sum_fit[treatment, 2]^2))
    }
    return(c(yi = NA, vi = NA))
  })
  
  stage1_df <- as.data.frame(do.call(rbind, stage1_results))
  stage1_df <- stage1_df[!is.na(stage1_df$yi), ]
  
  # 2. Meta-analyze stage 1 results
  if (check_package("metafor", quietly = TRUE)) {
    fit <- metafor::rma(yi = stage1_df$yi, vi = stage1_df$vi, 
                        method = ifelse(model_type == "random", "REML", "FE"))
    
    result <- list(
      estimate = as.numeric(fit$beta),
      se = fit$se,
      p_value = fit$pval,
      ci_lower = fit$ci.lb,
      ci_upper = fit$ci.ub,
      model = fit,
      method = "two-stage",
      model_type = model_type
    )
  } else {
    # Simple weighted average fallback
    w <- 1 / stage1_df$vi
    est <- sum(w * stage1_df$yi) / sum(w)
    se <- sqrt(1 / sum(w))
    
    result <- list(
      estimate = est,
      se = se,
      p_value = 2 * pnorm(abs(est/se), lower.tail = FALSE),
      ci_lower = est - 1.96*se,
      ci_upper = est + 1.96*se,
      method = "two-stage",
      model_type = model_type
    )
  }
  
  class(result) <- c("cbamm_ipd", "list")
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
      study_id = factor(rep(i, n)),
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
