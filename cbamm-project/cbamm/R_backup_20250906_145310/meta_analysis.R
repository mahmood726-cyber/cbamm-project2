#' Robust Meta-Analysis with HKSJ Adjustment
#'
#' @param yi Effect sizes
#' @param sei Standard errors
#' @param data Optional data frame
#' @param method Variance estimator method
#' @param weights Optional study weights
#' @param mods Optional moderators
#' @param use_hksj Whether to use Hartung-Knapp-Sidik-Jonkman adjustment
#'
#' @return Meta-analysis results object
robust_rma <- function(yi, sei, data = NULL, method = "REML",
                       weights = NULL, mods = NULL, use_hksj = TRUE) {
  stopifnot(length(yi) == length(sei))
  
  if (any(!is.finite(yi)) || any(!is.finite(sei)) || any(sei <= 0)) {
    stop("Non-finite or non-positive inputs in yi/sei")
  }
  
  if (length(yi) < 3) stop("Need at least 3 studies for meta-analysis")
  
  build_args <- function(use_vi = FALSE) {
    a <- list(yi = yi, method = method)
    if (use_vi) a$vi <- sei^2 else a$sei <- sei
    if (!is.null(data))    a$data    <- data
    if (!is.null(weights)) a$weights <- weights
    if (!is.null(mods))    a$mods    <- mods
    if (!is.null(use_hksj) && isTRUE(use_hksj)) a$test <- "knha"
    a
  }
  
  # Try different combinations if first fails
  args <- build_args(FALSE)
  fit <- try(do.call(metafor::rma, args), silent = TRUE)
  
  if (inherits(fit, "try-error")) {
    args$test <- NULL
    fit <- try(do.call(metafor::rma, args), silent = TRUE)
  }
  
  if (inherits(fit, "try-error")) stop("Meta-analysis failed: ", as.character(fit))
  fit
}

#' Run Core Meta-Analysis
#'
#' @param data Study data
#' @param config CBAMM configuration
#'
#' @return List with meta-analysis results
#' @keywords internal
run_core_meta_analysis <- function(data, config) {
  results <- list()
  
  # Basic pooled analysis with transport weights
  if (config$output$verbose) message("  Computing transport-weighted analysis...")
  
  fit_transport <- robust_rma(
    yi = data$yi, 
    sei = data$se, 
    data = data,
    method = config$estimators[1],  # Use first estimator as default
    weights = data$analysis_weights,
    use_hksj = config$methods$hksj
  )
  
  results$transport <- fit_transport
  
  # GRADE-weighted analysis if available
  if ("analysis_weights_grade" %in% names(data)) {
    if (config$output$verbose) message("  Computing GRADE-weighted analysis...")
    
    fit_grade <- robust_rma(
      yi = data$yi, 
      sei = data$se, 
      data = data,
      method = config$estimators[1],
      weights = data$analysis_weights_grade,
      use_hksj = config$methods$hksj
    )
    
    results$grade <- fit_grade
  }
  
  # Robust variance estimation if duplicates present
  if (isTRUE(config$methods$robust_variance) && any(duplicated(data$study_id))) {
    if (config$output$verbose) message("  Computing robust variance estimates...")
    results$rve <- .run_rve_analysis(fit_transport, data)
  }
  
  return(results)
}

#' Run Robust Variance Estimation
#' @keywords internal
.run_rve_analysis <- function(fit, data) {
  if (!check_package("clubSandwich", quietly = TRUE)) {
    warning("clubSandwich not available for RVE analysis")
    return(NULL)
  }
  
  with_fallback(
    pkg = "clubSandwich",
    fun_call = {
      vc <- clubSandwich::vcovCR(fit, cluster = data$study_id, type = "CR2")
      rob <- clubSandwich::coef_test(fit, vcov = vc)
      list(vcov = vc, coef_test = rob)
    },
    fallback_fun = function() NULL,
    fallback_msg = "RVE analysis failed"
  )
}

#' Report Meta-Analysis Results
#'
#' @param result Meta-analysis result object
#' @param label Label for output
#' @param include_pi Whether to include prediction intervals
#'
#' @return Invisible NULL (prints results)
report_meta_result <- function(result, label = "", include_pi = TRUE) {
  if (inherits(result, "try-error")) { 
    message(label, ": failed")
    return(invisible(NULL)) 
  }
  
  pr <- try(metafor::predict(result, transf = exp), silent = TRUE)
  if (!inherits(pr, "try-error")) {
    have_pi <- include_pi && !any(is.na(c(pr$pi.lb, pr$pi.ub)))
    if (have_pi) {
      message(sprintf("%s: HR=%.3f (95%% CI %.3f–%.3f, PI %.3f–%.3f) | k=%d, τ²=%.4f, I²=%.1f%%",
                  label, pr$pred, pr$ci.lb, pr$ci.ub, pr$pi.lb, pr$pi.ub, result$k, result$tau2, result$I2))
    } else {
      message(sprintf("%s: HR=%.3f (95%% CI %.3f–%.3f) | k=%d, τ²=%.4f, I²=%.1f%%",
                  label, pr$pred, pr$ci.lb, pr$ci.ub, result$k, result$tau2, result$I2))
    }
  } else {
    message(sprintf("%s: HR=%.3f (CI unavailable) | k=%d, τ²=%.4f, I²=%.1f%%",
                label, exp(coef(result)), result$k, result$tau2, result$I2))
  }
}

#' Run Stratified Meta-Analysis
#'
#' @param data Study data
#' @param config CBAMM configuration
#'
#' @return List of stratified results
#' @keywords internal
run_stratified_analysis <- function(data, config) {
  if (config$output$verbose) message("Running stratified analysis by study type...")
  
  res <- list()
  for (type in unique(as.character(data$study_type))) {
    d <- dplyr::filter(data, study_type == type)
    if (nrow(d) < 3) {
      if (config$output$verbose) message(sprintf("%-8s: insufficient studies (n=%d)", type, nrow(d)))
      res[[type]] <- NULL
      next
    }
    
    w <- d$analysis_weights
    fit <- try(robust_rma(d$yi, d$se, data = d, method = config$estimators[1],
                          weights = w, use_hksj = config$methods$hksj), silent = TRUE)
    
    if (!inherits(fit, "try-error")) {
      if (config$output$verbose) report_meta_result(fit, sprintf("%-8s", type))
      res[[type]] <- fit
    } else {
      if (config$output$verbose) message(sprintf("%-8s: analysis failed", type))
      res[[type]] <- NULL
    }
  }
  
  res
}

#' Run Transport Analysis
#'
#' @param data Study data with transport weights
#' @param target_population Target population characteristics
#' @param config CBAMM configuration
#'
#' @return Transport analysis results
#' @keywords internal
run_transport_analysis <- function(data, target_population, config) {
  
  # Compute transport weights if not already present
  if (!"transport_weights" %in% names(data)) {
    if (is.null(target_population)) {
      data$transport_weights <- rep(1 / nrow(data), nrow(data))
    } else {
      data$transport_weights <- compute_transport_weights(
        data, target_population, config$transportability$truncation %||% 0.02
      )
    }
  }
  
  # Compute analysis weights
  data$analysis_weights <- compute_analysis_weights(data, data$transport_weights)
  
  # Apply GRADE weighting if enabled and available
  if (isTRUE(config$methods$grade_weighting) && "grade" %in% names(data)) {
    data$analysis_weights_grade <- apply_grade_weighting(data$analysis_weights, data$grade)
  }
  
  # Run the analysis
  results <- list()
  
  # Transport-only analysis
  fit_transport <- robust_rma(
    yi = data$yi, 
    sei = data$se,
    weights = data$analysis_weights,
    method = config$estimators[1],
    use_hksj = config$methods$hksj
  )
  results$transport <- fit_transport
  
  # Transport + GRADE analysis
  if ("analysis_weights_grade" %in% names(data)) {
    fit_grade <- robust_rma(
      yi = data$yi, 
      sei = data$se,
      weights = data$analysis_weights_grade,
      method = config$estimators[1],
      use_hksj = config$methods$hksj
    )
    results$grade <- fit_grade
  }
  
  # Report balance if target population provided
  if (!is.null(target_population) && all(c("age_mean", "female_pct") %in% names(data))) {
    pre_age  <- mean(data$age_mean, na.rm = TRUE)
    post_age <- stats::weighted.mean(data$age_mean, data$transport_weights, na.rm = TRUE)
    
    if (config$output$verbose) {
      message(sprintf("Age balance: %.1f -> %.1f (target %.1f)", 
                     pre_age, post_age, target_population$age_mean))
    }
    
    results$balance <- list(
      pre_age = pre_age,
      post_age = post_age,
      target_age = target_population$age_mean
    )
  }
  
  return(results)
}


#' Safe Meta-Analysis Wrapper
#' @keywords internal
safe_robust_rma <- function(yi, sei, data = NULL, method = "REML",
                           weights = NULL, mods = NULL, use_hksj = TRUE) {
  
  # Input validation
  if (length(yi) == 0 || length(sei) == 0) {
    stop("Empty effect sizes or standard errors")
  }
  
  if (is.null(method) || length(method) == 0) {
    method <- "REML"  # Safe default
  }
  
  # Ensure method is valid
  valid_methods <- c("REML", "DL", "PM", "ML", "EB", "SJ")
  if (!method[1] %in% valid_methods) {
    warning("Invalid method '", method[1], "', using REML")
    method <- "REML"
  }
  
  # Call original robust_rma with safe parameters
  robust_rma(yi = yi, sei = sei, data = data, method = method[1],
            weights = weights, mods = mods, use_hksj = use_hksj)
}

#' Safe Meta-Analysis Runner
#' @keywords internal  
safe_run_core_meta_analysis <- function(data, config) {
  
  results <- list()
  
  # Validate inputs
  if (is.null(config$estimators) || length(config$estimators) == 0) {
    config$estimators <- "REML"  # Safe fallback
  }
  
  # Basic pooled analysis with transport weights
  if (config$output$verbose) message("  Computing transport-weighted analysis...")
  
  tryCatch({
    fit_transport <- safe_robust_rma(
      yi = data$yi, 
      sei = data$se, 
      data = data,
      method = config$estimators[1],  # Use first estimator as default
      weights = if ("analysis_weights" %in% names(data)) data$analysis_weights else NULL,
      use_hksj = config$methods$hksj
    )
    
    results$transport <- fit_transport
    
  }, error = function(e) {
    warning("Transport-weighted analysis failed: ", e$message)
    
    # Try basic analysis without weights
    tryCatch({
      fit_basic <- safe_robust_rma(
        yi = data$yi, 
        sei = data$se, 
        data = data,
        method = "REML",  # Most reliable method
        use_hksj = FALSE  # Disable HKSJ if causing issues
      )
      results$transport <- fit_basic
      
    }, error = function(e2) {
      warning("Basic analysis also failed: ", e2$message)
      results$transport <- NULL
    })
  })
  
  return(results)
}
