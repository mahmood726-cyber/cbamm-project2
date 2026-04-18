#' Run Sensitivity Analyses
#'
#' @param data Study data
#' @param config CBAMM configuration
#'
#' @return List of sensitivity analysis results
#' @keywords internal
run_sensitivity_analyses <- function(data, config) {
  results <- list()
  
  # PET-PEESE
  if (isTRUE(config$methods$pet_peese)) {
    if (config$output$verbose) message("  Running PET-PEESE...")
    results$pet_peese <- pet_peese(data$yi, data$se)
  }
  
  # Publication bias tests
  if (any(config$sensitivity$publication_bias_methods %in% c("egger", "begg", "trim_fill"))) {
    if (config$output$verbose) message("  Running publication bias tests...")
    results$pub_bias <- run_publication_bias_sensitivity(data)
  }
  
  # Missing studies sensitivity
  if (isTRUE(config$methods$missing_studies)) {
    if (config$output$verbose) message("  Running missing studies sensitivity...")
    results$missing_studies <- run_missing_study_sensitivity(data, config)
  }
  
  return(results)
}

#' PET-PEESE Analysis
#'
#' @param yi Effect sizes
#' @param sei Standard errors
#'
#' @return Vector with PET and PEESE intercepts
#' @export
pet_peese <- function(yi, sei) {
  w <- 1 / (sei^2)
  df <- data.frame(yi = yi, sei = sei, w = w)
  
  fit_pet   <- lm(yi ~ sei, data = df, weights = w)
  fit_peese <- lm(yi ~ I(sei^2), data = df, weights = w)
  
  c(PET = unname(coef(fit_pet)[1]), PEESE = unname(coef(fit_peese)[1]))
}

#' Run Publication Bias Sensitivity Analysis
#'
#' @param data Study data
#'
#' @return List of publication bias test results
#' @keywords internal
run_publication_bias_sensitivity <- function(data) {
  if (nrow(data) < 5) {
    message("Too few studies for publication bias tests")
    return(NULL)
  }
  
  results <- list()
  
  # Trim-and-fill
  tf <- try(metafor::trimfill(metafor::rma(yi = data$yi, sei = data$se, method = "REML")), silent = TRUE)
  if (!inherits(tf, "try-error")) {
    results$trimfill <- list(
      k0 = tf$k0,
      estimate = as.numeric(coef(tf)),
      se = as.numeric(tf$se),
      ci_lb = tf$ci.lb,
      ci_ub = tf$ci.ub
    )
  }
  
  # Selection model (weightr)
  if (check_package("weightr", quietly = TRUE)) {
    pvals <- 2 * pnorm(-abs(data$yi / data$se))
    br <- .cbamm_build_weightr_breaks(pvals)
    
    sm <- try(weightr::weightfunct(effect = data$yi, v = data$se^2, steps = br, table = FALSE), silent = TRUE)
    if (!inherits(sm, "try-error")) {
      results$selection_model <- list(
        breaks = br,
        estimate = sm$adj_est,
        se = sm$adj_se,
        ci_lb = sm$adj_ci.lb,
        ci_ub = sm$adj_ci.ub
      )
    }
  }
  
  return(results)
}

#' Build Weightr Breaks
#' @keywords internal
.cbamm_build_weightr_breaks <- function(pvals) {
  cand <- list(
    c(0, .025, .05, 1),
    c(0, quantile(pvals, .025, na.rm = TRUE), .05, 1),
    c(0, .01, .05, 1),
    c(0, quantile(pvals, .10, na.rm = TRUE), .25, 1)
  )
  
  for (br in cand) {
    idx <- cut(pvals, breaks = br, include.lowest = TRUE, right = TRUE)
    if (all(table(idx) > 0)) return(br)
  }
  
  br <- unique(c(0, quantile(pvals, probs = c(.1,.2,.3,.4,.5,.6,.7,.8,.9), na.rm = TRUE), 1))
  if (length(br) < 3) br <- c(0, .05, 1)
  br
}

#' Run Missing Study Sensitivity Analysis
#'
#' @param data Study data
#' @param config CBAMM configuration
#' @param n_max Maximum number of missing studies to test
#' @param hr_grid Grid of hazard ratios to test for missing studies
#'
#' @return Data frame with sensitivity results
#' @keywords internal
# CBAMM AUTO-FIX: Safe wrapper for missing study sensitivity
#' Run Missing Study Sensitivity Analysis (FIXED VERSION)
#'
#' @param data Study data
#' @param config CBAMM configuration
#' @param n_max Maximum number of missing studies to test
#' @param hr_grid Grid of hazard ratios to test for missing studies
#'
#' @return Data frame with sensitivity results
#' @keywords internal
run_missing_study_sensitivity <- function(data, config, n_max = NULL, hr_grid = c(0.8, 1.0, 1.2)) {
  k <- nrow(data)
  if (is.null(n_max)) n_max <- min(3, max(1, round(0.15 * k)))  # Reduced for speed
  
  se_new <- median(data$se, na.rm = TRUE)
  w_base <- if ('analysis_weights' %in% names(data)) {
    median(data$analysis_weights, na.rm = TRUE)
  } else {
    1
  }
  
  # Simple loop instead of tidyr to avoid data.frame issues
  results <- list()
  idx <- 1
  
  for (n_missing in 0:n_max) {
    for (hr in hr_grid) {
      tryCatch({
        if (n_missing == 0) {
          # Use original data
          fit <- robust_rma(data$yi, data$se, data = data,
                            method = config$estimators[1],
                            weights = if ('analysis_weights' %in% names(data)) data$analysis_weights else NULL,
                            use_hksj = isTRUE(config$methods$hksj))
          pooled_hr <- exp(as.numeric(coef(fit)))
        } else {
          # Add hypothetical studies
          add_studies <- data.frame(
            study_id = sprintf('HYP_%02d', seq_len(n_missing)),
            yi = rep(log(hr), n_missing),
            se = rep(se_new, n_missing),
            study_type = factor(rep('OBS', n_missing), levels = levels(data$study_type)),
            analysis_weights = rep(w_base, n_missing)
          )
          
          # Only use common columns
          common_cols <- intersect(names(data), names(add_studies))
          combined_data <- rbind(data[, common_cols], add_studies[, common_cols])
          
          fit <- robust_rma(combined_data$yi, combined_data$se,
                            method = config$estimators[1],
                            weights = if ('analysis_weights' %in% common_cols) combined_data$analysis_weights else NULL,
                            use_hksj = isTRUE(config$methods$hksj))
          pooled_hr <- exp(as.numeric(coef(fit)))
        }
        
        results[[idx]] <- data.frame(n_missing = n_missing, hr = hr, pooled_hr = pooled_hr)
        idx <- idx + 1
        
      }, error = function(e) {
        results[[idx]] <- data.frame(n_missing = n_missing, hr = hr, pooled_hr = NA)
        idx <<- idx + 1
      })
    }
  }
  
  # Combine results
  do.call(rbind, results)
}
#' Run RVE Analysis and Print Results
#'
#' @param res Meta-analysis result
#' @param cluster_vec Cluster vector for RVE
#' @param label Label for output
#'
#' @return RVE results (invisibly)
#' @keywords internal
rve_print <- function(res, cluster_vec, label = "[RVE-CR2]") {
  if (!check_package("clubSandwich", quietly = TRUE)) return(invisible(NULL))
  
  vc <- try(clubSandwich::vcovCR(res, cluster = cluster_vec, type = "CR2"), silent = TRUE)
  if (inherits(vc, "try-error")) return(invisible(NULL))
  
  rob <- try(clubSandwich::coef_test(res, vcov = vc), silent = TRUE)
  if (!inherits(rob, "try-error")) {
    message(label, " robust test for ", deparse(substitute(res)), ":")
    print(rob)
  }
  
  invisible(rob)
}
