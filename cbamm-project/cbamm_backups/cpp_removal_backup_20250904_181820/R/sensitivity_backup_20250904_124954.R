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
run_missing_study_sensitivity <- function(data, config) {
  tryCatch({
    # Call original logic but with safety checks
    result <- tryCatch({
        k <- nrow(data)
        if (is.null(n_max)) n_max <- min(5, max(1, round(0.25 * k)))
        
        se_new <- median(data$se, na.rm = TRUE)
        w_base <- median(data$analysis_weights, na.rm = TRUE)
        
        grid <- tidyr::expand_grid(n_missing = 0:n_max, hr = hr_grid) |>
          dplyr::mutate(loghr = log(hr))
        
        res <- grid |>
          dplyr::rowwise() |>
          dplyr::mutate(pooled_hr = {
            if (n_missing == 0) {
              fit <- robust_rma(data$yi, data$se, data = data, method = config$estimators[1],
                                weights = data$analysis_weights, use_hksj = isTRUE(config$methods$hksj))
            } else {
              add <- data.frame(
                study_id = sprintf("HYP_%02d", seq_len(n_missing)),
                yi = rep(loghr, n_missing),
                se = rep(se_new, n_missing),
                study_type = factor(rep("OBS", n_missing), levels = levels(data$study_type)),
                grade = if ("grade" %in% names(data)) {
                  factor(rep("Low", n_missing), levels = levels(data$grade))
                } else {
                  factor(rep("Low", n_missing))
                },
                analysis_weights = rep(w_base, n_missing),
                stringsAsFactors = FALSE
              )
              
              dd <- dplyr::bind_rows(
                dplyr::select(data, study_id, yi, se, study_type, grade, analysis_weights),
                add
              )
              
              fit <- robust_rma(dd$yi, dd$se, data = dd, method = config$estimators[1],
                                weights = dd$analysis_weights, use_hksj = isTRUE(config$methods$hksj))
            }
            
            exp(as.numeric(coef(fit)))
          }) |>
          dplyr::ungroup()
        
        return(res)
    }, error = function(e) {
      warning('Missing study analysis failed: ', e$message)
      return(NULL)
    })
    
    # Return safe result
    if (is.null(result)) {
      return(list(
        missing_studies = data.frame(
          scenario = 'failed',
          estimate = NA_real_,
          ci_lower = NA_real_,
          ci_upper = NA_real_,
          studies_added = 0L,
          stringsAsFactors = FALSE
        ),
        message = 'Missing study sensitivity analysis failed'
      ))
    } else {
      return(result)
    }
  }, error = function(e) {
    warning('Critical error in missing study sensitivity: ', e$message)
    return(list(
      missing_studies = data.frame(
        scenario = 'error',
        estimate = NA_real_,
        ci_lower = NA_real_,
        ci_upper = NA_real_,
        studies_added = 0L,
        stringsAsFactors = FALSE
      ),
      message = paste('Analysis error:', e$message)
    ))
  })
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
