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
  
  # Multiverse analysis
  if (isTRUE(config$methods$multiverse)) {
    if (config$output$verbose) message("  Running multiverse analysis...")
    results$multiverse <- run_multiverse_analysis(data, config)
  }
  
  return(results)
}

#' PET-PEESE Analysis
#'
#' @param yi Effect sizes
#' @param sei Standard errors
#'
#' @return Vector with PET and PEESE intercepts
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
run_missing_study_sensitivity <- function(data, config, n_max = NULL, hr_grid = c(0.8, 1.0, 1.2)) {
  k <- nrow(data)
  if (is.null(n_max)) n_max <- min(3, max(1, round(0.15 * k)))
  
  se_new <- median(data$se, na.rm = TRUE)
  w_base <- if ('analysis_weights' %in% names(data)) {
    median(data$analysis_weights, na.rm = TRUE)
  } else {
    1
  }
  
  results <- list()
  idx <- 1
  
  for (n_missing in 0:n_max) {
    for (hr in hr_grid) {
      tryCatch({
        if (n_missing == 0) {
          fit <- robust_rma(data$yi, data$se, data = data,
                            method = config$estimators[1],
                            weights = if ('analysis_weights' %in% names(data)) data$analysis_weights else NULL,
                            use_hksj = isTRUE(config$methods$hksj))
          pooled_hr <- exp(as.numeric(coef(fit)))
        } else {
          add_studies <- data.frame(
            study_id = sprintf('HYP_%02d', seq_len(n_missing)),
            yi = rep(log(hr), n_missing),
            se = rep(se_new, n_missing),
            study_type = factor(rep('OBS', n_missing), levels = levels(data$study_type)),
            analysis_weights = rep(w_base, n_missing)
          )
          
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
  
  do.call(rbind, results)
}

#' Run Multiverse Analysis
#'
#' @param data Study data
#' @param config CBAMM configuration
#'
#' @return Data frame with multiverse results
#' @export
run_multiverse_analysis <- function(data, config) {
  
  # Multiverse factors
  estimators <- c("REML", "DL", "PM", "EB")
  subsets <- list(
    all = 1:nrow(data),
    rct_only = if("study_type" %in% names(data)) which(data$study_type == "RCT") else NULL
  )
  subsets <- subsets[!sapply(subsets, is.null)]
  
  hksj_options <- c(TRUE, FALSE)
  
  # Grid of configurations
  grid <- expand.grid(
    estimator = estimators,
    subset_name = names(subsets),
    hksj = hksj_options,
    stringsAsFactors = FALSE
  )
  
  # Run meta-analysis for each grid point
  # Parallelization point
  if (isTRUE(config$methods$parallel) && check_package("future.apply", quietly = TRUE)) {
    results_list <- future.apply::future_lapply(1:nrow(grid), function(i) {
      .run_multiverse_step(grid[i, ], subsets, data, config)
    }, future.seed = TRUE)
  } else {
    results_list <- lapply(1:nrow(grid), function(i) {
      .run_multiverse_step(grid[i, ], subsets, data, config)
    })
  }
  
  results <- do.call(rbind, results_list)
  class(results) <- c("cbamm_multiverse", "data.frame")
  return(results)
}

#' Run Single Multiverse Step
#' @keywords internal
.run_multiverse_step <- function(row, subsets, data, config) {
  idx <- subsets[[row$subset_name]]
  
  if (length(idx) < 3) return(NULL)
  
  d_sub <- data[idx, ]
  
  fit <- try(robust_rma(
    yi = d_sub$yi, 
    sei = d_sub$se, 
    method = row$estimator,
    weights = if("analysis_weights" %in% names(d_sub)) d_sub$analysis_weights else NULL,
    use_hksj = row$hksj
  ), silent = TRUE)
  
  if (inherits(fit, "try-error")) return(NULL)
  
  data.frame(
    estimator = row$estimator,
    subset = row$subset_name,
    hksj = row$hksj,
    estimate = as.numeric(coef(fit)),
    se = fit$se,
    ci_lb = fit$ci.lb,
    ci_ub = fit$ci.ub,
    tau2 = fit$tau2,
    I2 = fit$I2,
    k = fit$k
  )
}

#' Run RVE Analysis and Print Results
#'
#' @param res Meta-analysis result
#' @param cluster_vec Cluster vector for RVE
#' @param label Label for output
#'
#' @return RVE results (invisibly)
#' @export
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
