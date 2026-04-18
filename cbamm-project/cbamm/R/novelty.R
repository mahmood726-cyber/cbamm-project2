#' Run Adaptive Advisor
#'
#' @param data Study data
#' @param pooled_results Pooled meta-analysis results
#' @param config CBAMM configuration
#'
#' @return List with advisor recommendations
#' @keywords internal
run_adaptive_advisor <- function(data, pooled_results, config) {
  k <- nrow(data)
  dup <- any(duplicated(data$study_id))
  
  fit <- pooled_results$transport
  I2 <- if (!is.null(fit)) fit$I2 else NA_real_
  
  eg <- try(metafor::regtest(fit, model = "lm"), silent = TRUE)
  bias_flag <- (!inherits(eg, "try-error")) && is.finite(eg$pval) && eg$pval < 0.10
  
  rec <- c()
  if (k < 5) {
    rec <- c(rec, "Fixed-effects is defensible (very small k); otherwise keep REML+HKSJ.")
  } else {
    rec <- c(rec, "Random-effects (REML+HKSJ) appropriate.")
  }
  
  if (is.finite(I2) && I2 > 75) {
    rec <- c(rec, "High heterogeneity: inspect moderators/time trends.")
  }
  
  if (dup) {
    rec <- c(rec, "Potential dependent effects: report RVE (CR2) robust tests.")
  }
  
  if (bias_flag) {
    rec <- c(rec, "Funnel asymmetry: report PET–PEESE & selection-model sensitivities.")
  }
  
  if (all(c("age_mean","female_pct","bmi_mean","charlson") %in% names(data))) {
    rec <- c(rec, "Covariates present: keep transport weighting (WeightIt/fallback).")
  }
  
  # Check for HKSJ fallback in the results
  hksj_active <- config$methods$hksj
  hksj_failed <- hksj_active && !is.null(fit) && (is.null(fit$test) || fit$test != "knha")
  
  if (hksj_failed) {
    rec <- c(rec, "WARNING: HKSJ adjustment failed and reverted to Wald z. Confidence intervals may be anti-conservative.")
  }
  
  if (config$output$verbose) {
    message(sprintf("k = %d, I² ≈ %s%%", k, ifelse(is.finite(I2), sprintf("%.1f", I2), "NA")))
    if (!inherits(eg, "try-error")) message(sprintf("Egger test (lm) p = %.3f", eg$pval))
    
    inf_method <- if (hksj_failed) "Wald z (HKSJ Fallback)" else if (hksj_active) "HKSJ (Knapp–Hartung)" else "Wald z"
    message("Inference: ", inf_method)
    
    message("Recommendations:")
    for (r in rec) message(" - ", r)
  }
  
  list(
    k = k, 
    I2 = I2, 
    egger_p = ifelse(inherits(eg,"try-error"), NA, eg$pval), 
    recommendations = rec,
    has_duplicates = dup,
    bias_detected = bias_flag
  )
}

#' Generate Advisor Recommendations
#'
#' @param results Full CBAMM results object
#' @param config CBAMM configuration
#'
#' @return List of recommendations
#' @keywords internal
generate_advisor_recommendations <- function(results, config) {
  
  recommendations <- list(
    methodological = c(),
    interpretation = c(),
    reporting = c()
  )
  
  # Extract key information
  k <- results$data_summary$n_studies
  has_heterogeneity <- FALSE
  has_bias <- FALSE
  
  if (!is.null(results$meta_results$transport)) {
    I2 <- results$meta_results$transport$I2
    has_heterogeneity <- !is.na(I2) && I2 > 50
  }
  
  # Methodological recommendations
  if (k < 5) {
    recommendations$methodological <- c(
      recommendations$methodological,
      "Consider fixed-effects model given small number of studies"
    )
  }
  
  if (has_heterogeneity) {
    recommendations$methodological <- c(
      recommendations$methodological,
      "Investigate sources of heterogeneity using subgroup analysis or meta-regression"
    )
  }
  
  if (!is.null(results$conflict_detection)) {
    recommendations$methodological <- c(
      recommendations$methodological,
      "Study conflicts detected - consider sensitivity analyses"
    )
  }
  
  # Interpretation recommendations
  if (!is.null(results$transport_results)) {
    recommendations$interpretation <- c(
      recommendations$interpretation,
      "Results adjusted for transportability to target population"
    )
  }
  
  # Reporting recommendations
  recommendations$reporting <- c(
    "Report both fixed and random effects results",
    "Include prediction intervals for clinical interpretation",
    "Provide forest plot and funnel plot"
  )
  
  if (!is.null(results$sensitivity_results)) {
    recommendations$reporting <- c(
      recommendations$reporting,
      "Report sensitivity analyses for publication bias"
    )
  }
  
  return(recommendations)
}

#' Detect Study Conflicts using K-means Clustering
#'
#' @param data Study data
#' @param threshold Minimum difference threshold for conflict detection
#' @param k_candidates Vector of K values to try
#'
#' @return Conflict detection results
#' @keywords internal

detect_study_conflicts <- function(data, threshold = 0.15, k_candidates = 2:4) {
  if (nrow(data) < 6) {
    if (length(intersect(names(data), c("analysis_weights", "analysis_weights_grade"))) == 0) {
      data$analysis_weights <- rep(1, nrow(data))
    }
    return(NULL)
  }
  
  set.seed(123)
  yi <- data$yi
  valid_k <- intersect(as.integer(k_candidates), 2:min(10L, floor(nrow(data) / 2)))
  if (length(valid_k) == 0) {
    return(NULL)
  }
  
  # Dynamic K selection using silhouette analysis
  pick_k <- valid_k[1]
  # Optional silhouette scoring is only helpful when there is enough signal to
  # discriminate between candidate solutions. On tiny datasets it adds overhead
  # and has been brittle in headless Windows Rscript runs.
  if (nrow(data) >= 20L && length(valid_k) > 1L && check_package("cluster", quietly = TRUE)) {
    best_k <- pick_k
    best_s <- -Inf
    
    for (K in valid_k) {
      km_try <- try(kmeans(yi, centers = K, nstart = 20), silent = TRUE)
      if (inherits(km_try, "try-error")) next
      
      sil <- try(cluster::silhouette(km_try$cluster, dist(yi)), silent = TRUE)
      if (!inherits(sil, "try-error")) {
        s <- mean(sil[, 3], na.rm = TRUE)
        if (is.finite(s) && s > best_s) { 
          best_s <- s
          best_k <- K 
        }
      }
    }
    pick_k <- best_k
  }
  
  km <- try(kmeans(yi, centers = pick_k, nstart = 20), silent = TRUE)
  if (inherits(km, "try-error")) return(NULL)
  
  # Add cluster information to data
  df <- data
  df$cluster <- factor(km$cluster)
  
  # Ensure weights are available
  w_all <- if ("analysis_weights_grade" %in% names(df)) {
    with(df, (1 / (se^2)) * ifelse(is.na(analysis_weights_grade), 
                                   ifelse(is.na(analysis_weights), 1, analysis_weights),
                                   analysis_weights_grade))
  } else if ("analysis_weights" %in% names(df)) {
    with(df, (1 / (se^2)) * analysis_weights)
  } else {
    1 / (df$se^2)
  }
  
  df$w <- w_all
  
  # Cluster summary using base R
  clusters <- split(df, df$cluster)
  comp_list <- lapply(clusters, function(cluster_data) {
    data.frame(
      cluster = cluster_data$cluster[1],
      n = nrow(cluster_data),
      mean_logHR = weighted.mean(cluster_data$yi, cluster_data$w, na.rm = TRUE),
      HR = exp(weighted.mean(cluster_data$yi, cluster_data$w, na.rm = TRUE)),
      mean_year = if ("year" %in% names(cluster_data)) mean(cluster_data$year, na.rm = TRUE) else NA_real_,
      prop_RCT = if ("study_type" %in% names(cluster_data)) mean(cluster_data$study_type == "RCT", na.rm = TRUE) else NA_real_,
      prop_OBS = if ("study_type" %in% names(cluster_data)) mean(cluster_data$study_type == "OBS", na.rm = TRUE) else NA_real_,
      prop_MR = if ("study_type" %in% names(cluster_data)) mean(cluster_data$study_type == "MR", na.rm = TRUE) else NA_real_,
      mean_grade = if ("grade" %in% names(cluster_data)) mean(as.numeric(cluster_data$grade), na.rm = TRUE) else NA_real_
    )
  })
  
  comp <- do.call(rbind, comp_list)
  centers <- comp$mean_logHR
  delta <- if (length(centers) >= 2) diff(range(centers)) else 0
  
  if (length(centers) < 2 || !is.finite(delta) || delta < threshold) {
    return(NULL)
  }
  
  # Explanation model for binary clustering
  mod <- NULL
  importance <- NULL
  
  if (pick_k == 2L) {
    df$cl_bin <- as.integer(df$cluster == levels(df$cluster)[2])
    
    # Build formula based on available variables
    predictors <- c()
    has_variation <- function(x) length(unique(x[!is.na(x)])) > 1L
    if ("study_type" %in% names(df) && has_variation(df$study_type)) predictors <- c(predictors, "study_type")
    if ("year" %in% names(df) && has_variation(df$year)) predictors <- c(predictors, "year")
    if ("grade" %in% names(df) && has_variation(df$grade)) predictors <- c(predictors, "as.numeric(grade)")
    
    if (length(unique(df$cl_bin)) > 1L && length(predictors) > 0) {
      # Use scale(year) in GLM for stability
      formula_glm <- as.formula(paste("cl_bin ~", 
                                     paste(gsub("^year$", "scale(year)", predictors), 
                                           collapse = " + ")))
      
      mod <- try(glm(formula_glm, data = df, family = binomial()), silent = TRUE)
      if (inherits(mod, "try-error")) mod <- NULL
      
      # Random forest importance is optional; skip it for small datasets where
      # it adds little value and can slow headless verification.
      if (nrow(df) >= 20L && length(predictors) > 1L &&
          check_package("randomForest", quietly = TRUE)) {
        formula_rf <- as.formula(paste("factor(cl_bin) ~", paste(predictors, collapse = " + ")))
        rf_mod <- try(randomForest::randomForest(formula_rf, data = df, ntree = 500), silent = TRUE)
        if (!inherits(rf_mod, "try-error")) {
          importance <- randomForest::importance(rf_mod)
        }
      }
    }
  }
  
  list(
    kmeans = km,
    summary = comp,
    model = mod,
    importance = importance,
    K = pick_k,
    delta = delta,
    data = df,
    threshold_met = delta >= threshold
  )
}

#' Run Conflict Detection
#'
#' @param data Study data  
#' @param threshold Minimum difference threshold
#' @param k_candidates K values to try
#'
#' @return Conflict detection results
#' @keywords internal
run_conflict_detection <- function(data, threshold = 0.15, k_candidates = 2:4) {
  return(detect_study_conflicts(data, threshold, k_candidates))
}

#' Display CBAMM Novelty Notes
#'
#' @return Invisible NULL (prints notes)
cbamm_novelty_notes <- function() {
  message("
CBAMM v5.7 novelty (concise):

• Integration: Transport weighting + HKSJ + CR2 RVE + PET–PEESE + Bayesian stacking + multiverse in one pipeline.

• Advisor: Automated guidance from k, I², Egger's test, duplicates (RVE) — not common in meta tools.

• Conflict detection: Data-driven clustering with explanation model; dynamic K (2–4) via silhouette when available.

• Missing-study sensitivity: Transparent grid augmentation (+ optional trimfill/selection if packages present).

• Multiverse: Systematic τ² × subset × weighting sweep encourages analytic transparency.

• Transportability fallback: WeightIt entropy balancing → robust entropy optimizer.

• Bayesian stacking: LOO-based weights across brms models; JAGS fallback ensures portability.

v5.7 fixes:

• JAGS: use update(jm, ...) S3 (not rjags::update); use as.matrix(sm) or coda::as.matrix.mcmc.list(sm).

• Funnel + PET/PEESE/LOO/Cumulative/Conflict/Missing/Bayesian plots all hardened.

• Helper to FORCE-DISPLAY & SAVE every plot even in headless runs.
")
  
  invisible(NULL)
}
