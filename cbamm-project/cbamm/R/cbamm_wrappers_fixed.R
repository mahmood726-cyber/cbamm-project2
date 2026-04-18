.normalize_cbamm_public_data <- function(data,
                                         effect_col = "yi",
                                         se_col = "se",
                                         study_col = "study_id") {
  if (!is.data.frame(data)) {
    stop("Data must be a data.frame")
  }

  normalized <- data

  if (!identical(effect_col, "yi") &&
      effect_col %in% names(normalized) &&
      !"yi" %in% names(normalized)) {
    normalized$yi <- normalized[[effect_col]]
  }

  if (!identical(se_col, "se") &&
      se_col %in% names(normalized) &&
      !"se" %in% names(normalized)) {
    normalized$se <- normalized[[se_col]]
  }

  if (!identical(study_col, "study_id") &&
      study_col %in% names(normalized) &&
      !"study_id" %in% names(normalized)) {
    normalized$study_id <- normalized[[study_col]]
  }

  if (!inherits(normalized, "cbamm_data")) {
    normalized <- standardize_cbamm_data(normalized)
  }

  normalized
}

.normalize_cbamm_mode <- function(mode) {
  mode <- tolower(mode[1])

  if (mode %in% c("quick", "memory_efficient")) {
    return("fast")
  }

  if (mode == "standard") {
    return("balanced")
  }

  if (mode == "full") {
    return("comprehensive")
  }

  match.arg(mode, c("fast", "balanced", "comprehensive"))
}

#' CBAMM Fast Analysis with Auto-Standardization
#' 
#' Fast meta-analysis with automatic data standardization
#' 
#' @param data Data frame containing effect sizes and standard errors
#' @param method Estimation method (default "REML")
#' @param ... Additional arguments
#' @return cbamm object with meta-analysis results
#' @export
#' 
#' @examples
#' # Works with any column naming convention
#' data <- data.frame(TE = rnorm(10), seTE = runif(10, 0.1, 0.3))
#' result <- cbamm_fast(data)
#' summary(result)
#' 
cbamm_fast <- function(data,
                       config = cbamm_config(),
                       target_population = NULL,
                       method = NULL,
                       ...) {
  normalized_data <- .normalize_cbamm_public_data(data)

  if (!is.null(method)) {
    config$estimators[1] <- method
  }

  result <- cbamm(
    data = normalized_data,
    config = config,
    target_population = target_population,
    performance_mode = "fast"
  )

  class(result) <- unique(c("cbamm_fast", class(result)))
  result
}

.cbamm_fast_public <- cbamm_fast

#' Internal Fast Implementation
#' @noRd
.cbamm_fast_internal <- function(data, method = "REML", ...) {
  n <- nrow(data)
  
  if (n < 2) {
    warning("Less than 2 studies - heterogeneity cannot be estimated")
  }
  
  yi <- data$yi
  vi <- if ("vi" %in% names(data)) data$vi else data$se^2
  wi <- 1/vi
  
  # Weighted mean effect
  theta <- sum(wi * yi) / sum(wi)
  
  # Q statistic
  Q <- sum(wi * (yi - theta)^2)
  df <- n - 1
  
  # Between-study variance (tau-squared)
  if (method == "DL") {
    # DerSimonian-Laird
    tau2 <- max(0, (Q - df) / (sum(wi) - sum(wi^2)/sum(wi)))
  } else {
    # Simplified REML
    tau2 <- max(0, (Q - df) / (sum(wi) - sum(wi^2)/sum(wi)))
    
    # Iterative refinement for REML
    for (i in 1:10) {
      wi_star <- 1/(vi + tau2)
      tau2_new <- max(0, sum(wi_star^2 * ((yi - sum(wi_star * yi)/sum(wi_star))^2 - vi)) / sum(wi_star^2))
      if (abs(tau2_new - tau2) < 0.0001) break
      tau2 <- tau2_new
    }
  }
  
  # Final estimates with tau2
  wi_star <- 1/(vi + tau2)
  theta_final <- sum(wi_star * yi) / sum(wi_star)
  se_theta <- sqrt(1/sum(wi_star))
  
  # I-squared
  I2 <- if (Q > df) 100 * (Q - df) / Q else 0
  
  # Create result object
  result <- list(
    b = theta_final,
    beta = theta_final,  # Alternative name
    se = se_theta,
    tau2 = tau2,
    tau = sqrt(tau2),
    k = n,
    data = data,
    method = method,
    QE = Q,
    QEp = if (df > 0) pchisq(Q, df, lower.tail = FALSE) else NA,
    I2 = I2,
    H2 = Q / df,
    ci.lb = theta_final - 1.96 * se_theta,
    ci.ub = theta_final + 1.96 * se_theta,
    pi.lb = theta_final - 1.96 * sqrt(se_theta^2 + tau2),
    pi.ub = theta_final + 1.96 * sqrt(se_theta^2 + tau2),
    weights = wi_star / sum(wi_star) * 100
  )
  
  class(result) <- c("cbamm_fast", "cbamm")
  return(result)
}

#' CBAMM Optimized Analysis
#' 
#' @param data Data frame with meta-analysis data
#' @param mode Analysis mode: "fast", "balanced", or "comprehensive"
#' @param ... Additional arguments
#' @export
cbamm_optimized <- function(data,
                            mode = "fast",
                            config = cbamm_config(),
                            target_population = NULL,
                            effect_col = "yi",
                            se_col = "se",
                            study_col = "study_id",
                            ...) {
  performance_mode <- .normalize_cbamm_mode(mode)
  normalized_data <- .normalize_cbamm_public_data(
    data = data,
    effect_col = effect_col,
    se_col = se_col,
    study_col = study_col
  )

  result <- cbamm(
    data = normalized_data,
    config = config,
    target_population = target_population,
    performance_mode = performance_mode
  )

  result$diagnostics$requested_mode <- mode
  class(result) <- unique(c("cbamm_optimized", class(result)))
  result
}

#' Balanced Implementation
#' @noRd
.cbamm_balanced_internal <- function(data, ...) {
  # Start with fast analysis
  result <- .cbamm_fast_internal(data, ...)
  
  # Add studentized residuals
  yi <- data$yi
  vi <- if ("vi" %in% names(data)) data$vi else data$se^2
  resid <- (yi - result$b) / sqrt(vi + result$tau2)
  result$residuals <- resid
  
  # Add Cook's distances
  n <- nrow(data)
  cooks <- numeric(n)
  for (i in 1:n) {
    temp <- .cbamm_fast_internal(data[-i, ], ...)
    cooks[i] <- (result$b - temp$b)^2 / result$se^2
  }
  result$cooks.distance <- cooks
  
  class(result) <- c("cbamm_balanced", "cbamm")
  return(result)
}

#' Comprehensive Implementation
#' @noRd
.cbamm_comprehensive_internal <- function(data, ...) {
  # Start with balanced analysis
  result <- .cbamm_balanced_internal(data, ...)
  
  # Leave-one-out analysis
  n <- nrow(data)
  loo_results <- data.frame(
    study = data$study_id,
    estimate = numeric(n),
    se = numeric(n),
    tau2 = numeric(n)
  )
  
  for (i in 1:n) {
    loo <- .cbamm_fast_internal(data[-i, ], ...)
    loo_results$estimate[i] <- loo$b
    loo_results$se[i] <- loo$se
    loo_results$tau2[i] <- loo$tau2
  }
  
  result$leave1out <- loo_results
  
  # Baujat plot data
  yi <- data$yi
  vi <- if ("vi" %in% names(data)) data$vi else data$se^2
  wi <- 1/(vi + result$tau2)
  contrib_Q <- wi * (yi - result$b)^2
  contrib_overall <- wi * abs(yi - result$b)
  
  result$baujat <- data.frame(
    x = contrib_overall,
    y = contrib_Q,
    study = data$study_id
  )
  
  class(result) <- c("cbamm_comprehensive", "cbamm")
  return(result)
}

#' Print method for cbamm objects
#' @export
print.cbamm <- function(x, digits = 4, ...) {
  cat("
Meta-Analysis Results (CBAMM)
")
  cat("==============================
")
  cat("Method:", x$method, "
")
  cat("Studies:", x$k, "

")
  
  cat("Summary Effect:
")
  cat("  Estimate:", round(x$b, digits), "
")
  cat("  SE:", round(x$se, digits), "
")
  cat("  95% CI: [", round(x$ci.lb, digits), ", ", round(x$ci.ub, digits), "]
")
  cat("  z =", round(x$b/x$se, 2), ", p =", format.pval(2*pnorm(-abs(x$b/x$se)), digits = 3), "

")
  
  cat("Heterogeneity:
")
  cat("  tau²:", round(x$tau2, digits), "
")
  cat("  tau:", round(x$tau, digits), "
")
  cat("  I²:", round(x$I2, 1), "%
")
  cat("  H²:", round(x$H2, 2), "
")
  cat("  Q(df=", x$k-1, ") = ", round(x$QE, 2), ", p = ", format.pval(x$QEp, digits = 3), "

", sep="")
  
  cat("Prediction Interval:
")
  cat("  95% PI: [", round(x$pi.lb, digits), ", ", round(x$pi.ub, digits), "]
")
  
  invisible(x)
}

#' Summary method for cbamm objects
#' @export
summary.cbamm <- function(object, ...) {
  print(object, ...)
  invisible(object)
}

