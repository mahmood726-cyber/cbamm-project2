
#' Enhanced Multivariate Meta-Analysis
#' 
#' @param data Data frame with multiple outcomes
#' @param outcomes Character vector of outcome column names
#' @param se_columns Character vector of SE column names
#' @param correlation Correlation matrix between outcomes
#' @param method Estimation method
#' @export
multivariate_meta_enhanced <- function(data, 
                                       outcomes,
                                       se_columns,
                                       correlation = NULL,
                                       method = "REML") {
  
  # Input validation
  if (length(outcomes) != length(se_columns)) {
    stop("Number of outcomes must match number of SE columns", call. = FALSE)
  }
  
  if (!all(outcomes %in% names(data))) {
    stop("Not all outcome columns found in data", call. = FALSE)
  }
  
  if (!all(se_columns %in% names(data))) {
    stop("Not all SE columns found in data", call. = FALSE)
  }
  
  k <- nrow(data)
  p <- length(outcomes)
  
  if (k < 3) {
    stop("At least 3 studies required for multivariate meta-analysis", call. = FALSE)
  }
  
  # Extract outcome data
  Y <- as.matrix(data[, outcomes])
  SE <- as.matrix(data[, se_columns])
  
  # Initialize between-study covariance
  Tau <- diag(rep(0.05, p))
  
  # Simple estimation
  wi <- 1/rowMeans(SE^2)
  
  # Pooled estimates per outcome
  theta <- numeric(p)
  se_theta <- numeric(p)
  
  for (j in 1:p) {
    wj <- 1/SE[,j]^2
    theta[j] <- sum(wj * Y[,j]) / sum(wj)
    
    # Heterogeneity
    Qj <- sum(wj * (Y[,j] - theta[j])^2)
    tau2j <- max(0, (Qj - (k-1)) / (sum(wj) - sum(wj^2)/sum(wj)))
    Tau[j,j] <- tau2j
    
    # SE of pooled estimate
    wj_star <- 1/(SE[,j]^2 + tau2j)
    se_theta[j] <- sqrt(1/sum(wj_star))
  }
  
  # Add correlation structure if provided
  if (!is.null(correlation) && p > 1) {
    for (i in 1:(p-1)) {
      for (j in (i+1):p) {
        Tau[i,j] <- Tau[j,i] <- correlation[i,j] * sqrt(Tau[i,i] * Tau[j,j])
      }
    }
  }
  
  # Calculate I-squared
  avg_vi <- mean(SE^2)
  I2 <- 100 * mean(diag(Tau)) / (mean(diag(Tau)) + avg_vi)
  
  # Results
  result <- list(
    coefficients = theta,
    se = se_theta,
    Tau = Tau,
    outcomes = outcomes,
    k = k,
    p = p,
    method = method,
    I2 = I2,
    ci.lb = theta - 1.96 * se_theta,
    ci.ub = theta + 1.96 * se_theta,
    data = data
  )
  
  class(result) <- c("cbamm_multivariate", "cbamm")
  return(result)
}

#' Print Multivariate Results
#' @export
print.cbamm_multivariate <- function(x, digits = 4, ...) {
  cat("
Multivariate Meta-Analysis Results
")
  cat("===================================
")
  cat("Studies:", x$k, "
")
  cat("Outcomes:", x$p, "
")
  cat("Method:", x$method, "

")
  
  cat("Pooled Estimates:
")
  for (i in 1:x$p) {
    cat(sprintf("  %s: %.3f (SE=%.3f, 95%% CI: [%.3f, %.3f])
",
                x$outcomes[i], x$coefficients[i], x$se[i], 
                x$ci.lb[i], x$ci.ub[i]))
  }
  
  cat("
Heterogeneity:
")
  cat(sprintf("  I² = %.1f%%
", x$I2))
  
  cat("
Between-study covariance:
")
  print(round(x$Tau, digits))
  
  invisible(x)
}

