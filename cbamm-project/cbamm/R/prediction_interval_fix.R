
#' Calculate Prediction Interval Wrapper
#' @export
calculate_prediction_interval_fixed <- function(result, level = 0.95) {
  # Extract what we need from result
  if (!is.null(result$core_results)) {
    est <- result$core_results$estimate
    se <- result$core_results$se
    tau <- result$core_results$tau
  } else if (!is.null(result$summary)) {
    est <- result$summary$estimate
    se <- result$summary$se
    tau <- result$summary$tau
  } else {
    stop("Cannot extract estimates from result")
  }
  
  # If original function expects different args
  if (length(formals(calculate_prediction_interval)) == 3) {
    return(calculate_prediction_interval(est, se, tau))
  } else {
    # Calculate manually
    z <- qnorm(1 - (1 - level) / 2)
    pred_se <- sqrt(se^2 + tau^2)
    return(c(est - z * pred_se, est + z * pred_se))
  }
}

