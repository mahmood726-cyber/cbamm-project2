
# Cumulative Meta-Analysis Implementation for CBAMM Package
# This module provides comprehensive cumulative meta-analysis capabilities

#' Cumulative Meta-Analysis
#'
#' Performs cumulative meta-analysis to assess how evidence accumulates over time
#' as studies are added sequentially to the analysis.
#'
#' @param data Data frame containing study data with required columns
#' @param order_by Variable to order studies by. Options: "year", "sample_size", 
#'   "precision", "quality", "custom"
#' @param order_direction Direction of ordering: "ascending" or "descending"
#' @param custom_order Optional vector specifying custom study ordering
#' @param method Meta-analysis method: "random" (default), "fixed", or "bayesian"
#' @param confidence_level Confidence level for intervals (default: 0.95)
#' @param minimum_studies Minimum number of studies required for analysis (default: 2)
#' @param track_statistics Which statistics to track: "all", "basic", or custom vector
#' @return Object of class "cbamm_cumulative" containing cumulative results
#' @export
#' @examples
#' \dontrun{
#' # Basic cumulative analysis by publication year
#' result <- cumulative_meta_analysis(data, order_by = "year")
#' 
#' # Cumulative analysis by sample size (largest first)
#' result <- cumulative_meta_analysis(data, 
#'                                   order_by = "sample_size", 
#'                                   order_direction = "descending")
#' 
#' # Plot cumulative results
#' plot(result)
#' summary(result)
#' }
cumulative_meta_analysis <- function(data,
                                    order_by = "year",
                                    order_direction = "ascending",
                                    custom_order = NULL,
                                    method = "random",
                                    confidence_level = 0.95,
                                    minimum_studies = 2,
                                    track_statistics = "all") {
  
  # Input validation
  validate_cumulative_inputs(data, order_by, order_direction, method)
  
  # Prepare and order data
  ordered_data <- prepare_cumulative_data(data, order_by, order_direction, custom_order)
  
  # Perform cumulative analysis
  cumulative_results <- run_cumulative_analysis(
    ordered_data, method, confidence_level, minimum_studies, track_statistics
  )
  
  # Calculate stability metrics
  stability_metrics <- calculate_stability_metrics(cumulative_results)
  
  # Assess evidence sufficiency
  sufficiency_assessment <- assess_evidence_sufficiency(cumulative_results)
  
  # Create cumulative analysis object
  result <- list(
    data = ordered_data,
    results = cumulative_results,
    stability = stability_metrics,
    sufficiency = sufficiency_assessment,
    order_by = order_by,
    order_direction = order_direction,
    method = method,
    confidence_level = confidence_level,
    call = match.call()
  )
  
  class(result) <- "cbamm_cumulative"
  return(result)
}

#' Validate Cumulative Meta-Analysis Inputs
#'
#' Validates input parameters for cumulative meta-analysis
#'
#' @param data Input data frame
#' @param order_by Ordering variable
#' @param order_direction Ordering direction
#' @param method Analysis method
validate_cumulative_inputs <- function(data, order_by, order_direction, method) {
  
  # Check data frame
  if (!is.data.frame(data)) {
    stop("Data must be a data frame")
  }
  
  if (nrow(data) < 2) {
    stop("Cumulative meta-analysis requires at least 2 studies")
  }
  
  # Check required columns
  required_cols <- c("study_id", "yi", "se")
  missing_cols <- setdiff(required_cols, names(data))
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Check ordering variable
  valid_order_vars <- c("year", "sample_size", "precision", "quality", "custom")
  if (!order_by %in% valid_order_vars && !order_by %in% names(data)) {
    stop("order_by must be one of: ", paste(valid_order_vars, collapse = ", "), 
         " or a column name in data")
  }
  
  # Check ordering direction
  if (!order_direction %in% c("ascending", "descending")) {
    stop("order_direction must be 'ascending' or 'descending'")
  }
  
  # Check method
  if (!method %in% c("random", "fixed", "bayesian")) {
    stop("method must be 'random', 'fixed', or 'bayesian'")
  }
  
  # Check for missing values in effect sizes and standard errors
  if (any(is.na(data$yi))) {
    stop("Effect sizes (yi) cannot contain missing values")
  }
  
  if (any(is.na(data$se))) {
    stop("Standard errors (se) cannot contain missing values")
  }
  
  if (any(data$se <= 0)) {
    stop("Standard errors must be positive")
  }
}

#' Prepare Cumulative Data
#'
#' Orders data according to specified criteria for cumulative analysis
#'
#' @param data Input data frame
#' @param order_by Ordering variable
#' @param order_direction Ordering direction
#' @param custom_order Custom ordering vector
#' @return Ordered data frame
prepare_cumulative_data <- function(data, order_by, order_direction, custom_order) {
  
  # Create ordering variable if needed
  if (order_by == "precision") {
    data$precision <- 1 / data$se^2
    order_var <- "precision"
  } else if (order_by == "quality" && !"quality" %in% names(data)) {
    # Default quality assessment based on sample size and precision
    data$quality <- scale(log(data$sample_size)) + scale(1/data$se^2)
    order_var <- "quality"
  } else if (order_by == "custom") {
    if (is.null(custom_order)) {
      stop("custom_order must be provided when order_by = 'custom'")
    }
    if (length(custom_order) != nrow(data)) {
      stop("custom_order length must equal number of studies")
    }
    data$custom_order <- custom_order
    order_var <- "custom_order"
  } else {
    order_var <- order_by
  }
  
  # Check if ordering variable exists
  if (!order_var %in% names(data)) {
    stop("Ordering variable '", order_var, "' not found in data")
  }
  
  # Order data
  if (order_direction == "ascending") {
    ordered_data <- data[order(data[[order_var]]), ]
  } else {
    ordered_data <- data[order(data[[order_var]], decreasing = TRUE), ]
  }
  
  # Add cumulative study number
  ordered_data$cumulative_n <- 1:nrow(ordered_data)
  
  # Reset row names
  rownames(ordered_data) <- NULL
  
  return(ordered_data)
}

#' Run Cumulative Analysis
#'
#' Performs the sequential meta-analysis calculations
#'
#' @param ordered_data Ordered study data
#' @param method Analysis method
#' @param confidence_level Confidence level
#' @param minimum_studies Minimum studies required
#' @param track_statistics Statistics to track
#' @return Cumulative results data frame
run_cumulative_analysis <- function(ordered_data, method, confidence_level, 
                                   minimum_studies, track_statistics) {
  
  n_studies <- nrow(ordered_data)
  z_value <- qnorm(1 - (1 - confidence_level) / 2)
  
  # Initialize results storage
  cumulative_results <- data.frame(
    step = integer(),
    n_studies = integer(),
    estimate = numeric(),
    se = numeric(),
    ci_lower = numeric(),
    ci_upper = numeric(),
    z_value = numeric(),
    p_value = numeric(),
    tau2 = numeric(),
    i2 = numeric(),
    h2 = numeric(),
    q_statistic = numeric(),
    q_p_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Perform cumulative analysis
  for (i in minimum_studies:n_studies) {
    
    current_data <- ordered_data[1:i, ]
    
    # Run meta-analysis for current subset
    if (method == "random") {
      meta_result <- run_random_effects_ma(current_data)
    } else if (method == "fixed") {
      meta_result <- run_fixed_effects_ma(current_data)
    } else {
      meta_result <- run_bayesian_ma(current_data)
    }
    
    # Calculate additional statistics
    heterogeneity_stats <- calculate_heterogeneity_statistics(current_data, meta_result)
    
    # Store results
    cumulative_results <- rbind(cumulative_results, data.frame(
      step = i - minimum_studies + 1,
      n_studies = i,
      estimate = meta_result$estimate,
      se = meta_result$se,
      ci_lower = meta_result$estimate - z_value * meta_result$se,
      ci_upper = meta_result$estimate + z_value * meta_result$se,
      z_value = meta_result$estimate / meta_result$se,
      p_value = 2 * pnorm(abs(meta_result$estimate / meta_result$se), lower.tail = FALSE),
      tau2 = if(!is.null(meta_result$tau2)) meta_result$tau2 else 0,
      i2 = heterogeneity_stats$i2,
      h2 = heterogeneity_stats$h2,
      q_statistic = heterogeneity_stats$q_statistic,
      q_p_value = heterogeneity_stats$q_p_value,
      stringsAsFactors = FALSE
    ))
  }
  
  return(cumulative_results)
}

#' Run Random Effects Meta-Analysis
#'
#' Performs random effects meta-analysis for current data subset
#'
#' @param data Current data subset
#' @return Meta-analysis results
run_random_effects_ma <- function(data) {
  
  # Use metafor if available
  if (requireNamespace("metafor", quietly = TRUE)) {
    
    tryCatch({
      fit <- metafor::rma(yi = data$yi, sei = data$se, method = "REML")
      
      result <- list(
        estimate = as.numeric(fit$beta),
        se = fit$se,
        tau2 = fit$tau2,
        method = "metafor_random"
      )
      
    }, error = function(e) {
      # Fallback to basic implementation
      result <- run_basic_random_effects(data)
    })
    
  } else {
    # Basic random effects implementation
    result <- run_basic_random_effects(data)
  }
  
  return(result)
}

#' Run Fixed Effects Meta-Analysis
#'
#' Performs fixed effects meta-analysis for current data subset
#'
#' @param data Current data subset
#' @return Meta-analysis results
run_fixed_effects_ma <- function(data) {
  
  # Inverse variance weights
  weights <- 1 / data$se^2
  
  # Weighted mean
  estimate <- sum(weights * data$yi) / sum(weights)
  
  # Standard error
  se <- sqrt(1 / sum(weights))
  
  result <- list(
    estimate = estimate,
    se = se,
    tau2 = 0,
    method = "fixed_effects"
  )
  
  return(result)
}

#' Run Basic Random Effects Meta-Analysis
#'
#' Basic implementation of random effects meta-analysis
#'
#' @param data Current data subset
#' @return Meta-analysis results
run_basic_random_effects <- function(data) {
  
  # Start with fixed effects
  weights <- 1 / data$se^2
  fixed_estimate <- sum(weights * data$yi) / sum(weights)
  
  # Calculate Q statistic for tau2 estimation
  q_stat <- sum(weights * (data$yi - fixed_estimate)^2)
  df <- nrow(data) - 1
  
  # Estimate tau2 using method of moments
  if (df > 0) {
    c_val <- sum(weights) - sum(weights^2) / sum(weights)
    tau2 <- max(0, (q_stat - df) / c_val)
  } else {
    tau2 <- 0
  }
  
  # Random effects weights
  random_weights <- 1 / (data$se^2 + tau2)
  
  # Random effects estimate
  estimate <- sum(random_weights * data$yi) / sum(random_weights)
  se <- sqrt(1 / sum(random_weights))
  
  result <- list(
    estimate = estimate,
    se = se,
    tau2 = tau2,
    method = "basic_random"
  )
  
  return(result)
}

#' Run Bayesian Meta-Analysis
#'
#' Simple Bayesian implementation for cumulative analysis
#'
#' @param data Current data subset
#' @return Meta-analysis results
run_bayesian_ma <- function(data) {
  
  # Simple Bayesian approximation using normal priors
  # This is a simplified implementation
  
  # Weakly informative priors
  prior_mean <- 0
  prior_precision <- 0.01
  
  # Likelihood precisions
  precisions <- 1 / data$se^2
  
  # Posterior precision
  posterior_precision <- prior_precision + sum(precisions)
  
  # Posterior mean
  posterior_mean <- (prior_precision * prior_mean + sum(precisions * data$yi)) / posterior_precision
  
  # Posterior standard error
  posterior_se <- sqrt(1 / posterior_precision)
  
  result <- list(
    estimate = posterior_mean,
    se = posterior_se,
    tau2 = 0.1, # Fixed value for simplicity
    method = "bayesian_approx"
  )
  
  return(result)
}

#' Calculate Heterogeneity Statistics
#'
#' Calculates I-squared, H-squared, and Q statistics
#'
#' @param data Current data subset
#' @param meta_result Meta-analysis results
#' @return Heterogeneity statistics
calculate_heterogeneity_statistics <- function(data, meta_result) {
  
  if (nrow(data) < 2) {
    return(list(i2 = 0, h2 = 1, q_statistic = 0, q_p_value = 1))
  }
  
  # Calculate Q statistic
  weights <- 1 / data$se^2
  q_statistic <- sum(weights * (data$yi - meta_result$estimate)^2)
  df <- nrow(data) - 1
  q_p_value <- pchisq(q_statistic, df, lower.tail = FALSE)
  
  # Calculate I-squared
  i2 <- max(0, (q_statistic - df) / q_statistic) * 100
  
  # Calculate H-squared
  h2 <- q_statistic / df
  
  return(list(
    i2 = i2,
    h2 = h2,
    q_statistic = q_statistic,
    q_p_value = q_p_value
  ))
}

#' Calculate Stability Metrics
#'
#' Assesses how stable the cumulative estimates become
#'
#' @param cumulative_results Cumulative analysis results
#' @return Stability assessment
calculate_stability_metrics <- function(cumulative_results) {
  
  n_steps <- nrow(cumulative_results)
  
  if (n_steps < 3) {
    return(list(
      stability_achieved = FALSE,
      stability_point = NA,
      final_change = NA,
      trend_direction = "insufficient_data"
    ))
  }
  
  estimates <- cumulative_results$estimate
  
  # Calculate changes between consecutive estimates
  changes <- abs(diff(estimates))
  
  # Calculate relative changes
  relative_changes <- changes / abs(estimates[-1])
  
  # Define stability as relative change < 5% for last 3 steps
  if (n_steps >= 5) {
    last_changes <- tail(relative_changes, 3)
    stability_achieved <- all(last_changes < 0.05)
    
    if (stability_achieved) {
      # Find first point where stability was achieved
      stability_point <- NA
      for (i in 3:(n_steps-1)) {
        if (all(relative_changes[i:(i+2)] < 0.05)) {
          stability_point <- cumulative_results$n_studies[i+2]
          break
        }
      }
    } else {
      stability_point <- NA
    }
  } else {
    stability_achieved <- FALSE
    stability_point <- NA
  }
  
  # Final change magnitude
  final_change <- tail(relative_changes, 1)
  
  # Trend direction
  if (n_steps >= 3) {
    recent_trend <- estimates[n_steps] - estimates[n_steps-2]
    trend_direction <- if (recent_trend > 0) "increasing" else "decreasing"
  } else {
    trend_direction <- "insufficient_data"
  }
  
  return(list(
    stability_achieved = stability_achieved,
    stability_point = stability_point,
    final_change = final_change,
    trend_direction = trend_direction,
    relative_changes = relative_changes
  ))
}

#' Assess Evidence Sufficiency
#'
#' Determines if accumulated evidence is sufficient
#'
#' @param cumulative_results Cumulative analysis results
#' @return Sufficiency assessment
assess_evidence_sufficiency <- function(cumulative_results) {
  
  final_result <- tail(cumulative_results, 1)
  
  # Statistical significance
  is_significant <- final_result$p_value < 0.05
  
  # Confidence interval precision (width relative to estimate)
  ci_width <- final_result$ci_upper - final_result$ci_lower
  relative_precision <- ci_width / abs(final_result$estimate)
  is_precise <- relative_precision < 0.5
  
  # Sample size adequacy (heuristic)
  adequate_sample <- final_result$n_studies >= 10
  
  # Heterogeneity assessment
  low_heterogeneity <- final_result$i2 < 50
  
  # Overall sufficiency
  sufficient <- is_significant && is_precise && adequate_sample
  
  return(list(
    sufficient = sufficient,
    is_significant = is_significant,
    is_precise = is_precise,
    adequate_sample = adequate_sample,
    low_heterogeneity = low_heterogeneity,
    final_n_studies = final_result$n_studies,
    final_estimate = final_result$estimate,
    final_ci_width = ci_width,
    relative_precision = relative_precision
  ))
}

# S3 Methods for cbamm_cumulative class

#' Print Method for Cumulative Meta-Analysis
#'
#' @param x cbamm_cumulative object
#' @param ... Additional arguments
#' @export
print.cbamm_cumulative <- function(x, ...) {
  
  cat("Cumulative Meta-Analysis Results\n")
  cat("================================\n\n")
  
  cat("Ordering:", x$order_by, "(", x$order_direction, ")\n")
  cat("Method:", x$method, "\n")
  cat("Total studies:", max(x$results$n_studies), "\n")
  cat("Analysis steps:", nrow(x$results), "\n\n")
  
  # Final results
  final_result <- tail(x$results, 1)
  cat("Final Results:\n")
  cat("Estimate:", round(final_result$estimate, 4), "\n")
  cat("95% CI: [", round(final_result$ci_lower, 4), ", ", 
      round(final_result$ci_upper, 4), "]\n")
  cat("P-value:", format.pval(final_result$p_value), "\n")
  cat("I-squared:", round(final_result$i2, 1), "%\n\n")
  
  # Stability assessment
  if (x$stability$stability_achieved) {
    cat("Evidence stability: Achieved")
    if (!is.na(x$stability$stability_point)) {
      cat(" (at", x$stability$stability_point, "studies)")
    }
    cat("\n")
  } else {
    cat("Evidence stability: Not yet achieved\n")
  }
  
  # Sufficiency assessment
  cat("Evidence sufficiency:", if(x$sufficiency$sufficient) "Sufficient" else "Insufficient", "\n")
}

#' Summary Method for Cumulative Meta-Analysis
#'
#' @param object cbamm_cumulative object
#' @param ... Additional arguments
#' @export
summary.cbamm_cumulative <- function(object, ...) {
  
  structure(object, class = c("summary.cbamm_cumulative", class(object)))
}

#' Plot Method for Cumulative Meta-Analysis
#'
#' @param x cbamm_cumulative object
#' @param type Plot type: "estimate", "ci", "precision", "heterogeneity"
#' @param ... Additional arguments
#' @export
plot.cbamm_cumulative <- function(x, type = "estimate", ...) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required for plotting")
  }
  
  if (type == "estimate") {
    plot_cumulative_estimates(x, ...)
  } else if (type == "ci") {
    plot_cumulative_ci(x, ...)
  } else if (type == "precision") {
    plot_cumulative_precision(x, ...)
  } else if (type == "heterogeneity") {
    plot_cumulative_heterogeneity(x, ...)
  } else {
    stop("Plot type must be: estimate, ci, precision, or heterogeneity")
  }
}

#' Plot Cumulative Estimates
#'
#' @param x cbamm_cumulative object
#' @param ... Additional arguments
plot_cumulative_estimates <- function(x, ...) {
  
  data <- x$results
  
  p <- ggplot2::ggplot(data, ggplot2::aes(x = n_studies, y = estimate)) +
    ggplot2::geom_line(size = 1, color = "blue") +
    ggplot2::geom_point(size = 2, color = "blue") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_lower, ymax = ci_upper), 
                        alpha = 0.3, fill = "blue") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    ggplot2::labs(
      title = "Cumulative Meta-Analysis: Effect Size Evolution",
      subtitle = paste("Ordering by", x$order_by),
      x = "Number of Studies",
      y = "Cumulative Effect Size"
    ) +
    ggplot2::theme_minimal()
  
  print(p)
  return(p)
}

#' Plot Cumulative Confidence Intervals
#'
#' @param x cbamm_cumulative object
#' @param ... Additional arguments
plot_cumulative_ci <- function(x, ...) {
  
  data <- x$results
  data$ci_width <- data$ci_upper - data$ci_lower
  
  p <- ggplot2::ggplot(data, ggplot2::aes(x = n_studies)) +
    ggplot2::geom_line(ggplot2::aes(y = ci_width), size = 1, color = "red") +
    ggplot2::geom_point(ggplot2::aes(y = ci_width), size = 2, color = "red") +
    ggplot2::labs(
      title = "Cumulative Meta-Analysis: Confidence Interval Width",
      x = "Number of Studies",
      y = "95% CI Width"
    ) +
    ggplot2::theme_minimal()
  
  print(p)
  return(p)
}

#' Plot Cumulative Precision
#'
#' @param x cbamm_cumulative object
#' @param ... Additional arguments
plot_cumulative_precision <- function(x, ...) {
  
  data <- x$results
  data$precision <- 1 / data$se^2
  
  p <- ggplot2::ggplot(data, ggplot2::aes(x = n_studies, y = precision)) +
    ggplot2::geom_line(size = 1, color = "green") +
    ggplot2::geom_point(size = 2, color = "green") +
    ggplot2::labs(
      title = "Cumulative Meta-Analysis: Precision Evolution",
      x = "Number of Studies",
      y = "Precision (1/SE²)"
    ) +
    ggplot2::theme_minimal()
  
  print(p)
  return(p)
}

#' Plot Cumulative Heterogeneity
#'
#' @param x cbamm_cumulative object
#' @param ... Additional arguments
plot_cumulative_heterogeneity <- function(x, ...) {
  
  data <- x$results
  
  p <- ggplot2::ggplot(data, ggplot2::aes(x = n_studies, y = i2)) +
    ggplot2::geom_line(size = 1, color = "orange") +
    ggplot2::geom_point(size = 2, color = "orange") +
    ggplot2::geom_hline(yintercept = c(25, 50, 75), 
                       linetype = "dashed", color = "gray", alpha = 0.7) +
    ggplot2::labs(
      title = "Cumulative Meta-Analysis: Heterogeneity Evolution",
      x = "Number of Studies",
      y = "I² (%)"
    ) +
    ggplot2::theme_minimal()
  
  print(p)
  return(p)
}

