#' Utility Functions for CBAMM Package.
#' 
#' Internal utility functions used throughout CBAMM

#' Check if Package is Available.
#' 
#' Safely check if a package is available and can be loaded
#' 
#' @param pkg Package name
#' @param quietly Whether to suppress messages
#' @return Logical indicating if package is available
#' @keywords internal
check_package <- function(pkg, quietly = TRUE) {
  requireNamespace(pkg, quietly = quietly)
}

#' Safe Package Loading with Fallback.
#' 
#' Attempt to use a package function with fallback behavior
#' 
#' @param pkg Package name
#' @param fun_call Function call to attempt
#' @param fallback_fun Fallback function if package unavailable
#' @param fallback_msg Message to show when using fallback
#' @keywords internal
with_fallback <- function(pkg, fun_call, fallback_fun = NULL, fallback_msg = NULL) {
  
  if (check_package(pkg, quietly = TRUE)) {
    return(fun_call)
  } else {
    if (!is.null(fallback_msg)) {
      message(fallback_msg)
    }
    if (!is.null(fallback_fun)) {
      return(fallback_fun())
    } else {
      stop("Package '", pkg, "' is required but not available. Please install it.")
    }
  }
}

#' Transform effect size metrics. Metrics.
#' 
#' Convert between different effect size metrics
#' 
#' @param yi Effect sizes
#' @param from Source metric
#' @param to Target metric  
#' @param ... Additional parameters for conversion
#' @return Converted effect sizes
convert_effect_size <- function(yi, from = "Cohen_d", to = "Hedges_g", ...) {
  
  if (from == to) return(yi)
  
  # Cohen's d to Hedges' g
  if (from == "Cohen_d" && to == "Hedges_g") {
    n <- list(...)$n
    if (is.null(n)) stop("Sample size 'n' required for Cohen's d to Hedges' g conversion")
    j <- 1 - (3 / (4 * n - 9))
    return(yi * j)
  }
  
  # Log odds ratio to Cohen's d
  if (from == "log_OR" && to == "Cohen_d") {
    return(yi * pi / sqrt(3))
  }
  
  # Add more conversions as needed
  stop("Conversion from '", from, "' to '", to, "' not implemented")
}

#' Compute Variance from Standard Error
#' 
#' @param se Standard errors
#' @return Variances
se_to_var <- function(se) {
  se^2
}

#' Compute Standard Error from Variance
#' 
#' @param vi Variances
#' @return Standard errors  
var_to_se <- function(vi) {
  sqrt(vi)
}

#' Robust Variance Calculation.
#' 
#' Calculate robust variance estimates for clustered data
#' 
#' @param data Data frame
#' @param cluster_var Clustering variable
#' @param method Robust variance method
#' @return Robust variance estimates
#' @keywords internal
compute_robust_variance <- function(data, cluster_var = NULL, method = "CR2") {
  
  if (is.null(cluster_var)) {
    return(data$se^2)  # No clustering
  }
  
  # Use clubSandwich if available
  with_fallback(
    pkg = "clubSandwich",
    fun_call = {
      # Implementation would go here
      # This is a placeholder for the actual robust variance calculation
      data$se^2
    },
    fallback_fun = function() data$se^2,
    fallback_msg = "clubSandwich not available, using standard variances"
  )
}

#' Format p-values for display.-values.
#' 
#' Format p-values for consistent display
#' 
#' @param p P-values
#' @param digits Number of decimal places
#' @return Formatted p-values
format_p <- function(p, digits = 3) {
  
  if (any(is.na(p))) {
    p[is.na(p)] <- NA
  }
  
  formatted <- ifelse(p < 0.001, "< 0.001",
                     ifelse(p < 0.01, sprintf("%.3f", p),
                           sprintf(paste0("%.", digits, "f"), p)))
  
  return(formatted)
}

#' Format Confidence Intervals
#' 
#' @param estimate Point estimate
#' @param lower Lower bound
#' @param upper Upper bound  
#' @param digits Number of decimal places
#' @return Formatted CI string
format_ci <- function(estimate, lower, upper, digits = 2) {
  sprintf("%.*f [%.*f, %.*f]", digits, estimate, digits, lower, digits, upper)
}

#' Calculate I-squared Statistic
#' 
#' @param Q Q statistic
#' @param df Degrees of freedom
#' @return I-squared value
calculate_i_squared <- function(Q, df) {
  i2 <- ((Q - df) / Q) * 100
  pmax(0, i2)  # I-squared cannot be negative
}

#' Generate Study Weights Summary
#' 
#' @param weights Study weights
#' @param study_ids Study identifiers
#' @return Data frame with weight summary
summarize_weights <- function(weights, study_ids = NULL) {
  
  if (is.null(study_ids)) {
    study_ids <- paste0("Study_", seq_along(weights))
  }
  
  weight_df <- data.frame(
    study_id = study_ids,
    weight = weights,
    weight_pct = (weights / sum(weights)) * 100,
    stringsAsFactors = FALSE
  )
  
  weight_df[order(weight_df$weight, decreasing = TRUE), ]
}

#' Protected division operation..
#' 
#' Division with handling for zero denominators
#' 
#' @param numerator Numerator
#' @param denominator Denominator
#' @param na_value Value to return when denominator is zero
#' @return Result of division
#' @keywords internal
safe_divide <- function(numerator, denominator, na_value = NA) {
  ifelse(denominator == 0 | is.na(denominator), na_value, numerator / denominator)
}

#' Detect Outliers Using IQR Method
#' 
#' @param x Numeric vector
#' @param k Multiplier for IQR (default 1.5)
#' @return Logical vector indicating outliers
detect_outliers_iqr <- function(x, k = 1.5) {
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  
  lower_bound <- q1 - k * iqr
  upper_bound <- q3 + k * iqr
  
  x < lower_bound | x > upper_bound
}

#' Generate summary tables..
#' 
#' Create a formatted summary table for results
#' 
#' @param results_list List of results to summarize
#' @param round_digits Number of digits to round to
#' @return Data frame with formatted results
create_summary_table <- function(results_list, round_digits = 3) {
  
  # Extract key statistics from different result types
  summary_df <- data.frame(
    Method = character(),
    Estimate = numeric(),
    SE = numeric(),
    Lower_CI = numeric(),
    Upper_CI = numeric(),
    P_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Add results from different methods
  for (method_name in names(results_list)) {
    result <- results_list[[method_name]]
    
    if (!is.null(result) && is.list(result)) {
      summary_df <- rbind(summary_df, data.frame(
        Method = method_name,
        Estimate = round(result$beta %||% result$estimate %||% NA, round_digits),
        SE = round(result$se %||% NA, round_digits),
        Lower_CI = round(result$ci.lb %||% result$lower %||% NA, round_digits),
        Upper_CI = round(result$ci.ub %||% result$upper %||% NA, round_digits),
        P_value = round(result$pval %||% result$p_value %||% NA, round_digits),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(summary_df)
}

#' Null Coalescing Operator.
#' 
#' Returns first non-NULL value
#' 
#' @param x First value
#' @param y Second value  
#' @return First non-NULL value
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
