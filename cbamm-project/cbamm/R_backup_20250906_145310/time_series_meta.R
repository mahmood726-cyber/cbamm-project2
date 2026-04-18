# time_series_meta.R
# Time-Series and Longitudinal Meta-Analysis for CBAMM

#' Time-Series Meta-Analysis
#' 
#' @param data Data with multiple time points per study
#' @param time_var Name of time variable
#' @param effect_var Name of effect variable
#' @param study_var Name of study variable
#' @param method Analysis method
#' @export
time_series_meta <- function(data,
                            time_var = "time",
                            effect_var = "effect",
                            study_var = "study",
                            method = c("multilevel", "growth_curve")) {
  
  method <- match.arg(method)
  
  # Extract variables
  times <- data[[time_var]]
  effects <- data[[effect_var]]
  studies <- data[[study_var]]
  
  # Calculate time-specific effects
  unique_times <- sort(unique(times))
  time_effects <- data.frame(
    time = unique_times,
    effect = numeric(length(unique_times)),
    se = numeric(length(unique_times))
  )
  
  for (i in seq_along(unique_times)) {
    time_data <- effects[times == unique_times[i]]
    if (length(time_data) > 0) {
      time_effects$effect[i] <- mean(time_data, na.rm = TRUE)
      time_effects$se[i] <- sd(time_data, na.rm = TRUE) / sqrt(length(time_data))
    }
  }
  
  # Analyze trend
  trend_model <- lm(effect ~ time, data = time_effects)
  slope <- coef(trend_model)[2]
  slope_p <- summary(trend_model)$coefficients[2, 4]
  
  result <- list(
    time_effects = time_effects,
    trend = list(
      slope = slope,
      p_value = slope_p,
      direction = ifelse(slope > 0, "increasing", "decreasing"),
      significant = slope_p < 0.05
    ),
    method = method,
    data = data
  )
  
  class(result) <- c("cbamm_timeseries", class(result))
  return(result)
}

#' Print method for time-series meta-analysis
#' @export
print.cbamm_timeseries <- function(x, ...) {
  cat("\nTime-Series Meta-Analysis\n")
  cat("=========================\n")
  cat("Time points:", nrow(x$time_effects), "\n\n")
  
  cat("Temporal Trend:\n")
  cat("  Direction:", x$trend$direction, "\n")
  cat("  Slope:", round(x$trend$slope, 4), "\n")
  cat("  p-value:", format.pval(x$trend$p_value), "\n")
  cat("  Significant:", x$trend$significant, "\n")
  
  invisible(x)
}

