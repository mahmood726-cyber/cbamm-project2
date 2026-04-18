
#' Bayesian Meta-Analysis Wrapper
#' @export
bayesian_meta_analysis_fixed <- function(data, ...) {
  # Ensure data has required columns
  if (!all(c("effect", "variance") %in% names(data))) {
    # Try to create from yi/se or effect_size/standard_error
    if (all(c("yi", "se") %in% names(data))) {
      data$effect <- data$yi
      data$variance <- data$se^2
    } else if (all(c("effect_size", "standard_error") %in% names(data))) {
      data$effect <- data$effect_size
      data$variance <- data$standard_error^2
    }
  }
  
  # Remove any non-numeric values
  data <- data[!is.na(data$effect) & !is.na(data$variance), ]
  
  # Call original function
  return(bayesian_meta_analysis(data, ...))
}

