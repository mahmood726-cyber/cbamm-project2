# living_systematic_review.R
# Living Systematic Review and Sequential Meta-Analysis for CBAMM Package

#' Living Systematic Review Framework
#' 
#' @param initial_data Initial data frame
#' @param effect_size Column name for effect sizes
#' @param variance Column name for variances
#' @param study_id Column name for study IDs
#' @param se Whether variance column contains standard errors
#' @param method Sequential method
#' @param boundaries Boundary type
#' @param alpha Significance level
#' @param power Statistical power
#' @param minimum_studies Minimum number of studies
#' @param target_effect Target effect size
#' @param ... Additional arguments
#' @return Object of class cbamm_living
#' @export
living_systematic_review <- function(initial_data, 
                                    effect_size = "effect",
                                    variance = "variance",
                                    study_id = "study_id",
                                    se = FALSE,
                                    method = c("sequential", "cumulative"),
                                    boundaries = c("obrien-fleming", "pocock"),
                                    alpha = 0.05,
                                    power = 0.80,
                                    minimum_studies = 3,
                                    target_effect = NULL,
                                    ...) {
  
  method <- match.arg(method)
  boundaries <- match.arg(boundaries)
  
  # Extract data
  effects <- initial_data[[effect_size]]
  
  if (se) {
    ses <- initial_data[[variance]]
    variances <- ses^2
  } else {
    variances <- initial_data[[variance]]
    ses <- sqrt(variances)
  }
  
  if (study_id %in% names(initial_data)) {
    studies <- initial_data[[study_id]]
  } else {
    studies <- seq_len(nrow(initial_data))
  }
  
  n_studies <- length(studies)
  
  # Calculate cumulative Z-scores
  cumulative_z <- numeric(n_studies)
  cumulative_effects <- numeric(n_studies)
  
  for (i in 1:n_studies) {
    weights <- 1 / variances[1:i]
    pooled_effect <- sum(effects[1:i] * weights) / sum(weights)
    pooled_se <- sqrt(1 / sum(weights))
    cumulative_z[i] <- pooled_effect / pooled_se
    cumulative_effects[i] <- pooled_effect
  }
  
  # Calculate boundaries
  alpha_one_sided <- alpha / 2
  
  if (boundaries == "obrien-fleming") {
    info_frac <- (1:n_studies) / n_studies
    upper <- qnorm(1 - alpha_one_sided) / sqrt(info_frac)
    lower <- -upper
  } else {  # Pocock
    upper <- rep(qnorm(1 - alpha_one_sided), n_studies)
    lower <- -upper
  }
  
  # Check for boundary crossing
  stopped <- any(cumulative_z > upper | cumulative_z < lower)
  
  # Final meta-analysis
  weights <- 1 / ses^2
  pooled_effect <- sum(effects * weights) / sum(weights)
  pooled_se <- sqrt(1 / sum(weights))
  
  result <- list(
    results = list(
      effect = pooled_effect,
      se = pooled_se,
      ci_lower = pooled_effect - 1.96 * pooled_se,
      ci_upper = pooled_effect + 1.96 * pooled_se,
      p_value = 2 * pnorm(-abs(pooled_effect / pooled_se))
    ),
    sequential = list(
      cumulative_z = cumulative_z,
      cumulative_effects = cumulative_effects,
      boundaries = list(upper = upper, lower = lower),
      stopped = stopped
    ),
    settings = list(
      method = method,
      boundaries = boundaries,
      alpha = alpha,
      power = power
    ),
    data = initial_data
  )
  
  class(result) <- c("cbamm_living", class(result))
  return(result)
}

#' Update Living Review
#' 
#' @param living_review Existing living review object
#' @param new_data New data to add
#' @param ... Additional arguments
#' @return Updated living review object
#' @export
update_living_review <- function(living_review, new_data, ...) {
  
  # Combine data
  combined_data <- rbind(living_review$data, new_data)
  
  # Get column names from original data
  effect_col <- names(living_review$data)[2]  # Assume second column
  variance_col <- names(living_review$data)[3]  # Assume third column
  study_col <- names(living_review$data)[1]  # Assume first column
  
  # Re-run analysis
  updated_review <- living_systematic_review(
    initial_data = combined_data,
    effect_size = effect_col,
    variance = variance_col,
    study_id = study_col,
    method = living_review$settings$method,
    boundaries = living_review$settings$boundaries,
    alpha = living_review$settings$alpha,
    power = living_review$settings$power,
    ...
  )
  
  return(updated_review)
}

#' Simulate Sequential Data
#' 
#' @param n_studies Number of studies
#' @param true_effect True effect size
#' @param start_id Starting study ID
#' @return Data frame with sequential data
#' @export
simulate_sequential_data <- function(n_studies = 10, 
                                    true_effect = 0.3,
                                    start_id = 1) {
  
  study_ids <- start_id:(start_id + n_studies - 1)
  
  # Generate effects with some variation
  effects <- rnorm(n_studies, true_effect, 0.1)
  
  # Generate standard errors
  ses <- runif(n_studies, 0.05, 0.2)
  
  return(data.frame(
    study_id = paste0("Study_", study_ids),
    effect = effects,
    variance = ses^2,
    se = ses
  ))
}

#' Print method for living review
#' @export
print.cbamm_living <- function(x, ...) {
  cat("Living Systematic Review\n")
  cat("-----------------------\n")
  cat("Method:", x$settings$method, "\n")
  cat("Boundaries:", x$settings$boundaries, "\n")
  cat("Studies:", nrow(x$data), "\n")
  cat("Stopped:", x$sequential$stopped, "\n")
  cat("Current effect:", round(x$results$effect, 3), "\n")
  cat("Current Z:", round(tail(x$sequential$cumulative_z, 1), 3), "\n")
  invisible(x)
}

#' Summary method for living review
#' @export
summary.cbamm_living <- function(object, ...) {
  print(object)
}
