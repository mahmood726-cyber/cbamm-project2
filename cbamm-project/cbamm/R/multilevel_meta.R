# multilevel_meta.R
# Advanced Multilevel Meta-Analysis Module for CBAMM Package

#' Run Multilevel Meta-Analysis (3-Level Models)
#'
#' Handles hierarchical data structures where effect sizes are nested 
#' (e.g., multiple outcomes within studies).
#'
#' @param data Data frame with study-level and effect-size data
#' @param outer_level Column name for outer grouping (e.g., "study_id")
#' @param inner_level Column name for inner nesting (e.g., "effect_id")
#' @param yi Column name for effect sizes
#' @param se Column name for standard errors
#' @param method Estimation method (default "REML")
#' @return Object of class "cbamm_multilevel"
#' @export
run_multilevel_meta <- function(data, 
                               outer_level = "study_id", 
                               inner_level = "effect_id",
                               yi = "yi", 
                               se = "se",
                               method = "REML") {
  
  if (!check_package("metafor", quietly = TRUE)) {
    stop("metafor package required for multilevel meta-analysis")
  }
  
  # Ensure inner level ID exists
  if (!inner_level %in% names(data)) {
    data[[inner_level]] <- 1:nrow(data)
  }
  
  # Fit 3-level model
  # random = ~ 1 | outer / inner
  random_formula <- as.formula(paste("~ 1 |", outer_level, "/", inner_level))
  
  fit <- try(metafor::rma.mv(yi = data[[yi]], 
                            V = data[[se]]^2, 
                            random = random_formula,
                            data = data,
                            method = method), silent = TRUE)
  
  if (inherits(fit, "try-error")) {
    stop("Multilevel model fitting failed: ", as.character(fit))
  }
  
  # Calculate distribution of variance across levels
  total_var <- sum(fit$sigma2)
  var_dist <- if(total_var > 0) fit$sigma2 / total_var * 100 else rep(0, length(fit$sigma2))
  names(var_dist) <- c(paste("Level 3 (", outer_level, ")", sep=""), 
                      paste("Level 2 (", inner_level, ")", sep=""))
  
  result <- list(
    model = fit,
    variance_distribution = var_dist,
    data = data,
    levels = list(outer = outer_level, inner = inner_level)
  )
  
  class(result) <- "cbamm_multilevel"
  return(result)
}

#' Print method for cbamm_multilevel
#' @export
print.cbamm_multilevel <- function(x, ...) {
  cat("Multilevel (3-Level) Meta-Analysis Results
")
  cat("------------------------------------------
")
  cat("Outer Level:", x$levels$outer, "
")
  cat("Inner Level:", x$levels$inner, "

")
  
  cat("Variance Distribution:
")
  for(i in seq_along(x$variance_distribution)) {
    cat(sprintf("  %s: %.1f%%
", names(x$variance_distribution)[i], x$variance_distribution[i]))
  }
  cat("
")
  
  print(summary(x$model))
}
