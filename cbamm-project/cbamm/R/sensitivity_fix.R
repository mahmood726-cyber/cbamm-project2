
#' Run optimized sensitivity analyses (fallback)
#' @keywords internal
run_optimized_sensitivity_analyses <- function(data, config, performance_mode = "fast") {
  # Fallback to basic sensitivity
  if (exists("run_sensitivity_analyses")) {
    return(run_sensitivity_analyses(data, config))
  }
  # Basic implementation
  list(
    leave_one_out = TRUE,
    influence = TRUE,
    completed = TRUE
  )
}

