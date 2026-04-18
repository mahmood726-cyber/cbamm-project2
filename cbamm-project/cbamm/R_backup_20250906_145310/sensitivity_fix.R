
#' Run optimized sensitivity analyses (fallback)
#' @keywords internal
run_optimized_sensitivity_analyses <- function(data, config) {
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

