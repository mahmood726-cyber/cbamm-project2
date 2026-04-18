
optimize_config_for_performance <- function(config, performance_mode = "standard", data = NULL) {
  if (is.null(config)) {
    config <- list(
      method = "fixed",
      conf_level = 0.95,
      continuity_correction = "auto",
      show_forest = TRUE,
      show_funnel = TRUE
    )
  }
  return(config)
}