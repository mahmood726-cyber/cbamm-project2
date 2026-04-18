optimize_config_for_performance <- function(config, performance_mode = "standard", data = NULL) {
  performance_mode <- match.arg(
    performance_mode,
    c("standard", "fast", "balanced", "comprehensive")
  )

  config <- if (is.null(config)) cbamm_config() else config
  n_studies <- if (is.null(data) || is.null(nrow(data))) 0L else nrow(data)

  if (performance_mode == "standard") {
    config$performance_mode <- performance_mode
    return(config)
  }

  config <- create_optimized_config(
    base_config = config,
    performance_mode = performance_mode,
    n_studies = n_studies
  )

  if (is.null(config$methods)) config$methods <- list()
  if (is.null(config$output)) config$output <- list()
  if (is.null(config$bayesian)) config$bayesian <- list()
  if (is.null(config$sensitivity)) config$sensitivity <- list()

  if (isTRUE(config$skip_bayesian)) {
    config$methods$bayesian <- FALSE
  }

  if (isTRUE(config$skip_complex_plots)) {
    config$output$plots <- FALSE
  }

  if (!is.null(config$bayesian_method) && !identical(config$bayesian_method, "none")) {
    config$bayesian$method <- config$bayesian_method
  }

  if (isTRUE(config$minimal_sensitivity)) {
    config$methods$pet_peese <- FALSE
    config$methods$missing_studies <- FALSE
    config$methods$multiverse <- FALSE
    config$sensitivity$publication_bias_methods <- character(0)
  }

  config
}
