#' Generate CBAMM Plots
#'
#' @param results CBAMM results object
#' @param config CBAMM configuration
#'
#' @return List of plots
#' @keywords internal
generate_cbamm_plots <- function(results, config) {
  plots <- list()
  
  # Helper function to safely add plots
  add_plot <- function(lst, key, plot_obj) {
    if (!is.null(plot_obj)) lst[[key]] <- plot_obj
    lst
  }
  
  # Forest plot
  if (!is.null(results$data_summary)) {
    plots <- add_plot(plots, "forest", 
      try(.create_forest_plot(results, config), silent = TRUE) |>
        (\(x) if (inherits(x, "try-error")) NULL else x)())
  }
  
  # Funnel plot  
  if (!is.null(results$meta_results)) {
    plots <- add_plot(plots, "funnel",
      try(.create_funnel_plot(results$meta_results, config), silent = TRUE) |>
        (\(x) if (inherits(x, "try-error")) NULL else x)())
  }
  
  # Multiverse plot
  if (!is.null(results$multiverse)) {
    plots <- add_plot(plots, "multiverse",
      try(.create_multiverse_plot(results$multiverse, config), silent = TRUE) |>
        (\(x) if (inherits(x, "try-error")) NULL else x)())
  }
  
  return(plots)
}

#' Create Forest Plot
#' @keywords internal  
.create_forest_plot <- function(results, config) {
  # This is a simplified version - full implementation would extract data from results
  # For now, return a placeholder
  p <- ggplot2::ggplot() + 
    ggplot2::labs(title = "Forest Plot", subtitle = "Implementation in progress") +
    ggplot2::theme_minimal()
  
  if (config$output$interactive && check_package("plotly", quietly = TRUE)) {
    return(plotly::ggplotly(p))
  }
  
  p
}

#' Create Funnel Plot
#' @keywords internal
.create_funnel_plot <- function(meta_results, config) {
  fit <- meta_results$transport
  if (is.null(fit)) return(NULL)
  
  plot_data <- data.frame(
    yi = as.numeric(fit$yi), 
    sei = sqrt(as.numeric(fit$vi))
  )
  plot_data <- dplyr::filter(plot_data, is.finite(yi) & is.finite(sei))
  
  if (nrow(plot_data) == 0) return(NULL)
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = yi, y = sei)) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::geom_vline(xintercept = as.numeric(coef(fit)), 
                       linetype = "dashed", color = "red") +
    gg