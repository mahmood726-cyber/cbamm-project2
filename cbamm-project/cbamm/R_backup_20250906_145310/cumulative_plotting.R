
#' Create Comprehensive Cumulative Analysis Dashboard
#'
#' Creates a four-panel dashboard showing effect evolution, precision, 
#' statistical significance, and heterogeneity over cumulative steps.
#'
#' @param x cbamm_cumulative object
#' @param save_plot Logical, whether to save the plot to file
#' @param filename Character, filename for saving (if save_plot = TRUE)
#' @param width Numeric, plot width in inches (default: 12)
#' @param height Numeric, plot height in inches (default: 10)
#' @return Creates dashboard plot (invisible return)
#' @export
#' @examples
#' \dontrun{
#' result <- cumulative_meta_analysis(cumulative_example_data, order_by = "year")
#' create_cumulative_dashboard(result)
#' create_cumulative_dashboard(result, save_plot = TRUE, filename = "dashboard.png")
#' }
create_cumulative_dashboard <- function(x, save_plot = FALSE, filename = "cumulative_dashboard.png", 
                                       width = 12, height = 10) {
  
  if (!inherits(x, "cbamm_cumulative")) {
    stop("Input must be a cbamm_cumulative object")
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required for dashboard plotting")
  }
  
  data <- x$results
  
  # Plot 1: Effect size evolution with stability indicators
  p1 <- ggplot2::ggplot(data, ggplot2::aes(x = n_studies, y = estimate)) +
    ggplot2::geom_line(linewidth = 1.2, color = "blue", alpha = 0.8) +
    ggplot2::geom_point(size = 2.5, color = "blue", alpha = 0.7) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_lower, ymax = ci_upper), 
                        alpha = 0.3, fill = "blue") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
    {if(!is.na(x$stability$stability_point)) 
      ggplot2::geom_vline(xintercept = x$stability$stability_point, 
                         linetype = "dotted", color = "green", linewidth = 1)} +
    ggplot2::labs(
      title = "Effect Size Evolution",
      x = "Number of Studies",
      y = "Cumulative Effect Size",
      subtitle = if(!is.na(x$stability$stability_point)) 
        paste("Stability achieved at", x$stability$stability_point, "studies") else "Stability not yet achieved"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  
  # Plot 2: Precision evolution
  data$precision <- 1 / data$se^2
  p2 <- ggplot2::ggplot(data, ggplot2::aes(x = n_studies, y = precision)) +
    ggplot2::geom_line(linewidth = 1.2, color = "green", alpha = 0.8) +
    ggplot2::geom_point(size = 2.5, color = "green", alpha = 0.7) +
    ggplot2::labs(
      title = "Precision Evolution",
      x = "Number of Studies",
      y = "Precision (1/SE²)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  
  # Plot 3: Statistical significance tracking
  data$significant <- data$p_value < 0.05
  p3 <- ggplot2::ggplot(data, ggplot2::aes(x = n_studies, y = -log10(p_value))) +
    ggplot2::geom_line(linewidth = 1.2, color = "orange", alpha = 0.8) +
    ggplot2::geom_point(ggplot2::aes(color = significant), size = 2.5, alpha = 0.7) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    ggplot2::scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "orange")) +
    ggplot2::labs(
      title = "Statistical Significance Evolution",
      x = "Number of Studies", 
      y = "-log10(p-value)",
      color = "Significant"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  
  # Plot 4: Heterogeneity with interpretation bands
  p4 <- ggplot2::ggplot(data, ggplot2::aes(x = n_studies, y = i2)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = 25), fill = "green", alpha = 0.1) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = 25, ymax = 50), fill = "yellow", alpha = 0.1) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = 50, ymax = 75), fill = "orange", alpha = 0.1) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = 75, ymax = 100), fill = "red", alpha = 0.1) +
    ggplot2::geom_line(linewidth = 1.2, color = "darkred", alpha = 0.8) +
    ggplot2::geom_point(size = 2.5, color = "darkred", alpha = 0.7) +
    ggplot2::geom_hline(yintercept = c(25, 50, 75), 
                       linetype = "dashed", color = "gray", alpha = 0.7) +
    ggplot2::labs(
      title = "Heterogeneity Evolution",
      x = "Number of Studies",
      y = "I² (%)",
      subtitle = "Green: Low, Yellow: Moderate, Orange: Substantial, Red: Considerable"
    ) +
    ggplot2::ylim(0, 100) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  
  # Combine plots using available packages
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    combined_plot <- gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
  } else if (requireNamespace("patchwork", quietly = TRUE)) {
    combined_plot <- (p1 + p2) / (p3 + p4)
    print(combined_plot)
  } else {
    # Fallback: print individual plots
    print(p1)
    print(p2)
    print(p3)
    print(p4)
  }
  
  if (save_plot) {
    if (requireNamespace("gridExtra", quietly = TRUE)) {
      ggplot2::ggsave(filename, combined_plot, width = width, height = height, dpi = 300)
    } else {
      cat("Note: Saving individual plots since gridExtra/patchwork not available\n")
      ggplot2::ggsave(gsub("\\.([^.]+)$", "_p1.\\1", filename), p1, width = width/2, height = height/2, dpi = 300)
      ggplot2::ggsave(gsub("\\.([^.]+)$", "_p2.\\1", filename), p2, width = width/2, height = height/2, dpi = 300)
      ggplot2::ggsave(gsub("\\.([^.]+)$", "_p3.\\1", filename), p3, width = width/2, height = height/2, dpi = 300)
      ggplot2::ggsave(gsub("\\.([^.]+)$", "_p4.\\1", filename), p4, width = width/2, height = height/2, dpi = 300)
    }
    cat("Dashboard saved as:", filename, "\n")
  }
  
  return(invisible())
}

#' Plot Evidence Accumulation Trajectory
#'
#' Creates an advanced visualization showing how evidence accumulates with
#' line thickness representing cumulative weight and colors showing change magnitude.
#'
#' @param x cbamm_cumulative object
#' @param highlight_studies Numeric vector of study numbers to highlight
#' @return ggplot object
#' @export
#' @examples
#' \dontrun{
#' result <- cumulative_meta_analysis(cumulative_example_data, order_by = "year")
#' plot_evidence_trajectory(result)
#' plot_evidence_trajectory(result, highlight_studies = c(10, 15, 20))
#' }
plot_evidence_trajectory <- function(x, highlight_studies = NULL) {
  
  if (!inherits(x, "cbamm_cumulative")) {
    stop("Input must be a cbamm_cumulative object")
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required for trajectory plotting")
  }
  
  data <- x$results
  data$trajectory_change <- c(NA, abs(diff(data$estimate)))
  data$cumulative_weight <- cumsum(1/data$se^2)
  
  # Create trajectory plot
  p <- ggplot2::ggplot(data, ggplot2::aes(x = n_studies)) +
    # Confidence interval ribbon
    ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_lower, ymax = ci_upper), 
                        alpha = 0.2, fill = "blue") +
    # Effect size line with varying thickness based on precision
    ggplot2::geom_line(ggplot2::aes(y = estimate, linewidth = cumulative_weight), 
                      color = "blue", alpha = 0.8) +
    # Points colored by magnitude of change
    ggplot2::geom_point(ggplot2::aes(y = estimate, color = trajectory_change), 
                       size = 3, alpha = 0.8) +
    # Highlight specific studies if requested
    {if(!is.null(highlight_studies)) 
      ggplot2::geom_point(data = data[data$n_studies %in% highlight_studies, ],
                         ggplot2::aes(y = estimate), 
                         color = "red", size = 4, shape = 1, stroke = 2)} +
    # Reference line at zero
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
    # Styling
    ggplot2::scale_color_gradient(low = "lightblue", high = "darkblue", na.value = "gray") +
    ggplot2::scale_linewidth_continuous(range = c(0.5, 2)) +
    ggplot2::labs(
      title = "Evidence Accumulation Trajectory",
      subtitle = paste("Ordering by", x$order_by, "-", x$order_direction),
      x = "Number of Studies",
      y = "Cumulative Effect Size",
      color = "Change\nMagnitude",
      linewidth = "Cumulative\nWeight"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )
  
  print(p)
  return(p)
}

#' Plot Stability Assessment
#'
#' Visualizes how the stability of cumulative estimates changes over time,
#' showing relative changes between consecutive estimates.
#'
#' @param x cbamm_cumulative object
#' @return ggplot object
#' @export
#' @examples
#' \dontrun{
#' result <- cumulative_meta_analysis(cumulative_example_data, order_by = "year")
#' plot_stability_assessment(result)
#' }
plot_stability_assessment <- function(x) {
  
  if (!inherits(x, "cbamm_cumulative")) {
    stop("Input must be a cbamm_cumulative object")
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required for stability plotting")
  }
  
  data <- x$results
  
  # Calculate rolling stability metrics
  if (length(x$stability$relative_changes) > 0) {
    stability_data <- data.frame(
      n_studies = data$n_studies[-1],
      relative_change = x$stability$relative_changes,
      stable = x$stability$relative_changes < 0.05
    )
    
    p <- ggplot2::ggplot(stability_data, ggplot2::aes(x = n_studies, y = relative_change)) +
      ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 1) +
      ggplot2::geom_line(linewidth = 1, color = "blue", alpha = 0.8) +
      ggplot2::geom_point(ggplot2::aes(color = stable), size = 3, alpha = 0.8) +
      {if(!is.na(x$stability$stability_point)) 
        ggplot2::geom_vline(xintercept = x$stability$stability_point, 
                           linetype = "dotted", color = "green", linewidth = 1)} +
      ggplot2::scale_color_manual(values = c("FALSE" = "red", "TRUE" = "green")) +
      ggplot2::scale_y_continuous(labels = function(x) paste0(round(x*100, 1), "%")) +
      ggplot2::labs(
        title = "Stability Assessment",
        subtitle = if(x$stability$stability_achieved) 
          paste("Stability achieved at", x$stability$stability_point, "studies") else
          "Stability not yet achieved",
        x = "Number of Studies",
        y = "Relative Change from Previous Step (%)",
        color = "Stable\n(< 5%)"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
    
  } else {
    p <- ggplot2::ggplot() + 
      ggplot2::geom_text(ggplot2::aes(x = 0.5, y = 0.5, label = "Insufficient data for stability assessment"),
                        size = 6) +
      ggplot2::theme_void()
  }
  
  print(p)
  return(p)
}

#' Create Evolving Forest Plot
#'
#' Creates a forest plot showing cumulative effect sizes at selected analysis steps,
#' useful for visualizing how estimates evolve over time.
#'
#' @param x cbamm_cumulative object
#' @param steps_to_show Numeric vector of steps to display (default: every 5th step)
#' @return ggplot object
#' @export
#' @examples
#' \dontrun{
#' result <- cumulative_meta_analysis(cumulative_example_data, order_by = "year")
#' plot_forest_evolution(result)
#' plot_forest_evolution(result, steps_to_show = c(1, 5, 10, 15, 20, 24))
#' }
plot_forest_evolution <- function(x, steps_to_show = NULL) {
  
  if (!inherits(x, "cbamm_cumulative")) {
    stop("Input must be a cbamm_cumulative object")
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required for forest evolution plotting")
  }
  
  data <- x$results
  
  if (is.null(steps_to_show)) {
    steps_to_show <- unique(c(seq(1, nrow(data), by = 5), nrow(data)))
  }
  
  forest_data <- data[data$step %in% steps_to_show, ]
  forest_data$step_label <- paste("Step", forest_data$step, "(n =", forest_data$n_studies, ")")
  
  p <- ggplot2::ggplot(forest_data, ggplot2::aes(y = factor(step_label, levels = rev(step_label)))) +
    ggplot2::geom_point(ggplot2::aes(x = estimate), size = 3, color = "blue") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = ci_lower, xmax = ci_upper), 
                           height = 0.2, color = "blue", alpha = 0.7) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
    ggplot2::labs(
      title = "Forest Plot Evolution",
      subtitle = "Cumulative effect sizes at selected steps",
      x = "Effect Size",
      y = "Analysis Step"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(size = 8)
    )
  
  print(p)
  return(p)
}

