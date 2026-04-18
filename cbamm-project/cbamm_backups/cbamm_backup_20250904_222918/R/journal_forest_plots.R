
# Journal and Software-Specific Forest Plots

#' NEJM Style Forest Plot
#' @param data Data frame with study data
#' @param pooled_est Pooled estimate (optional)
#' @export
forest_plot_nejm <- function(data, pooled_est = NULL) {
  plot_data <- data
  plot_data$study_order <- nrow(plot_data):1
  
  if (is.null(pooled_est)) {
    pooled_est <- list(or = 0.55, ci_lower = 0.42, ci_upper = 0.73)
  }
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(y = study_order)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "solid", color = "black", size = 0.5) +
    ggplot2::geom_point(ggplot2::aes(x = or), size = 2.5, shape = 15, color = "black") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = ci_lower, xmax = ci_upper),
                           height = 0.15, color = "black", size = 0.5) +
    ggplot2::geom_point(ggplot2::aes(x = pooled_est$or, y = 0),
                       size = 4, shape = 18, color = "black") +
    ggplot2::scale_x_log10(limits = c(0.1, 10), 
                          breaks = c(0.1, 0.5, 1, 2, 5, 10)) +
    ggplot2::labs(title = "Effect of Treatment vs Control", x = "Odds Ratio (log scale)") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"))
  
  print(p)
  return(p)
}

#' Lancet Style Forest Plot
#' @param data Data frame with study data
#' @param pooled_est Pooled estimate (optional)
#' @export
forest_plot_lancet <- function(data, pooled_est = NULL) {
  plot_data <- data
  plot_data$study_order <- nrow(plot_data):1
  
  if (is.null(pooled_est)) {
    pooled_est <- list(or = 0.55, ci_lower = 0.42, ci_upper = 0.73)
  }
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(y = study_order)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "solid", color = "black", size = 0.8) +
    ggplot2::geom_point(ggplot2::aes(x = or, size = weight), 
                       shape = 15, color = "#2E75B6", alpha = 0.8) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = ci_lower, xmax = ci_upper),
                           height = 0.2, color = "#2E75B6", size = 0.6) +
    ggplot2::geom_point(ggplot2::aes(x = pooled_est$or, y = 0),
                       size = 6, shape = 18, color = "black") +
    ggplot2::scale_size_continuous(range = c(2, 8), guide = "none") +
    ggplot2::scale_x_log10(limits = c(0.1, 5), breaks = c(0.1, 0.2, 0.5, 1, 2, 5)) +
    ggplot2::labs(title = "Meta-analysis of Treatment Effect", x = "Odds ratio (95% CI)") +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"))
  
  print(p)
  return(p)
}

#' RevMan Style Forest Plot
#' @param data Data frame with study data
#' @param pooled_est Pooled estimate (optional)
#' @export
forest_plot_revman <- function(data, pooled_est = NULL) {
  plot_data <- data
  plot_data$study_order <- nrow(plot_data):1
  
  if (is.null(pooled_est)) {
    pooled_est <- list(or = 0.55, ci_lower = 0.42, ci_upper = 0.73)
  }
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(y = study_order)) +
    ggplot2::geom_rect(ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
                      fill = "white", color = "black", size = 1) +
    ggplot2::geom_vline(xintercept = 1, linetype = "solid", color = "black", size = 0.5) +
    ggplot2::geom_point(ggplot2::aes(x = or, size = weight), 
                       shape = 15, color = "black") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = ci_lower, xmax = ci_upper),
                           height = 0.15, color = "black", size = 0.5) +
    ggplot2::geom_point(ggplot2::aes(x = pooled_est$or, y = 0),
                       size = 8, shape = 18, color = "black", fill = "#FFD700") +
    ggplot2::scale_size_continuous(range = c(1.5, 4), guide = "none") +
    ggplot2::scale_x_log10(limits = c(0.01, 20), breaks = c(0.01, 0.1, 1, 10)) +
    ggplot2::labs(title = "Comparison: Treatment vs Control", x = "Odds Ratio") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0, size = 12, face = "bold"))
  
  print(p)
  return(p)
}

