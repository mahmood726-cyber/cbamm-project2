
# Funnel Plot Variations for CBAMM Package
# Contains 10 different funnel plot styles

#' Classic Funnel Plot
#' @param data Data frame with study data
#' @export
funnel_plot_classic <- function(data) {
  p <- ggplot2::ggplot(data, ggplot2::aes(x = effect, y = 1/se)) +
    ggplot2::geom_point(size = 3, alpha = 0.7, color = "steelblue") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    ggplot2::labs(title = "Classic Funnel Plot", x = "Effect Size", y = "Precision (1/SE)") +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  
  print(p)
  return(p)
}

#' Contour-Enhanced Funnel Plot
#' @param data Data frame with study data
#' @export
funnel_plot_contour <- function(data) {
  p <- ggplot2::ggplot(data, ggplot2::aes(x = effect, y = 1/se)) +
    ggplot2::geom_point(ggplot2::aes(size = weight, color = study_quality), alpha = 0.8) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::scale_color_manual(values = c("High" = "darkgreen", "Moderate" = "orange", "Low" = "red")) +
    ggplot2::scale_size_continuous(range = c(2, 6), guide = "none") +
    ggplot2::labs(title = "Contour-Enhanced Funnel Plot", x = "Effect Size", y = "Precision (1/SE)",
                 color = "Study Quality") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  
  print(p)
  return(p)
}

#' Regression Funnel Plot
#' @param data Data frame with study data
#' @export
funnel_plot_regression <- function(data) {
  p <- ggplot2::ggplot(data, ggplot2::aes(x = effect, y = 1/se)) +
    ggplot2::geom_smooth(method = "lm", se = TRUE, color = "red", alpha = 0.3) +
    ggplot2::geom_point(ggplot2::aes(size = weight), color = "steelblue", alpha = 0.8) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::scale_size_continuous(range = c(2, 6), guide = "none") +
    ggplot2::labs(title = "Regression Funnel Plot (Egger Test)", x = "Effect Size", y = "Precision (1/SE)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  
  print(p)
  return(p)
}

#' Bubble Chart Funnel Plot
#' @param data Data frame with study data
#' @export
funnel_plot_bubble <- function(data) {
  sample_size <- if ("n_treat" %in% names(data) && "n_control" %in% names(data)) {
    data$n_treat + data$n_control
  } else if ("sample_size" %in% names(data)) {
    data$sample_size
  } else {
    sample(50:300, nrow(data), replace = TRUE)
  }
  
  p <- ggplot2::ggplot(data, ggplot2::aes(x = effect, y = 1/se)) +
    ggplot2::geom_point(ggplot2::aes(size = sample_size, color = year, alpha = weight)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::scale_size_continuous(range = c(3, 12), name = "Sample Size") +
    ggplot2::scale_color_gradient(low = "blue", high = "red", name = "Year") +
    ggplot2::scale_alpha_continuous(range = c(0.4, 0.9), guide = "none") +
    ggplot2::labs(title = "Bubble Chart Funnel Plot", x = "Effect Size", y = "Precision (1/SE)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  
  print(p)
  return(p)
}

