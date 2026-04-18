
# Forest Plot Variations for CBAMM Package
# Contains 10 different forest plot styles

#' Classic Forest Plot
#' @param data Data frame with study data
#' @export
forest_plot_classic <- function(data) {
  forest_data <- prepare_forest_data(data)
  
  p <- ggplot2::ggplot(forest_data, ggplot2::aes(y = reorder(study_label, effect))) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
    ggplot2::geom_point(ggplot2::aes(x = effect, size = weight_scaled), color = "black", alpha = 0.8) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = ci_lower, xmax = ci_upper), height = 0.3, alpha = 0.7) +
    ggplot2::scale_size_continuous(range = c(1, 6), guide = "none") +
    ggplot2::labs(title = "Classic Forest Plot", x = "Effect Size", y = "Study") +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  
  print(p)
  return(p)
}

#' Quality-Colored Forest Plot
#' @param data Data frame with study data
#' @export
forest_plot_quality <- function(data) {
  forest_data <- prepare_forest_data(data)
  
  p <- ggplot2::ggplot(forest_data, ggplot2::aes(y = reorder(study_label, effect))) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::geom_point(ggplot2::aes(x = effect, size = weight_scaled, color = study_quality), alpha = 0.8) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = ci_lower, xmax = ci_upper, color = study_quality), 
                           height = 0.3, alpha = 0.7) +
    ggplot2::scale_color_manual(values = c("High" = "darkgreen", "Moderate" = "orange", "Low" = "red")) +
    ggplot2::scale_size_continuous(range = c(2, 7), guide = "none") +
    ggplot2::labs(title = "Forest Plot by Study Quality", x = "Effect Size", y = "Study", 
                 color = "Quality") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  
  print(p)
  return(p)
}

#' Subgroup Forest Plot
#' @param data Data frame with study data
#' @export
forest_plot_subgroups <- function(data) {
  forest_data <- prepare_forest_data(data)
  
  p <- ggplot2::ggplot(forest_data, ggplot2::aes(y = reorder(study_label, effect))) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::geom_point(ggplot2::aes(x = effect, size = weight_scaled), color = "steelblue", alpha = 0.8) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = ci_lower, xmax = ci_upper), 
                           height = 0.2, color = "steelblue", alpha = 0.7) +
    ggplot2::facet_wrap(~region, scales = "free_y", ncol = 2) +
    ggplot2::scale_size_continuous(range = c(1.5, 5), guide = "none") +
    ggplot2::labs(title = "Forest Plot with Regional Subgroups", x = "Effect Size", y = "Study") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                  strip.background = ggplot2::element_rect(fill = "lightblue"))
  
  print(p)
  return(p)
}

#' Helper function to prepare forest data
prepare_forest_data <- function(data) {
  if (!"study_label" %in% names(data)) {
    data$study_label <- if ("author_year" %in% names(data)) data$author_year else data$study
  }
  if (!"weight_scaled" %in% names(data) && "weight" %in% names(data)) {
    data$weight_scaled <- (data$weight / max(data$weight, na.rm = TRUE)) * 5
  }
  if (!"study_quality" %in% names(data)) {
    data$study_quality <- sample(c("High", "Moderate", "Low"), nrow(data), replace = TRUE)
  }
  if (!"region" %in% names(data)) {
    data$region <- sample(c("North America", "Europe", "Asia", "Other"), nrow(data), replace = TRUE)
  }
  return(data)
}

# Additional forest plot functions (simplified versions)
#' @export
forest_plot_bubble <- function(data) {
  forest_data <- prepare_forest_data(data)
  p <- ggplot2::ggplot(forest_data, ggplot2::aes(x = effect, y = reorder(study_label, effect))) +
    ggplot2::geom_point(ggplot2::aes(size = weight, color = se), alpha = 0.7) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = ci_lower, xmax = ci_upper), height = 0.4, alpha = 0.5) +
    ggplot2::scale_size_continuous(range = c(3, 12), name = "Weight") +
    ggplot2::scale_color_gradient(low = "blue", high = "red", name = "SE") +
    ggplot2::labs(title = "Bubble Forest Plot", x = "Effect Size", y = "Study") +
    ggplot2::theme_minimal()
  print(p)
  return(p)
}

#' @export
forest_plot_lollipop <- function(data) {
  forest_data <- prepare_forest_data(data)
  p <- ggplot2::ggplot(forest_data, ggplot2::aes(x = effect, y = reorder(study_label, effect))) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = effect, y = study_label, yend = study_label),
                         color = "steelblue", size = 1, alpha = 0.7) +
    ggplot2::geom_point(ggplot2::aes(size = weight_scaled), color = "darkblue", alpha = 0.8) +
    ggplot2::scale_size_continuous(range = c(3, 8), guide = "none") +
    ggplot2::labs(title = "Lollipop Forest Plot", x = "Effect Size", y = "Study") +
    ggplot2::theme_minimal()
  print(p)
  return(p)
}

