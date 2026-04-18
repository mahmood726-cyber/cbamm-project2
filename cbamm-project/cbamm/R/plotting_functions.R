
# Fixed Plotting Functions for CBAMM Package
# These replace the broken functions with working versions

#' Forest Plot - Fixed Version
#' @param data Data frame with effect_size, se, study columns
#' @param title Plot title
#' @param show_pooled Show pooled estimate
#' @export
forest_plot_classic <- function(data, title = "Forest Plot", show_pooled = TRUE) {
  
  # Input validation
  if (!is.data.frame(data)) {
    stop("Input must be a data frame")
  }
  
  required_cols <- c("effect_size", "se", "study")
  missing_cols <- required_cols[!required_cols %in% names(data)]
  
  if (length(missing_cols) > 0) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Data cleaning
  data <- as.data.frame(data)
  data$effect_size <- as.numeric(data$effect_size)
  data$se <- as.numeric(data$se)
  data$study <- as.character(data$study)
  
  # Remove incomplete cases
  complete_rows <- complete.cases(data[, required_cols])
  data <- data[complete_rows, ]
  
  if (nrow(data) == 0) {
    stop("No complete cases found")
  }
  
  # Calculate derived values
  data$lower_ci <- data$effect_size - 1.96 * data$se
  data$upper_ci <- data$effect_size + 1.96 * data$se
  data$weight <- 1 / (data$se^2)
  data$size <- sqrt(data$weight)
  data$size_scaled <- scales::rescale(data$size, to = c(2, 6))
  
  # Calculate pooled estimate
  pooled_effect <- sum(data$effect_size * data$weight) / sum(data$weight)
  
  # Create plot
  p <- ggplot2::ggplot(data, ggplot2::aes(x = effect_size, y = reorder(study, effect_size))) +
    ggplot2::geom_point(ggplot2::aes(size = size_scaled), color = "darkblue", alpha = 0.7) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lower_ci, xmax = upper_ci), 
                   height = 0.2, color = "black") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
    ggplot2::scale_size_continuous(range = c(2, 6), guide = "none") +
    ggplot2::labs(x = "Effect Size", y = "Study", title = title) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 10),
      plot.title = ggplot2::element_text(size = 12, hjust = 0.5)
    )
  
  if (show_pooled) {
    p <- p + ggplot2::geom_vline(xintercept = pooled_effect, linetype = "solid", 
                        color = "blue", alpha = 0.7)
  }
  
  print(p)
  return(p)
}

#' Funnel Plot - Fixed Version
#' @param data Data frame with effect_size and se columns
#' @param title Plot title
#' @param show_funnel Show funnel boundaries
#' @export
funnel_plot_classic <- function(data, title = "Funnel Plot", show_funnel = TRUE) {
  
  # Input validation
  if (!is.data.frame(data)) {
    stop("Input must be a data frame")
  }
  
  required_cols <- c("effect_size", "se")
  missing_cols <- required_cols[!required_cols %in% names(data)]
  
  if (length(missing_cols) > 0) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Data cleaning
  data <- as.data.frame(data)
  data$effect_size <- as.numeric(data$effect_size)
  data$se <- as.numeric(data$se)
  
  # Remove incomplete cases
  complete_rows <- complete.cases(data[, required_cols])
  data <- data[complete_rows, ]
  
  if (nrow(data) == 0) {
    stop("No complete cases found")
  }
  
  # Calculate precision and pooled estimate
  data$precision <- 1 / data$se
  data$weight <- 1 / (data$se^2)
  pooled_effect <- sum(data$effect_size * data$weight) / sum(data$weight)
  
  # Start with basic plot
  p <- ggplot2::ggplot(data, ggplot2::aes(x = effect_size, y = precision)) +
    ggplot2::geom_point(size = 3, color = "blue", alpha = 0.7) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    ggplot2::geom_vline(xintercept = pooled_effect, linetype = "solid", color = "blue") +
    ggplot2::labs(x = "Effect Size", y = "Precision (1/SE)", title = title) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 10),
      plot.title = ggplot2::element_text(size = 12, hjust = 0.5)
    )
  
  # Add funnel boundaries if requested
  if (show_funnel) {
    se_max <- max(data$se) * 1.2
    se_range <- seq(0.01, se_max, length.out = 100)
    precision_range <- 1 / se_range
    
    funnel_data <- data.frame(
      x = c(pooled_effect - 1.96 * se_range, 
            pooled_effect + 1.96 * rev(se_range)),
      y = c(precision_range, rev(precision_range))
    )
    
    p <- p + ggplot2::geom_polygon(data = funnel_data, ggplot2::aes(x = x, y = y), 
                          fill = "lightgray", alpha = 0.3, color = "gray")
  }
  
  print(p)
  return(p)
}

#' Contour Funnel Plot - Fixed Version
#' @param data Data frame with effect_size and se columns
#' @param title Plot title
#' @export
funnel_plot_contour <- function(data, title = "Contour-Enhanced Funnel Plot") {
  
  # Input validation
  if (!is.data.frame(data)) {
    stop("Input must be a data frame")
  }
  
  required_cols <- c("effect_size", "se")
  missing_cols <- required_cols[!required_cols %in% names(data)]
  
  if (length(missing_cols) > 0) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Clean data
  data$effect_size <- as.numeric(data$effect_size)
  data$se <- as.numeric(data$se)
  data <- data[complete.cases(data[, required_cols]), ]
  
  if (nrow(data) == 0) {
    stop("No complete cases found")
  }
  
  data$precision <- 1 / data$se
  data$weight <- 1 / data$se^2
  pooled_effect <- sum(data$effect_size * data$weight) / sum(data$weight)
  
  # Create basic contour plot
  x_range <- range(data$effect_size)
  x_extend <- diff(x_range) * 0.5
  x_seq <- seq(x_range[1] - x_extend, x_range[2] + x_extend, length.out = 50)
  y_seq <- seq(0.5, max(data$precision) * 1.1, length.out = 50)
  
  grid <- expand.grid(effect = x_seq, precision = y_seq)
  grid$se <- 1 / grid$precision
  grid$z_score <- abs(grid$effect) / grid$se
  grid$p_value <- 2 * (1 - pnorm(grid$z_score))
  
  p <- ggplot2::ggplot() +
    ggplot2::stat_contour_filled(data = grid, 
                        ggplot2::aes(x = effect, y = precision, z = p_value),
                        breaks = c(0, 0.01, 0.05, 0.1, 1),
                        alpha = 0.4) +
    ggplot2::scale_fill_manual(values = c("red", "orange", "yellow", "lightblue"),
                      name = "p-value",
                      labels = c("< 0.01", "0.01-0.05", "0.05-0.10", "> 0.10")) +
    ggplot2::geom_point(data = data, ggplot2::aes(x = effect_size, y = precision), 
               size = 3, color = "black", alpha = 0.8) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    ggplot2::geom_vline(xintercept = pooled_effect, linetype = "solid", color = "red", size = 1) +
    ggplot2::labs(x = "Effect Size", y = "Precision (1/SE)", title = title) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "right",
      axis.text = ggplot2::element_text(size = 10),
      plot.title = ggplot2::element_text(size = 12, hjust = 0.5)
    )
  
  print(p)
  return(p)
}

#' Simple Trim and Fill Plot - Fixed Version
#' @param data Data frame with effect_size and se columns
#' @param title Plot title
#' @export
funnel_plot_trim_fill <- function(data, title = "Trim-and-Fill Analysis") {
  
  # Input validation
  if (!is.data.frame(data)) {
    stop("Input must be a data frame")
  }
  
  required_cols <- c("effect_size", "se")
  missing_cols <- required_cols[!required_cols %in% names(data)]
  
  if (length(missing_cols) > 0) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Clean data
  data$effect_size <- as.numeric(data$effect_size)
  data$se <- as.numeric(data$se)
  data <- data[complete.cases(data[, required_cols]), ]
  
  if (nrow(data) == 0) {
    stop("No complete cases found")
  }
  
  if (nrow(data) < 3) {
    cat("Warning: Too few studies for trim-and-fill analysis\n")
    return(funnel_plot_classic(data, paste(title, "(insufficient data)")))
  }
  
  # Calculate precision
  data$precision <- 1 / data$se
  data$weight <- 1 / data$se^2
  
  # Simple asymmetry check
  pooled_effect <- sum(data$effect_size * data$weight) / sum(data$weight)
  
  # Create plot data - start with observed studies
  plot_data <- data.frame(
    effect_size = data$effect_size,
    precision = data$precision,
    type = "Observed"
  )
  
  # Simple imputation logic
  n_left <- sum(data$effect_size < pooled_effect)
  n_right <- sum(data$effect_size > pooled_effect)
  
  if (abs(n_left - n_right) > 1 && nrow(data) >= 4) {
    if (n_left < n_right) {
      # Impute on left
      right_studies <- data[data$effect_size > pooled_effect, ]
      if (nrow(right_studies) > 0) {
        n_impute <- min(2, nrow(right_studies))
        impute_idx <- sample(nrow(right_studies), n_impute)
        imputed_effects <- 2 * pooled_effect - right_studies$effect_size[impute_idx]
        imputed_precision <- right_studies$precision[impute_idx]
        
        imputed_data <- data.frame(
          effect_size = imputed_effects,
          precision = imputed_precision,
          type = "Imputed"
        )
        plot_data <- rbind(plot_data, imputed_data)
      }
    }
  }
  
  # Calculate adjusted effect
  all_weights <- 1 / (1/plot_data$precision)^2
  adj_effect <- sum(plot_data$effect_size * all_weights) / sum(all_weights)
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = effect_size, y = precision)) +
    ggplot2::geom_point(ggplot2::aes(color = type, shape = type), size = 3, alpha = 0.8) +
    ggplot2::scale_color_manual(values = c("Observed" = "blue", "Imputed" = "red")) +
    ggplot2::scale_shape_manual(values = c("Observed" = 16, "Imputed" = 17)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    ggplot2::geom_vline(xintercept = adj_effect, color = "red", size = 1) +
    ggplot2::labs(x = "Effect Size", y = "Precision (1/SE)", 
         title = title, color = "Study Type", shape = "Study Type") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "bottom",
      axis.text = ggplot2::element_text(size = 10),
      plot.title = ggplot2::element_text(size = 12, hjust = 0.5)
    )
  
  print(p)
  cat("Number of imputed studies:", sum(plot_data$type == "Imputed"), "\n")
  cat("Adjusted pooled effect:", round(adj_effect, 3), "\n")
  return(p)
}

#' P-Curve Plot - Working Version
#' @param data Data frame with effect_size and se columns
#' @param title Plot title
#' @export
p_curve_plot <- function(data, title = "P-Curve Analysis") {
  
  # Input validation
  if (!is.data.frame(data)) {
    stop("Input must be a data frame")
  }
  
  required_cols <- c("effect_size", "se")
  missing_cols <- required_cols[!required_cols %in% names(data)]
  
  if (length(missing_cols) > 0) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Calculate p-values
  data$effect_size <- as.numeric(data$effect_size)
  data$se <- as.numeric(data$se)
  data <- data[complete.cases(data[, required_cols]), ]
  
  if (nrow(data) == 0) {
    stop("No complete cases found")
  }
  
  data$z_score <- data$effect_size / data$se
  data$p_value <- 2 * (1 - pnorm(abs(data$z_score)))
  
  # Filter significant results
  sig_data <- data[data$p_value < 0.05 & !is.na(data$p_value), ]
  
  if (nrow(sig_data) == 0) {
    cat("No significant results found for p-curve analysis\n")
    return(NULL)
  }
  
  # Create histogram
  breaks <- seq(0, 0.05, by = 0.005)
  hist_data <- hist(sig_data$p_value, breaks = breaks, plot = FALSE)
  
  plot_df <- data.frame(
    p_mid = hist_data$mids,
    count = hist_data$counts
  )
  
  # Expected under null
  n_studies <- nrow(sig_data)
  uniform_expected <- n_studies / length(breaks[-1])
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = p_mid)) +
    ggplot2::geom_bar(ggplot2::aes(y = count), stat = "identity", 
             fill = "lightblue", color = "black", alpha = 0.7, width = 0.004) +
    ggplot2::geom_hline(yintercept = uniform_expected, 
               color = "red", linetype = "dashed", size = 1) +
    ggplot2::scale_x_continuous(limits = c(0, 0.05), 
                      breaks = seq(0, 0.05, 0.01)) +
    ggplot2::labs(x = "p-value", y = "Frequency", 
         title = title,
         subtitle = paste("Based on", nrow(sig_data), "significant results")) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 10),
      plot.title = ggplot2::element_text(size = 12, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 10, hjust = 0.5)
    )
  
  print(p)
  return(p)
}

#' Test all fixed plotting functions
#' @export
test_fixed_plotting <- function() {
  
  cat("Testing Fixed CBAMM Plotting Functions\n")
  cat("====================================\n\n")
  
  # Create test data
  set.seed(123)
  test_data <- data.frame(
    study = paste("Study", 1:6),
    effect_size = c(0.2, 0.4, 0.1, 0.6, 0.3, 0.5),
    se = c(0.15, 0.12, 0.18, 0.14, 0.16, 0.13)
  )
  
  functions_to_test <- list(
    list(name = "Forest Plot Classic", func = forest_plot_classic),
    list(name = "Funnel Plot Classic", func = funnel_plot_classic),
    list(name = "Funnel Plot Contour", func = funnel_plot_contour),
    list(name = "Funnel Plot Trim Fill", func = funnel_plot_trim_fill),
    list(name = "P-Curve Plot", func = p_curve_plot)
  )
  
  success_count <- 0
  for (test_func in functions_to_test) {
    cat("Testing", test_func$name, "... ")
    tryCatch({
      result <- test_func$func(test_data)
      cat("SUCCESS\n")
      success_count <- success_count + 1
    }, error = function(e) {
      cat("ERROR:", e$message, "\n")
    })
  }
  
  cat("\n")
  cat("Results:", success_count, "out of", length(functions_to_test), "functions working\n")
  return(success_count == length(functions_to_test))
}

