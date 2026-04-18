# dose_response.R
# Dose-Response Meta-Analysis Module for CBAMM Package

#' Dose-Response Meta-Analysis
#'
#' Performs dose-response meta-analysis using linear or non-linear models.
#'
#' @param data Data frame containing study data with dose and effect columns
#' @param dose Column name for dose
#' @param effect Column name for effect size (yi)
#' @param se Column name for standard error (se)
#' @param study_id Column name for study identifier
#' @param cases Optional: Column name for number of cases
#' @param n Optional: Column name for total number of subjects
#' @param type_dr Outcome type for dosresmeta: "log-relative", "ir", "cc" (default "log-relative")
#' @param covariance Covariance method: "gl", "h", "indep" (default "indep")
#' @param scale Character: "relative" (exponential scale) or "linear" (natural scale)
#' @param type Model type: "linear" (default) or "spline"
#' @param knots Number of knots for restricted cubic splines (default 3)
#' @return Object of class "cbamm_dr"
#' @export
run_dose_response <- function(data, dose = "dose", effect = "yi", se = "se", 
                             study_id = "study_id", cases = NULL, n = NULL,
                             type_dr = "log-relative",
                             covariance = "indep",
                             scale = c("relative", "linear"),
                             type = c("linear", "spline"), 
                             knots = 3) {
  
  type <- match.arg(type)
  scale <- match.arg(scale)
  
  if (!check_package("dosresmeta", quietly = TRUE)) {
    stop("Package 'dosresmeta' is required for dose-response meta-analysis")
  }
  
  # Ensure data is ordered by study and dose
  data <- data[order(data[[study_id]], data[[dose]]), ]
  
  # Map user columns to internal names to avoid symbol resolution issues
  df_internal <- data.frame(
    .id = data[[study_id]],
    .y = data[[effect]],
    .se = data[[se]],
    .dose = data[[dose]]
  )
  
  if (!is.null(cases)) df_internal$.cases <- data[[cases]]
  if (!is.null(n)) df_internal$.n <- data[[n]]
  
  formula_str <- if (type == "linear") {
    ".y ~ .dose"
  } else {
    paste0(".y ~ dosresmeta::rcs(.dose, ", knots, ")")
  }
  
  # Fit the model
  # We use the internal names and pass the internal data frame
  fit <- try(dosresmeta::dosresmeta(formula = as.formula(formula_str), 
                                   id = .id, 
                                   type = type_dr,
                                   se = .se, 
                                   covariance = covariance,
                                   data = df_internal), silent = TRUE)
  
  if (inherits(fit, "try-error")) {
    stop("Dose-response model fitting failed: ", as.character(fit))
  }
  
  result <- list(
    model = fit,
    data = data,
    type = type,
    scale = scale,
    dose_var = dose,
    effect_var = effect,
    knots = if(type == "spline") knots else NULL
  )
  
  class(result) <- "cbamm_dr"
  return(result)
}

#' Plot Dose-Response Curve
#'
#' @param x cbamm_dr object
#' @param dose_range Range of doses to predict
#' @param ... Additional arguments
#' @export
plot_dose_response <- function(x, dose_range = NULL, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required for plotting")
  }
  
  fit <- x$model
  data <- x$data
  
  if (is.null(dose_range)) {
    dose_range <- seq(min(data[[x$dose_var]]), max(data[[x$dose_var]]), length.out = 100)
  }
  
  # Predict values
  # Map prediction back to internal variable name '.dose'
  new_data <- data.frame(.dose = dose_range)
  
  # Predict based on scale
  is_expo <- x$scale == "relative"
  preds <- try(predict(fit, newdata = new_data, expo = is_expo), silent = TRUE)
  
  if (inherits(preds, "try-error")) {
    stop("Prediction failed for plotting: ", as.character(preds))
  }
  
  plot_df <- data.frame(
    dose = dose_range,
    estimate = preds$pred,
    ci_lb = preds$ci.lb,
    ci_ub = preds$ci.ub
  )
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = dose, y = estimate)) +
    ggplot2::geom_line(color = "blue", linewidth = 1) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_lb, ymax = ci_ub), alpha = 0.2, fill = "blue") +
    ggplot2::geom_hline(yintercept = if(is_expo) 1 else 0, linetype = "dashed", color = "red") +
    ggplot2::labs(
      title = paste("Dose-Response:", ifelse(x$type == "linear", "Linear", "Non-linear Spline")),
      x = x$dose_var,
      y = if(is_expo) "Relative Effect (Exp scale)" else "Mean Difference (Natural scale)"
    ) +
    ggplot2::theme_minimal()
  
  return(p)
}

#' Print method for cbamm_dr
#' @export
print.cbamm_dr <- function(x, ...) {
  cat("Dose-Response Meta-Analysis\n")
  cat("---------------------------\n")
  cat("Model Type:", x$type, "\n")
  if (!is.null(x$knots)) cat("Knots:", x$knots, "\n")
  print(summary(x$model))
}
