#' S3 Methods for CBAMM Objects
#' 
#' Print, summary, and plot methods for CBAMM results and configuration objects

#' Print Method for CBAMM Results
#' @param x CBAMM results object
#' @param ... Additional arguments (unused)
#' @return Invisible x
#' @export
print.cbamm_results <- function(x, ...) {
  cat("CBAMM Meta-Analysis Results\n")
  cat("===========================\n\n")
  
  if (!is.null(x$data_summary)) {
    cat("Studies:", x$data_summary$n_studies, "\n")
  }
  
  if (!is.null(x$meta_results$transport)) {
    pred <- try(metafor::predict(x$meta_results$transport, transf = exp), silent = TRUE)
    if (!inherits(pred, "try-error")) {
      cat(sprintf("Transport HR: %.3f (95%% CI: %.3f-%.3f)\n", 
                 pred$pred, pred$ci.lb, pred$ci.ub))
    }
  }
  
  invisible(x)
}

#' Summary Method for CBAMM Results
#' @param object CBAMM results object  
#' @param ... Additional arguments (unused)
#' @return Summary data frame
#' @export
summary.cbamm_results <- function(object, ...) {
  if (!is.null(object$meta_results$transport)) {
    fit <- object$meta_results$transport
    return(data.frame(
      Method = "Transport-weighted",
      Estimate = as.numeric(coef(fit)),
      SE = as.numeric(fit$se),
      I_squared = fit$I2,
      Studies = fit$k
    ))
  }
  return(data.frame())
}

#' Plot Method for CBAMM Results
#' @param x CBAMM results object
#' @param which Which plot to show
#' @param ... Additional arguments  
#' @return Plot object
#' @export
plot.cbamm_results <- function(x, which = "all", ...) {
  if (length(x$plots) == 0) {
    warning("No plots available")
    return(NULL)
  }
  
  if (which == "all") {
    return(x$plots)
  } else if (which %in% names(x$plots)) {
    return(x$plots[[which]])
  } else {
    warning("Plot not found: ", which)
    return(NULL)
  }
}

