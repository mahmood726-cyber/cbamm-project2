
#' Enhanced Network Meta-Analysis
#' 
#' @param data Data frame with treatment comparisons
#' @param treat1 Column name for treatment 1
#' @param treat2 Column name for treatment 2
#' @param yi Column name for effect size
#' @param se Column name for standard error
#' @param reference Reference treatment
#' @export
network_meta_enhanced <- function(data, 
                                  treat1 = "treat1",
                                  treat2 = "treat2",
                                  yi = "yi",
                                  se = "se",
                                  reference = NULL) {
  
  # Get unique treatments
  treatments <- sort(unique(c(data[[treat1]], data[[treat2]])))
  n_treat <- length(treatments)
  
  if (n_treat < 2) {
    stop("At least 2 treatments required", call. = FALSE)
  }
  
  if (is.null(reference)) {
    reference <- treatments[1]
  }
  
  # Create contrasts matrix
  n_studies <- nrow(data)
  X <- matrix(0, nrow = n_studies, ncol = n_treat - 1)
  
  # Set column names (all treatments except reference)
  other_treats <- setdiff(treatments, reference)
  colnames(X) <- other_treats
  
  # Fill design matrix
  for (i in 1:n_studies) {
    t1 <- data[[treat1]][i]
    t2 <- data[[treat2]][i]
    
    if (t1 != reference && t1 %in% other_treats) {
      X[i, which(other_treats == t1)] <- 1
    }
    if (t2 != reference && t2 %in% other_treats) {
      X[i, which(other_treats == t2)] <- -1
    }
  }
  
  # Weighted least squares
  w <- 1 / data[[se]]^2
  XtW <- t(X * w)
  
  # Check if solvable
  if (qr(XtW %*% X)$rank < ncol(X)) {
    stop("Network is not connected or has insufficient data", call. = FALSE)
  }
  
  # Estimates
  beta <- solve(XtW %*% X) %*% (XtW %*% data[[yi]])
  se_beta <- sqrt(diag(solve(XtW %*% X)))
  
  # Calculate P-scores for ranking
  p_scores <- numeric(n_treat)
  names(p_scores) <- treatments
  
  # Reference gets compared to all others
  ref_effects <- -beta[,1]
  p_scores[reference] <- mean(pnorm(ref_effects / se_beta))
  
  # Other treatments
  for (i in 1:length(other_treats)) {
    treat_effects <- beta[i,1]
    p_scores[other_treats[i]] <- pnorm(treat_effects / se_beta[i])
  }
  
  # Normalize P-scores
  p_scores <- p_scores / max(p_scores)
  
  result <- list(
    treatments = treatments,
    reference = reference,
    coefficients = as.vector(beta),
    se = se_beta,
    p_scores = sort(p_scores, decreasing = TRUE),
    ranking = names(sort(p_scores, decreasing = TRUE)),
    k = n_studies,
    n_treatments = n_treat
  )
  
  names(result$coefficients) <- other_treats
  
  class(result) <- c("cbamm_network", "cbamm")
  return(result)
}

#' Print Network Meta-Analysis
#' @export
print.cbamm_network <- function(x, digits = 3, ...) {
  cat("
Network Meta-Analysis Results
")
  cat("==============================
")
  cat("Studies:", x$k, "
")
  cat("Treatments:", x$n_treatments, "
")
  cat("Reference:", x$reference, "

")
  
  cat("Treatment Rankings (P-scores):
")
  for (i in 1:length(x$ranking)) {
    cat(sprintf("  %d. %s (P=%.3f)
", i, x$ranking[i], x$p_scores[x$ranking[i]]))
  }
  
  cat("
Effects vs Reference (", x$reference, "):
", sep="")
  for (i in 1:length(x$coefficients)) {
    cat(sprintf("  %s: %.3f (SE=%.3f)
",
                names(x$coefficients)[i], x$coefficients[i], x$se[i]))
  }
  
  invisible(x)
}

