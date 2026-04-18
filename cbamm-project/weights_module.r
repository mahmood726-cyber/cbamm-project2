#' Compute Transport Weights for Generalizability
#'
#' @param data Study data with covariate information
#' @param target_population List with target population characteristics
#' @param truncation Quantile level for weight truncation (default 0.02)
#'
#' @return Vector of transport weights
#' @export
compute_transport_weights <- function(data, target_population, truncation = 0.02) {
  req <- c("age_mean", "female_pct", "bmi_mean", "charlson")
  miss <- setdiff(req, names(data))
  
  if (length(miss)) {
    warning("Missing transportability variables: ", paste(miss, collapse = ", "), ". Using uniform weights.")
    return(rep(1 / nrow(data), nrow(data)))
  }
  
  X <- data |>
    dplyr::select(dplyr::all_of(req)) |>
    dplyr::mutate(intercept = 1) |>
    dplyr::select(intercept, dplyr::everything()) |>
    as.data.frame()
  
  target_moments <- c(1,
                      target_population$age_mean,
                      target_population$female_pct,
                      target_population$bmi_mean,
                      target_population$charlson)
  
  # Try WeightIt first (if available)
  if (check_package("WeightIt", quietly = TRUE)) {
    weights <- .compute_weightit_weights(X, target_population, truncation)
    if (!is.null(weights)) return(weights)
  }
  
  # Fallback: robust entropy optimizer
  .compute_entropy_weights(X, target_moments, truncation, data)
}

#' Compute WeightIt-based Transport Weights
#' @keywords internal
.compute_weightit_weights <- function(X, target_population, truncation) {
  with_fallback(
    pkg = "WeightIt",
    fun_call = {
      set.seed(123)
      M <- 1000
      tgt <- data.frame(
        intercept  = 1,
        age_mean   = rnorm(M, target_population$age_mean, 3),
        female_pct = pmin(pmax(rnorm(M, target_population$female_pct, 0.05), 0), 1),
        bmi_mean   = rnorm(M, target_population$bmi_mean, 1.2),
        charlson   = pmax(rnorm(M, target_population$charlson, 0.4), 0)
      )
      
      pool <- dplyr::bind_rows(
        X   |> dplyr::mutate(treat = 0L),
        tgt |> dplyr::mutate(treat = 1L)
      )
      
      fml <- as.formula("treat ~ age_mean + female_pct + bmi_mean + charlson")
      wobj <- try(WeightIt::weightit(fml, data = pool, method = "ebal"), silent = TRUE)
      
      if (!inherits(wobj, "try-error")) {
        w <- as.numeric(wobj$weights[pool$treat == 0])
        w <- w / sum(w)
        
        if (truncation > 0) {
          lo <- stats::quantile(w, truncation)
          hi <- stats::quantile(w, 1 - truncation)
          w <- pmax(pmin(w, hi), lo)
          w <- w / sum(w)
        }
        
        # Print balance if cobalt available
        if (check_package("cobalt", quietly = TRUE)) {
          suppressMessages(try(print(cobalt::bal.tab(wobj)), silent = TRUE))
        }
        
        return(w)
      }
      NULL
    },
    fallback_fun = function() NULL,
    fallback_msg = "WeightIt method failed, using entropy fallback"
  )
}

#' Compute Entropy-based Transport Weights
#' @keywords internal
.compute_entropy_weights <- function(X, target_moments, truncation, data) {
  init <- rep(1 / nrow(X), nrow(X))
  
  obj <- function(lambda) {
    eta <- as.numeric(as.matrix(X) %*% lambda)
    eta <- pmin(eta, 700)  # Prevent overflow
    w <- init * exp(eta)
    w <- as.numeric(w / sum(w))
    sum((colSums(as.matrix(X) * w) - target_moments)^2)
  }
  
  opt <- try(stats::optim(par = rep(0, ncol(X)), fn = obj,
                          control = list(maxit = 2000, reltol = 1e-10)), silent = TRUE)
  
  if (!inherits(opt, "try-error")) {
    w <- init * exp(as.matrix(X) %*% opt$par)
    w <- as.numeric(w / sum(w))
    
    if (truncation > 0) {
      lo <- stats::quantile(w, truncation)
      hi <- stats::quantile(w, 1 - truncation)
      w <- pmax(pmin(w, hi), lo)
      w <- w / sum(w)
    }
    
    return(w)
  }
  
  warning("Transport optimization failed; using uniform weights.")
  rep(1 / nrow(data), nrow(data))
}

#' Compute Analysis Weights from Transport Weights
#'
#' @param data Study data
#' @param transport_weights Transport weights (optional)
#'
#' @return Vector of analysis weights
#' @export
compute_analysis_weights <- function(data, transport_weights = NULL) {
  w <- if (is.null(transport_weights)) rep(1, nrow(data)) else as.numeric(transport_weights)
  w <- ifelse(is.finite(w) & w > 0, w, 0)
  w / mean(w)
}

#' Apply GRADE-based Weighting
#'
#' @param weights Base weights to modify
#' @param grade_scores GRADE quality scores
#'
#' @return Modified weights incorporating GRADE scores
#' @export
apply_grade_weighting <- function(weights, grade_scores) {
  grade_multiplier <- dplyr::case_when(
    grade_scores == "High"     ~ 1.0,
    grade_scores == "Moderate" ~ 0.8,
    grade_scores == "Low"      ~ 0.6,
    grade_scores == "Very low" ~ 0.4,
    TRUE ~ 0.8
  )
  
  z <- as.numeric(weights) * as.numeric(grade_multiplier)
  z / mean(z)
}