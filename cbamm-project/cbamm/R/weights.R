#' Compute Transport Weights for Generalizability
#'
#' @param data Study data with covariate information
#' @param target_population List with target population characteristics
#' @param truncation Quantile level for weight truncation (default 0.02)
#'
#' @return Vector of transport weights
#' @export
compute_transport_weights <- function(data, target_population, truncation = 0.02) {
  # Default required numeric variables
  req_num <- c("age_mean", "female_pct", "bmi_mean", "charlson")
  
  # Identify target variables (numeric or categorical)
  target_vars <- names(target_population)
  miss <- setdiff(target_vars, names(data))
  
  if (length(miss)) {
    warning("Missing transportability variables in data: ", paste(miss, collapse = ", "), ". Using uniform weights.")
    return(rep(1 / nrow(data), nrow(data)))
  }
  
  # Prepare model matrix X and target moments
  # This handles both numeric and categorical variables
  prep <- .prepare_transport_matrix(data, target_population)
  X <- prep$X
  target_moments <- prep$target_moments
  
  # Try WeightIt first (if available)
  if (check_package("WeightIt", quietly = TRUE)) {
    weights <- .compute_weightit_weights_v2(data, target_population, truncation)
    if (!is.null(weights)) return(weights)
  }
  
  # Fallback: robust entropy optimizer
  .compute_entropy_weights(X, target_moments, truncation, data)
}

#' Prepare Transport Matrix and Moments
#' @keywords internal
.prepare_transport_matrix <- function(data, target_population) {
  vars <- names(target_population)
  
  # Create a dummy data frame for the target to use with model.matrix
  target_df <- as.data.frame(target_population)
  
  # Combine data and target to ensure same factor levels
  combined <- dplyr::bind_rows(
    data |> dplyr::select(dplyr::all_of(vars)) |> dplyr::mutate(.is_target = 0),
    target_df |> dplyr::mutate(.is_target = 1)
  )
  
  # Convert character to factor
  combined <- combined |> 
    dplyr::mutate(across(where(is.character), as.factor))
  
  # Generate model matrix
  fml <- as.formula(paste("~", paste(vars, collapse = " + ")))
  mm <- model.matrix(fml, data = combined)
  
  # Split back
  X <- mm[combined$.is_target == 0, , drop = FALSE]
  target_moments <- mm[combined$.is_target == 1, , drop = FALSE][1, ]
  
  list(X = X, target_moments = target_moments)
}

#' Updated WeightIt-based Transport Weights
#' @keywords internal
.compute_weightit_weights_v2 <- function(data, target_population, truncation) {
  with_fallback(
    pkg = "WeightIt",
    fun_call = {
      vars <- names(target_population)
      
      # Create target data (simple replication or simulation if numeric)
      M <- 1000
      tgt_list <- lapply(names(target_population), function(v) {
        val <- target_population[[v]]
        if (is.numeric(val)) {
          # Heuristic: simulate around the mean
          if (v == "female_pct") return(pmin(pmax(rnorm(M, val, 0.05), 0), 1))
          return(rnorm(M, val, abs(val)*0.1 + 1))
        } else {
          # Categorical: repeat the value
          return(rep(val, M))
        }
      })
      names(tgt_list) <- vars
      tgt <- as.data.frame(tgt_list)
      
      pool <- dplyr::bind_rows(
        data |> dplyr::select(dplyr::all_of(vars)) |> dplyr::mutate(treat = 0L),
        tgt  |> dplyr::mutate(treat = 1L)
      )
      
      fml <- as.formula(paste("treat ~", paste(vars, collapse = " + ")))
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
        
        return(w)
      }
      NULL
    },
    fallback_fun = function() NULL,
    fallback_msg = "WeightIt method failed"
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
compute_analysis_weights <- function(data, transport_weights = NULL) {
  w <- if (is.null(transport_weights)) rep(1, nrow(data)) else as.numeric(transport_weights)
  w <- ifelse(is.finite(w) & w > 0, w, 0)
  w / mean(w)
}

#' Compute Effective Variance from Standard Errors and Analysis Weights
#' @keywords internal
compute_effective_variance <- function(se, analysis_weights = NULL) {
  vi <- se^2

  if (is.null(analysis_weights)) {
    return(vi)
  }

  w <- as.numeric(analysis_weights)
  valid <- is.finite(w) & w > 0
  vi[valid] <- vi[valid] / w[valid]
  vi[!valid] <- NA_real_
  vi
}

#' Apply GRADE-based Weighting
#'
#' @param weights Base weights to modify
#' @param grade_scores GRADE quality scores
#'
#' @return Modified weights incorporating GRADE scores
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
