
#' Standardize Meta-Analysis Data
#'
#' Standardizes column names and formats for consistent use across CBAMM functions
#'
#' @param data Input data frame
#' @param verbose Logical; print standardization messages
#' @return Standardized data frame with consistent column names
#' @export
#'
#' @details
#' This function standardizes various column naming conventions:
#' - study, study_id, studlab -> study_id
#' - yi, TE, effect, effect_size -> yi  
#' - sei, seTE, se, standard_error -> se
#' - vi, variance, var -> vi
#' - ni, n, sample_size -> ni
#'
#' @examples
#' data <- data.frame(TE = rnorm(10), seTE = runif(10))
#' std_data <- standardize_cbamm_data(data)
#'
standardize_cbamm_data <- function(data, verbose = FALSE) {
  
  if (!is.data.frame(data)) {
    stop("Input must be a data frame", call. = FALSE)
  }
  
  original_names <- names(data)
  
  # Create standardized copy
  std_data <- data
  
  # Standardize study identifier
  study_cols <- c("study", "studlab", "study_label", "trial", "author")
  for (col in study_cols) {
    if (col %in% names(std_data) && !"study_id" %in% names(std_data)) {
      std_data$study_id <- std_data[[col]]
      if (verbose) message("Renamed ", col, " to study_id")
      break
    }
  }
  
  # If no study column exists, create one
  if (!"study_id" %in% names(std_data)) {
    std_data$study_id <- paste0("Study_", seq_len(nrow(std_data)))
    if (verbose) message("Created study_id column")
  }
  
  # Standardize effect size
  effect_cols <- c("TE", "effect", "effect_size", "estimate", "theta", "ES")
  for (col in effect_cols) {
    if (col %in% names(std_data) && !"yi" %in% names(std_data)) {
      std_data$yi <- std_data[[col]]
      if (verbose) message("Renamed ", col, " to yi")
      break
    }
  }
  
  # Standardize standard error
  se_cols <- c("seTE", "sei", "standard_error", "std_error", "SE")
  for (col in se_cols) {
    if (col %in% names(std_data) && !"se" %in% names(std_data)) {
      std_data$se <- std_data[[col]]
      if (verbose) message("Renamed ", col, " to se")
      break
    }
  }
  
  # Create vi from se if missing
  if ("se" %in% names(std_data) && !"vi" %in% names(std_data)) {
    std_data$vi <- std_data$se^2
    if (verbose) message("Created vi from se")
  }
  
  # Standardize variance
  var_cols <- c("variance", "var", "v")
  for (col in var_cols) {
    if (col %in% names(std_data) && !"vi" %in% names(std_data)) {
      std_data$vi <- std_data[[col]]
      if (verbose) message("Renamed ", col, " to vi")
      break
    }
  }
  
  # Create se from vi if missing
  if ("vi" %in% names(std_data) && !"se" %in% names(std_data)) {
    std_data$se <- sqrt(std_data$vi)
    if (verbose) message("Created se from vi")
  }
  
  # Standardize sample size
  n_cols <- c("n", "sample_size", "N", "total")
  for (col in n_cols) {
    if (col %in% names(std_data) && !"ni" %in% names(std_data)) {
      std_data$ni <- std_data[[col]]
      if (verbose) message("Renamed ", col, " to ni")
      break
    }
  }
  
  # Add class
  class(std_data) <- c("cbamm_data", class(std_data))
  
  attr(std_data, "original_names") <- original_names
  attr(std_data, "standardized") <- TRUE
  
  return(std_data)
}

