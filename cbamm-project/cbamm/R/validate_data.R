
#' Validate CBAMM Data
#'
#' Validates meta-analysis data and provides helpful error messages
#'
#' @param data Input data frame
#' @param require_study Logical; require study identifiers
#' @param require_n Logical; require sample sizes
#' @param auto_fix Logical; attempt to fix common issues
#' @return Validated data or informative error
#' @export
#'
validate_cbamm_data <- function(data, 
                               require_study = TRUE,
                               require_n = FALSE,
                               auto_fix = TRUE) {
  
  # Check if data frame
  if (!is.data.frame(data)) {
    stop("Input must be a data frame.\n",
         "  You provided: ", class(data)[1],
         call. = FALSE)
  }
  
  # Check for empty data
  if (nrow(data) == 0) {
    stop("Data frame is empty.\n",
         "  Please provide data with at least one study.",
         call. = FALSE)
  }
  
  # Auto-standardize only when the input uses aliases or is missing canonical
  # fields. Canonical CBAMM inputs should validate quietly.
  has_canonical_core <- all(c("yi", "se") %in% names(data))
  needs_study_id <- isTRUE(require_study) && !"study_id" %in% names(data)
  needs_standardization <- auto_fix &&
    !inherits(data, "cbamm_data") &&
    (!has_canonical_core || needs_study_id)

  if (needs_standardization) {
    message("Auto-standardizing data format...")
    data <- standardize_cbamm_data(data, verbose = TRUE)
  }
  
  # Required columns
  required <- c("yi", "se")
  missing <- setdiff(required, names(data))
  
  if (length(missing) > 0) {
    # Try to help user identify the issue
    available <- names(data)
    suggestions <- character()
    
    if ("yi" %in% missing) {
      possible_yi <- c("TE", "effect", "effect_size", "estimate")
      found <- intersect(possible_yi, available)
      if (length(found) > 0) {
        suggestions <- c(suggestions,
          paste0("  Found column '", found[1], "' - rename to 'yi' or use standardize_cbamm_data()"))
      }
    }
    
    if ("se" %in% missing) {
      possible_se <- c("seTE", "sei", "standard_error")
      found <- intersect(possible_se, available)
      if (length(found) > 0) {
        suggestions <- c(suggestions,
          paste0("  Found column '", found[1], "' - rename to 'se' or use standardize_cbamm_data()"))
      }
    }
    
    error_msg <- paste0(
      "Missing required columns: ", paste(missing, collapse = ", "), "\n",
      "Available columns: ", paste(available, collapse = ", ")
    )
    
    if (length(suggestions) > 0) {
      error_msg <- paste0(error_msg, "\n\nSuggestions:\n", 
                         paste(suggestions, collapse = "\n"))
    }
    
    stop(error_msg, call. = FALSE)
  }
  
  # Check for NA values
  na_in_yi <- sum(is.na(data$yi))
  na_in_se <- sum(is.na(data$se))
  
  if (na_in_yi > 0 || na_in_se > 0) {
    warning("Found missing values:\n",
            "  ", na_in_yi, " missing in yi\n",
            "  ", na_in_se, " missing in se\n",
            "Consider removing with: data <- data[complete.cases(data[c('yi', 'se')]), ]",
            call. = FALSE)
  }
  
  # Check for invalid values
  if (any(data$se <= 0, na.rm = TRUE)) {
    n_invalid <- sum(data$se <= 0, na.rm = TRUE)
    stop("Found ", n_invalid, " non-positive standard errors.\n",
         "  Standard errors must be positive.\n",
         "  Check rows: ", paste(which(data$se <= 0), collapse = ", "),
         call. = FALSE)
  }
  
  # Check study identifiers if required
  if (require_study && !"study_id" %in% names(data)) {
    data$study_id <- paste0("Study_", seq_len(nrow(data)))
    message("Created study_id column")
  }
  
  # Check sample sizes if required
  if (require_n && !"ni" %in% names(data)) {
    stop("Sample sizes (ni) required but not found.\n",
         "  Add a column 'ni' with sample sizes.",
         call. = FALSE)
  }
  
  return(invisible(data))
}
