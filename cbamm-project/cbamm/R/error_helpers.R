
#' CBAMM Error Helper
#' 
#' Provides helpful error messages with solutions
#' @export
cbamm_error <- function(error_type, ...) {
  
  errors <- list(
    missing_data = "No data provided.\nUsage: result <- cbamm_fast(your_data)",
    
    wrong_format = "Data is not in the correct format.\nUse: data <- standardize_cbamm_data(your_data)",
    
    missing_columns = "Required columns missing.\nCBAMM needs: yi (effect sizes) and se (standard errors)\nYour columns: {paste(names(data), collapse = ', ')}",
    
    no_studies = "No valid studies after filtering.\nCheck for missing values in yi and se columns.",
    
    convergence = "Model did not converge.\nTry: cbamm_fast(data, method = 'DL') for simpler method",
    
    memory = "Not enough memory for analysis.\nTry: cbamm_optimized(data, mode = 'fast')"
  )
  
  msg <- errors[[error_type]]
  if (!is.null(msg)) {
    stop(msg, call. = FALSE)
  }
}

#' Get CBAMM Help
#' 
#' Quick help for common CBAMM tasks
#' @export
cbamm_help <- function(topic = NULL) {
  
  if (is.null(topic)) {
    cat("CBAMM Quick Help\n")
    cat("================\n\n")
    cat("Common tasks:\n")
    cat("  cbamm_help('start')     - Getting started\n")
    cat("  cbamm_help('data')      - Data preparation\n")
    cat("  cbamm_help('analysis')  - Running analysis\n")
    cat("  cbamm_help('plotting')  - Creating plots\n")
    cat("  cbamm_help('errors')    - Troubleshooting\n")
    return(invisible())
  }
  
  help_topics <- list(
    start = "Getting Started:\n1. Load your data\n2. data <- standardize_cbamm_data(your_data)\n3. result <- cbamm_fast(data)\n4. summary(result)",
    
    data = "Data Preparation:\nYour data needs:\n- Effect sizes (yi, TE, or effect_size)\n- Standard errors (se, sei, or seTE)\n\nUse standardize_cbamm_data() to auto-fix column names",
    
    analysis = "Running Analysis:\n- Quick: cbamm_fast(data)\n- Full: cbamm_optimized(data, mode = 'comprehensive')\n- Bayesian: bayesian_meta_analysis(data)\n- Network: network_meta_analysis(data)",
    
    plotting = "Creating Plots:\n- forest_plot_classic(result)\n- funnel_plot_classic(result)\n- plot_publication_bias_suite(result)",
    
    errors = "Common Errors:\n- Missing columns: Use standardize_cbamm_data()\n- No convergence: Try simpler method\n- Memory issues: Use mode = 'fast'"
  )
  
  msg <- help_topics[[topic]]
  if (!is.null(msg)) {
    cat(msg, "\n")
  } else {
    cat("Topic not found. Use cbamm_help() to see available topics.\n")
  }
}

