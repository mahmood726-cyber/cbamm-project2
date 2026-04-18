#!/usr/bin/env Rscript
# ============================================================================
# Script to Fix NAMESPACE and Integration Issues
# ============================================================================

cat("========================================\n")
cat("Fixing NAMESPACE and Integration Issues\n")
cat("========================================\n\n")

# -----------------------------------------------------------------------------
# STEP 1: Clean up the environment
# -----------------------------------------------------------------------------
cat("Step 1: Cleaning up environment...\n")

# Remove conflicting objects if they exist
if (exists("cross_design_synthesis")) {
  rm(cross_design_synthesis)
  cat("  Removed conflicting cross_design_synthesis\n")
}
if (exists("ipd_meta_analysis")) {
  rm(ipd_meta_analysis)
  cat("  Removed conflicting ipd_meta_analysis\n")
}

cat("\n")

# -----------------------------------------------------------------------------
# STEP 2: Fix the NAMESPACE file
# -----------------------------------------------------------------------------
cat("Step 2: Fixing NAMESPACE file...\n")

# Read current NAMESPACE
namespace_lines <- readLines("NAMESPACE")

# Remove problematic entries
namespace_lines <- namespace_lines[!grepl("advanced_meta_analysis", namespace_lines)]
namespace_lines <- namespace_lines[!grepl("print.cbamm_advanced", namespace_lines)]

# Write back cleaned NAMESPACE
writeLines(namespace_lines, "NAMESPACE")
cat("  ✓ Removed non-existent exports from NAMESPACE\n")

cat("\n")

# -----------------------------------------------------------------------------
# STEP 3: Add the missing optimize_config_for_performance function
# -----------------------------------------------------------------------------
cat("Step 3: Adding missing optimize_config_for_performance function...\n")

config_fix <- '
#\' Optimize Configuration for Performance
#\' 
#\' Internal function to optimize configuration settings
#\' @keywords internal
optimize_config_for_performance <- function(config, performance_mode = "standard", data = NULL) {
  # Simple implementation - just return the config as-is
  # This fixes the missing function error
  if (is.null(config)) {
    config <- list(
      method = "fixed",
      conf_level = 0.95,
      continuity_correction = "auto",
      show_forest = TRUE,
      show_funnel = TRUE
    )
  }
  
  # Apply performance mode if specified
  if (performance_mode == "fast") {
    config$show_forest <- FALSE
    config$show_funnel <- FALSE
  } else if (performance_mode == "complete") {
    config$show_forest <- TRUE
    config$show_funnel <- TRUE
  }
  
  return(config)
}
'

# Check if cbamm.R exists and add the function
if (file.exists("R/cbamm.R")) {
  cbamm_content <- readLines("R/cbamm.R")
  
  # Check if function already exists
  if (!any(grepl("optimize_config_for_performance", cbamm_content))) {
    # Add the function at the end
    cbamm_content <- c(cbamm_content, "", config_fix)
    writeLines(cbamm_content, "R/cbamm.R")
    cat("  ✓ Added optimize_config_for_performance to cbamm.R\n")
  } else {
    cat("  ✓ optimize_config_for_performance already exists\n")
  }
} else {
  # Create a new file with the function
  writeLines(config_fix, "R/optimize_config.R")
  cat("  ✓ Created optimize_config.R with missing function\n")
}

cat("\n")

# -----------------------------------------------------------------------------
# STEP 4: Fix the test file
# -----------------------------------------------------------------------------
cat("Step 4: Fixing test file...\n")

# Fix the diagnostic accuracy test
if (file.exists("tests/testthat/test-advanced.R")) {
  test_content <- readLines("tests/testthat/test-advanced.R")
  
  # Find and fix the diagnostic accuracy test
  for (i in seq_along(test_content)) {
    if (grepl("diagnostic_accuracy\\(dta_data\\)", test_content[i])) {
      # Fix the function call to include required arguments
      test_content[i] <- '  result <- diagnostic_accuracy(dta_data, TP = "TP", FP = "FP", FN = "FN", TN = "TN")'
    }
  }
  
  writeLines(test_content, "tests/testthat/test-advanced.R")
  cat("  ✓ Fixed diagnostic_accuracy test\n")
}

cat("\n")

# -----------------------------------------------------------------------------
# STEP 5: Create a clean test script
# -----------------------------------------------------------------------------
cat("Step 5: Creating clean test script...\n")

clean_test <- '# Clean test of advanced methods
# Run after fixing namespace issues

# Clear environment
rm(list = ls())

# Load the package
library(cbamm)

cat("Testing Advanced Methods (Clean)\\n")
cat("================================\\n\\n")

# Test 1: Prediction Intervals
tryCatch({
  cat("1. Prediction Interval: ")
  pi <- calculate_prediction_interval(effect = 0.5, se = 0.1, tau2 = 0.04)
  cat("✓ [", round(pi$prediction_interval, 3), "]\\n", sep = " ")
}, error = function(e) {
  cat("✗ Error:", e$message, "\\n")
})

# Test 2: IPD Meta-Analysis
tryCatch({
  cat("2. IPD Meta-Analysis: ")
  ipd_data <- simulate_ipd_data(n_studies = 3, n_patients = 50)
  ipd_result <- ipd_meta_analysis(ipd_data)
  cat("✓ Effect =", round(ipd_result$estimate, 3), "\\n")
}, error = function(e) {
  cat("✗ Error:", e$message, "\\n")
})

# Test 3: Diagnostic Accuracy
tryCatch({
  cat("3. Diagnostic Accuracy: ")
  dta_data <- simulate_dta_data(n_studies = 5)
  dta_result <- diagnostic_accuracy(dta_data, TP = "TP", FP = "FP", FN = "FN", TN = "TN")
  cat("✓ Sens =", round(dta_result$sensitivity, 3), 
      ", Spec =", round(dta_result$specificity, 3), "\\n")
}, error = function(e) {
  cat("✗ Error:", e$message, "\\n")
})

# Test 4: Living Review
tryCatch({
  cat("4. Living Review: ")
  seq_data <- simulate_sequential_data(n_studies = 5)
  living_result <- living_systematic_review(seq_data)
  cat("✓ Stopped =", living_result$stopped, "\\n")
}, error = function(e) {
  cat("✗ Error:", e$message, "\\n")
})

# Test 5: Cross-Design
tryCatch({
  cat("5. Cross-Design Synthesis: ")
  rct <- data.frame(effect = rnorm(5, 0.4, 0.05), variance = runif(5, 0.01, 0.03))
  obs <- data.frame(effect = rnorm(5, 0.45, 0.08), variance = runif(5, 0.02, 0.05))
  cross_result <- cross_design_synthesis(rct, obs)
  cat("✓ Combined =", round(cross_result$combined$effect, 3), "\\n")
}, error = function(e) {
  cat("✗ Error:", e$message, "\\n")
})

cat("\\nTest complete!\\n")
'

writeLines(clean_test, "inst/examples/clean_test.R")
cat("  ✓ Created inst/examples/clean_test.R\n")

cat("\n")

# -----------------------------------------------------------------------------
# Final instructions
# -----------------------------------------------------------------------------
cat("========================================\n")
cat("✓ Fixes Applied!\n")
cat("========================================\n\n")

cat("Issues fixed:\n")
cat("  ✓ Removed non-existent exports from NAMESPACE\n")
cat("  ✓ Added missing optimize_config_for_performance function\n")
cat("  ✓ Fixed diagnostic_accuracy test\n")
cat("  ✓ Created clean test script\n\n")

cat("Next steps:\n")
cat("1. Restart R session to clear environment:\n")
cat("   Session -> Restart R (or Ctrl+Shift+F10)\n\n")

cat("2. Rebuild the package:\n")
cat("   devtools::document()\n")
cat("   devtools::install()\n\n")

cat("3. Test the functions:\n")
cat("   source('inst/examples/clean_test.R')\n\n")

cat("Script completed!\n")
