# Clean test of advanced methods
# Run after fixing namespace issues

# Clear environment
rm(list = ls())

# Load the package
library(cbamm)

cat("Testing Advanced Methods (Clean)\n")
cat("================================\n\n")

# Test 1: Prediction Intervals
tryCatch({
  cat("1. Prediction Interval: ")
  pi <- calculate_prediction_interval(effect = 0.5, se = 0.1, tau2 = 0.04)
  cat("✓ [", round(pi$prediction_interval, 3), "]\n", sep = " ")
}, error = function(e) {
  cat("✗ Error:", e$message, "\n")
})

# Test 2: IPD Meta-Analysis
tryCatch({
  cat("2. IPD Meta-Analysis: ")
  ipd_data <- simulate_ipd_data(n_studies = 3, n_patients = 50)
  ipd_result <- ipd_meta_analysis(ipd_data)
  cat("✓ Effect =", round(ipd_result$estimate, 3), "\n")
}, error = function(e) {
  cat("✗ Error:", e$message, "\n")
})

# Test 3: Diagnostic Accuracy
tryCatch({
  cat("3. Diagnostic Accuracy: ")
  dta_data <- simulate_dta_data(n_studies = 5)
  dta_result <- diagnostic_accuracy(dta_data, TP = "TP", FP = "FP", FN = "FN", TN = "TN")
  cat("✓ Sens =", round(dta_result$sensitivity, 3), 
      ", Spec =", round(dta_result$specificity, 3), "\n")
}, error = function(e) {
  cat("✗ Error:", e$message, "\n")
})

# Test 4: Living Review
tryCatch({
  cat("4. Living Review: ")
  seq_data <- simulate_sequential_data(n_studies = 5)
  living_result <- living_systematic_review(seq_data)
  cat("✓ Stopped =", living_result$stopped, "\n")
}, error = function(e) {
  cat("✗ Error:", e$message, "\n")
})

# Test 5: Cross-Design
tryCatch({
  cat("5. Cross-Design Synthesis: ")
  rct <- data.frame(effect = rnorm(5, 0.4, 0.05), variance = runif(5, 0.01, 0.03))
  obs <- data.frame(effect = rnorm(5, 0.45, 0.08), variance = runif(5, 0.02, 0.05))
  cross_result <- cross_design_synthesis(rct, obs)
  cat("✓ Combined =", round(cross_result$combined$effect, 3), "\n")
}, error = function(e) {
  cat("✗ Error:", e$message, "\n")
})

cat("\nTest complete!\n")

