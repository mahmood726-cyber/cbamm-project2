# Test advanced methods
library(cbamm)

cat("Testing Advanced Methods\n")
cat("========================\n\n")

# Test 1: Prediction Intervals
cat("1. Prediction Interval:\n")
pi <- calculate_prediction_interval(effect = 0.5, se = 0.1, tau2 = 0.04)
cat("  Interval:", round(pi$prediction_interval, 3), "\n\n")

# Test 2: IPD Meta-Analysis
cat("2. IPD Meta-Analysis:\n")
ipd_data <- simulate_ipd_data(n_studies = 3, n_patients = 50)
ipd_result <- ipd_meta_analysis(ipd_data)
print(ipd_result)
cat("\n")

# Test 3: Diagnostic Accuracy
cat("3. Diagnostic Accuracy:\n")
dta_data <- simulate_dta_data(n_studies = 5)
dta_result <- diagnostic_accuracy(dta_data)
print(dta_result)
cat("\n")

# Test 4: Living Review
cat("4. Living Review:\n")
seq_data <- simulate_sequential_data(n_studies = 5)
living_result <- living_systematic_review(seq_data)
print(living_result)
cat("\n")

# Test 5: Cross-Design
cat("5. Cross-Design Synthesis:\n")
rct <- data.frame(effect = rnorm(5, 0.4, 0.05), variance = runif(5, 0.01, 0.03))
obs <- data.frame(effect = rnorm(5, 0.45, 0.08), variance = runif(5, 0.02, 0.05))
cross_result <- cross_design_synthesis(rct, obs)
cat("  RCT effect:", round(cross_result$rct$effect, 3), "\n")
cat("  Obs effect:", round(cross_result$obs$effect, 3), "\n")
cat("  Combined:", round(cross_result$combined$effect, 3), "\n")
cat("  Bias:", round(cross_result$bias, 3), "\n\n")

cat("All tests completed!\n")

