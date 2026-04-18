
# CBAMM Package Test Script
# Comprehensive testing of cumulative meta-analysis functionality

library(CBAMM)

# Load example data
data(cumulative_example_data)

cat("Testing CBAMM Cumulative Meta-Analysis\n")
cat("=====================================\n\n")

# Test 1: Basic cumulative analysis by year
cat("Test 1: Basic cumulative analysis by publication year\n")
result1 <- cumulative_meta_analysis(cumulative_example_data, order_by = "year")
print(result1)
cat("\n")

# Test 2: Cumulative analysis by sample size (largest first)
cat("Test 2: Cumulative analysis by sample size (descending)\n")
result2 <- cumulative_meta_analysis(cumulative_example_data, 
                                   order_by = "sample_size", 
                                   order_direction = "descending")
print(result2)
cat("\n")

# Test 3: Cumulative analysis by precision
cat("Test 3: Cumulative analysis by precision\n")
result3 <- cumulative_meta_analysis(cumulative_example_data, order_by = "precision")
print(result3)
cat("\n")

# Test 4: Fixed effects method
cat("Test 4: Fixed effects cumulative analysis\n")
result4 <- cumulative_meta_analysis(cumulative_example_data, 
                                   order_by = "year", 
                                   method = "fixed")
print(result4)
cat("\n")

# Test 5: Bayesian method
cat("Test 5: Bayesian cumulative analysis\n")
result5 <- cumulative_meta_analysis(cumulative_example_data, 
                                   order_by = "year", 
                                   method = "bayesian")
print(result5)
cat("\n")

# Test 6: Plotting functionality
if (requireNamespace("ggplot2", quietly = TRUE)) {
  cat("Test 6: Plotting functionality\n")
  
  # Plot cumulative estimates
  plot(result1, type = "estimate")
  
  # Plot confidence interval width
  plot(result1, type = "ci")
  
  # Plot precision evolution
  plot(result1, type = "precision")
  
  # Plot heterogeneity evolution
  plot(result1, type = "heterogeneity")
  
  cat("✓ All plotting tests completed\n")
} else {
  cat("Note: ggplot2 not available - skipping plot tests\n")
}

cat("\nAll tests completed successfully!\n")

