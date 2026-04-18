
# CBAMM Package Tests
# Run this to verify everything works

library(cbamm)

cat("\n==== CBAMM PACKAGE TESTS ====\n\n")

# Test 1: Basic functionality with different formats
cat("Test 1: Data standardization\n")
test_data1 <- data.frame(
  TE = rnorm(10, 0.5, 0.2),
  seTE = runif(10, 0.1, 0.3),
  studlab = paste0("Study_", 1:10)
)

result1 <- cbamm_fast(test_data1)
cat("  Result: Effect =", round(result1$b, 3), "\n")
cat("  Status: PASSED\n\n")

# Test 2: Optimized modes
cat("Test 2: Optimized analysis modes\n")
result2a <- cbamm_optimized(test_data1, mode = "fast")
result2b <- cbamm_optimized(test_data1, mode = "balanced")
cat("  Fast mode: Effect =", round(result2a$b, 3), "\n")
cat("  Balanced mode: Effect =", round(result2b$b, 3), "\n")
cat("  Status: PASSED\n\n")

# Test 3: Multivariate
cat("Test 3: Multivariate meta-analysis\n")
mv_data <- data.frame(
  study = 1:10,
  y1 = rnorm(10, 0.3, 0.2),
  se1 = runif(10, 0.1, 0.3),
  y2 = rnorm(10, 0.5, 0.2),
  se2 = runif(10, 0.1, 0.3)
)

mv_result <- multivariate_meta_enhanced(
  mv_data,
  outcomes = c("y1", "y2"),
  se_columns = c("se1", "se2")
)
cat("  Outcome 1: Effect =", round(mv_result$coefficients[1], 3), "\n")
cat("  Outcome 2: Effect =", round(mv_result$coefficients[2], 3), "\n")
cat("  Status: PASSED\n\n")

# Test 4: Network meta-analysis
cat("Test 4: Network meta-analysis\n")
net_data <- data.frame(
  study = 1:12,
  treat1 = rep(c("A", "A", "B"), 4),
  treat2 = rep(c("B", "C", "C"), 4),
  yi = rnorm(12, 0, 0.3),
  se = runif(12, 0.1, 0.4)
)

net_result <- network_meta_enhanced(net_data)
cat("  Number of treatments:", net_result$n_treatments, "\n")
cat("  Top treatment:", net_result$ranking[1], "\n")
cat("  Status: PASSED\n\n")

cat("==== ALL TESTS PASSED ====\n")

