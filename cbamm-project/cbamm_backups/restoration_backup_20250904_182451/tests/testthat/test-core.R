# Basic Tests for CBAMM Package

library(testthat)
library(cbamm)

# Test Configuration Functions
test_that("cbamm_config creates valid configuration", {
  config <- cbamm_config()
  
  expect_s3_class(config, "cbamm_config")
  expect_true(is.list(config))
  expect_true(all(c("methods", "estimators", "output") %in% names(config)))
})

test_that("cbamm_config validation works", {
  expect_error(
    cbamm_config(estimators = c("INVALID")),
    "Invalid estimators"
  )
  
  expect_error(
    cbamm_config(bayesian = list(method = "invalid", chains = 4, iter = 2000, warmup = 1000, cores = 1, prior = "weakly_informative")),
    "bayesian.*method must be 'brms' or 'jags'"
  )
})

test_that("config update works correctly", {
  config <- cbamm_config()
  updated <- update_config(config, methods = list(transport = FALSE))
  
  expect_false(updated$methods$transport)
  expect_true(updated$methods$hksj)  # Should remain unchanged
})

# Test Data Validation
test_that("data validation catches missing columns", {
  invalid_data <- data.frame(x = 1:5, y = 1:5)
  
  expect_error(
    validate_cbamm_data(invalid_data),
    "Missing required columns"
  )
})

test_that("data validation catches invalid data types", {
  invalid_data <- data.frame(
    study_id = 1:5,
    yi = letters[1:5],  # Should be numeric
    se = 1:5
  )
  
  expect_error(
    validate_cbamm_data(invalid_data),
    "yi.*must be numeric"
  )
})

test_that("data validation works with valid data", {
  valid_data <- data.frame(
    study_id = paste0("Study_", 1:5),
    yi = c(0.2, 0.5, -0.1, 0.8, 0.3),
    se = c(0.1, 0.15, 0.12, 0.18, 0.09)
  )
  
  expect_silent(validate_cbamm_data(valid_data))
})

# Test Utility Functions
test_that("effect size conversion works", {
  # Test identity conversion
  yi <- c(0.2, 0.5, 0.8)
  expect_equal(convert_effect_size(yi, "Cohen_d", "Cohen_d"), yi)
  
  # Test Cohen's d to Hedges' g
  n <- c(50, 60, 40)
  hedges_g <- convert_effect_size(yi, "Cohen_d", "Hedges_g", n = n)
  expect_true(all(hedges_g < yi))  # Hedges' g should be smaller
})

test_that("variance conversion functions work", {
  se <- c(0.1, 0.2, 0.15)
  vi <- se_to_var(se)
  
  expect_equal(vi, se^2)
  expect_equal(var_to_se(vi), se)
})

test_that("p-value formatting works correctly", {
  p_vals <- c(0.05, 0.001, 0.0001, 0.123)
  formatted <- format_p(p_vals)
  
  expect_equal(formatted[1], "0.050")
  expect_equal(formatted[2], "0.001")  
  expect_equal(formatted[3], "< 0.001")
  expect_equal(formatted[4], "0.123")
})

test_that("confidence interval formatting works", {
  ci_string <- format_ci(0.5, 0.2, 0.8, digits = 2)
  expect_equal(ci_string, "0.50 [0.20, 0.80]")
})

test_that("I-squared calculation works", {
  Q <- 10
  df <- 4
  i2 <- calculate_i_squared(Q, df)
  
  expected <- ((Q - df) / Q) * 100
  expect_equal(i2, expected)
  
  # Test negative I-squared handling
  Q_small <- 2
  i2_small <- calculate_i_squared(Q_small, df)
  expect_equal(i2_small, 0)  # Should not be negative
})

test_that("outlier detection works", {
  x <- c(1, 2, 3, 4, 5, 100)  # 100 is clearly an outlier
  outliers <- detect_outliers_iqr(x)
  
  expect_true(outliers[6])    # 100 should be detected
  expect_false(outliers[1])   # 1 should not be detected
})

test_that("safe division handles edge cases", {
  expect_equal(safe_divide(10, 2), 5)
  expect_true(is.na(safe_divide(10, 0)))
  expect_equal(safe_divide(10, 0, na_value = -999), -999)
})

# Test Data Simulation (if you have this function)
test_that("simulated data has correct structure", {
  skip_if_not(exists("simulate_cbamm_data"), "simulate_cbamm_data function not available")
  
  sim_data <- simulate_cbamm_data(n_rct = 5, n_obs = 3)
  
  expect_true(is.data.frame(sim_data))
  expect_equal(nrow(sim_data), 8)
  expect_true(all(c("study_id", "yi", "se") %in% names(sim_data)))
  
  # Validate the simulated data
  expect_silent(validate_cbamm_data(sim_data))
})

# Integration Tests (basic)
test_that("cbamm runs without errors on simulated data", {
  skip_if_not(exists("simulate_cbamm_data"), "simulate_cbamm_data function not available")
  
  # Simple test data
  test_data <- data.frame(
    study_id = paste0("Study_", 1:5),
    yi = c(0.2, 0.5, -0.1, 0.8, 0.3),
    se = c(0.1, 0.15, 0.12, 0.18, 0.09),
    study_type = rep("RCT", 5)
  )
  
  # Test with minimal configuration
  config <- cbamm_config(
    methods = list(
      transport = FALSE,
      hksj = TRUE,
      bayesian = FALSE,
      robust_variance = TRUE,
      pet_peese = FALSE,
      conflict_detection = FALSE,
      missing_studies = FALSE
    ),
    output = list(verbose = FALSE, plots = FALSE)
  )
  
  # This test will only pass once we implement the core functions
  skip_if_not(exists("cbamm"), "cbamm function not available")
  
  expect_silent(results <- cbamm(test_data, config = config))
  expect_s3_class(results, "cbamm_results")
})