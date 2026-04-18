# tests/testthat/test-comprehensive.R
# Comprehensive tests for all CBAMM functionality

library(testthat)
library(cbamm)

# Test 1: Basic Package Loading and Data Simulation
test_that("Package loads and data simulation works", {
  
  # Test data simulation
  expect_silent({
    sim_data <- simulate_cbamm_data(n_rct = 8, n_obs = 6, n_mr = 4)
  })
  
  expect_s3_class(sim_data, "data.frame")
  expect_equal(nrow(sim_data), 18)
  expect_true(all(c("study_id", "yi", "se", "study_type") %in% names(sim_data)))
  expect_silent(validate_cbamm_data(sim_data))
})

# Test 2: Conflict Detection (FIXED)
test_that("Conflict detection works with fixed parameters", {
  
  # Create data with potential conflicts
  conflicted_data <- data.frame(
    study_id = paste0("Study_", 1:10),
    yi = c(rep(0.2, 5), rep(0.8, 5)),  # Two distinct groups
    se = rep(0.1, 10),
    study_type = factor(rep("RCT", 10)),
    analysis_weights = rep(1, 10)
  )
  
  # Test conflict detection function directly
  expect_silent({
    conflicts <- detect_study_conflicts(
      data = conflicted_data,
      threshold = 0.15,
      k_candidates = 2:4
    )
  })
  
  expect_true(!is.null(conflicts))
  if (!is.null(conflicts)) {
    expect_true("delta" %in% names(conflicts))
    expect_true("threshold_met" %in% names(conflicts))
  }
})

# Test 3: Core Meta-Analysis
test_that("Core meta-analysis runs without errors", {
  
  # Simple test data
  test_data <- data.frame(
    study_id = paste0("Study_", 1:8),
    yi = c(0.2, 0.5, -0.1, 0.8, 0.3, 0.1, 0.4, 0.2),
    se = c(0.1, 0.15, 0.12, 0.18, 0.09, 0.11, 0.13, 0.10),
    study_type = factor(rep(c("RCT", "OBS"), each = 4), levels = c("RCT", "OBS"))
  )
  
  # Basic configuration
  config <- cbamm_config(
    methods = list(
      transport = FALSE,
      hksj = TRUE,
      bayesian = FALSE,
      robust_variance = FALSE,
      pet_peese = TRUE,
      conflict_detection = TRUE,
      missing_studies = FALSE
    ),
    output = list(verbose = FALSE, plots = FALSE)
  )
  
  expect_silent({
    results <- cbamm(test_data, config = config)
  })
  
  expect_s3_class(results, "cbamm_results")
  expect_true(!is.null(results$meta_results))
  expect_true(is.numeric(results$diagnostics$computation_time))
})

# Test 4: S3 Methods
test_that("S3 methods work correctly", {
  
  simple_data <- simulate_cbamm_data(n_rct = 5, n_obs = 5, seed = 789)
  
  config <- cbamm_config(
    methods = list(bayesian = FALSE, conflict_detection = FALSE),
    output = list(verbose = FALSE, plots = FALSE)
  )
  
  results <- cbamm(simple_data, config = config)
  
  # Test print method
  expect_output(print(results), "CBAMM Meta-Analysis Results")
  
  # Test summary method  
  expect_output(summary(results), "CBAMM Meta-Analysis Results")
  
  # Test plot method (should handle gracefully when no plots)
  expect_message(plot(results), "No plots available")
})
