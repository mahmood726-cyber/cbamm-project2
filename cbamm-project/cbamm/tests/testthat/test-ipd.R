context("IPD Meta-Analysis")

test_that("IPD meta-analysis works with simulated data", {
  # Generate test data
  ipd_data <- simulate_ipd_data(n_studies = 3, n_patients_per_study = 50)
  
  # Run analysis
  result <- ipd_meta_analysis(
    ipd_data = ipd_data,
    outcome = "outcome",
    treatment = "treatment",
    method = "two-stage",
    model_type = "fixed"
  )
  
  expect_s3_class(result, "cbamm_ipd")
  expect_true(!is.null(result$results))
  expect_true(!is.na(result$results$estimate))
})

