context("Diagnostic Test Accuracy")

test_that("Diagnostic accuracy analysis works", {
  # Generate test data
  dta_data <- simulate_dta_data(n_studies = 5)
  
  # Run analysis
  result <- diagnostic_accuracy(
    data = dta_data,
    TP = "TP", FP = "FP", FN = "FN", TN = "TN",
    study_id = "study_id",
    method = "bivariate"
  )
  
  expect_s3_class(result, "cbamm_dta")
  expect_true(result$summary$sensitivity >= 0 && result$summary$sensitivity <= 1)
  expect_true(result$summary$specificity >= 0 && result$summary$specificity <= 1)
})

