context("Advanced Methods")

test_that("prediction intervals work", {
  pi <- calculate_prediction_interval(effect = 0.5, se = 0.1, tau2 = 0.04)
  expect_true(is.list(pi))
  expect_true(pi$prediction_interval[1] < pi$prediction_interval[2])
})

test_that("IPD simulation works", {
  ipd_data <- simulate_ipd_data(n_studies = 2, n_patients = 20)
  expect_equal(length(ipd_data), 2)
  expect_true(is.data.frame(ipd_data[[1]]))
})

test_that("diagnostic accuracy works", {
  dta_data <- simulate_dta_data(n_studies = 3)
  result <- diagnostic_accuracy(dta_data, TP = "TP", FP = "FP", FN = "FN", TN = "TN")
  expect_true(result$sensitivity >= 0 && result$sensitivity <= 1)
  expect_true(result$specificity >= 0 && result$specificity <= 1)
})

