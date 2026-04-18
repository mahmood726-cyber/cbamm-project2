library(testthat)
library(cbamm)

context("New CBAMM Features v5.7")

test_that("categorical transportability weighting works", {
  # Create data with categorical covariate
  data <- data.frame(
    study_id = paste0("S", 1:10),
    yi = rnorm(10, 0.5, 0.1),
    se = runif(10, 0.05, 0.2),
    region = factor(sample(c("US", "EU", "Asia"), 10, replace = TRUE)),
    age_mean = rnorm(10, 60, 5),
    female_pct = runif(10, 0.4, 0.6),
    bmi_mean = rnorm(10, 25, 2),
    charlson = runif(10, 1, 3)
  )
  
  # Target population with same categorical levels
  target_pop <- list(
    region = "EU",
    age_mean = 65,
    female_pct = 0.55,
    bmi_mean = 27,
    charlson = 2.0
  )
  
  # Compute weights
  expect_silent(w <- compute_transport_weights(data, target_pop))
  
  expect_equal(length(w), 10)
  expect_true(all(w >= 0))
  expect_equal(sum(w), 1, tolerance = 1e-7)
})

test_that("multiverse analysis works", {
  data <- data.frame(
    study_id = paste0("S", 1:10),
    yi = rnorm(10, 0.5, 0.1),
    se = runif(10, 0.05, 0.2),
    study_type = factor(c(rep("RCT", 6), rep("OBS", 4)))
  )
  
  config <- cbamm_config(methods = list(multiverse = TRUE, parallel = FALSE))
  
  expect_silent(mv <- run_multiverse_analysis(data, config))
  
  expect_s3_class(mv, "cbamm_multiverse")
  expect_true(nrow(mv) > 0)
})

test_that("Dose-Response Meta-Analysis works", {
  skip_if_not_installed("dosresmeta")
  
  # Mock dose-response data
  data <- data.frame(
    study_id = rep(1:3, each = 3),
    dose = rep(c(0, 10, 20), 3),
    yi = c(0, 0.2, 0.5, 0, 0.3, 0.4, 0, 0.1, 0.6),
    se = rep(c(0, 0.1, 0.1), 3)
  )
  
  expect_silent(result <- run_dose_response(data, type = "linear"))
  expect_s3_class(result, "cbamm_dr")
  
  p <- plot_dose_response(result)
  expect_s3_class(p, "ggplot")
})

test_that("Multilevel Meta-Analysis (3-level) works", {
  # Mock nested data (outcomes within studies)
  data <- data.frame(
    study_id = rep(1:5, each = 2),
    effect_id = 1:10,
    yi = rnorm(10, 0.5, 0.1),
    se = 0.1
  )
  
  expect_silent(result <- run_multilevel_meta(data, outer_level = "study_id", inner_level = "effect_id"))
  expect_s3_class(result, "cbamm_multilevel")
  expect_true("Level 3 (study_id)" %in% names(result$variance_distribution))
})

test_that("Bias-Adjusted Cross-Design Synthesis works", {
  rct <- data.frame(effect = 0.5, variance = 0.01)
  obs <- data.frame(effect = 0.7, variance = 0.02)
  
  # With bias offset (expecting result closer to RCT)
  result <- cross_design_synthesis(rct, obs, bias_offset = 0.2)
  
  expect_equal(result$obs$effect, 0.5) # 0.7 - 0.2
  expect_equal(result$combined$effect, 0.5) # IVW of 0.5 and 0.5
})

test_that("Component NMA works", {
  data <- data.frame(
    study_id = rep(1:5, each = 2),
    treatment = c("A", "Placebo", "B", "Placebo", "A + B", "A", "A + B", "Placebo", "B", "A"),
    effect = c(0.5, 0, 0.3, 0, 0.8, 0.4, 0.9, 0, 0.4, 0.6),
    se = 0.1
  )
  
  expect_silent(result <- network_meta_analysis(data, studies="study_id", treatments="treatment",
                                               effect_size="effect", std_error="se",
                                               reference="Placebo", model="component"))
  
  expect_equal(result$model, "component")
})

test_that("One-stage IPD meta-analysis works", {
  skip_if_not_installed("lme4")
  ipd_list <- simulate_ipd_data(n_studies = 3, n_patients = 50)
  result <- ipd_meta_analysis(ipd_list, method = "one-stage", model_type = "random")
  expect_equal(result$method, "one-stage")
})
