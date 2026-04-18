library(testthat)
library(cbamm)

test_that("cbamm large-dataset path initializes results before thinning warning", {
  local_mocked_bindings(
    fast_robust_rma = function(...) list(mock_fit = TRUE),
    run_optimized_sensitivity_analyses = function(...) list(),
    generate_advisor_recommendations = function(...) list(),
    .package = "cbamm"
  )

  big_data <- data.frame(
    study_id = paste0("S", seq_len(5001)),
    yi = rnorm(5001),
    se = rep(0.1, 5001),
    study_type = factor(rep("RCT", 5001), levels = c("RCT", "OBS"))
  )

  config <- cbamm_config(
    methods = list(
      transport = FALSE,
      bayesian = FALSE,
      pet_peese = FALSE,
      cross_design = FALSE,
      conflict_detection = FALSE,
      missing_studies = FALSE,
      multiverse = FALSE
    ),
    output = list(verbose = FALSE, plots = FALSE),
    sensitivity = list(publication_bias_methods = character(0))
  )

  expect_silent(result <- suppressMessages(cbamm(big_data, config = config)))
  expect_s3_class(result, "cbamm_results")
  expect_equal(result$data_summary$n_studies, 5000)
  expect_true(any(grepl("Dataset thinned to 5,000 studies", result$diagnostics$warnings, fixed = TRUE)))
})

test_that("cbamm returns transport and cross-design results for mixed designs", {
  data <- data.frame(
    study_id = paste0("Study_", 1:8),
    yi = c(0.20, 0.35, 0.10, 0.30, 0.22, 0.28, 0.18, 0.26),
    se = c(0.10, 0.12, 0.11, 0.10, 0.12, 0.11, 0.10, 0.12),
    study_type = factor(rep(c("RCT", "OBS"), each = 4), levels = c("RCT", "OBS")),
    age_mean = c(58, 60, 62, 59, 63, 64, 61, 65),
    female_pct = c(0.45, 0.50, 0.48, 0.46, 0.55, 0.58, 0.57, 0.60),
    bmi_mean = c(26, 27, 27, 26, 29, 30, 29, 31),
    charlson = c(1.4, 1.6, 1.5, 1.4, 2.0, 2.2, 2.1, 2.3)
  )

  target_population <- list(
    age_mean = 65,
    female_pct = 0.55,
    bmi_mean = 28,
    charlson = 1.8
  )

  config <- cbamm_config(
    methods = list(
      transport = TRUE,
      bayesian = FALSE,
      pet_peese = FALSE,
      cross_design = TRUE,
      conflict_detection = FALSE,
      missing_studies = FALSE,
      multiverse = FALSE
    ),
    output = list(verbose = FALSE, plots = FALSE),
    sensitivity = list(publication_bias_methods = character(0))
  )

  expect_silent(result <- suppressMessages(cbamm(data, config = config, target_population = target_population)))
  expect_s3_class(result, "cbamm_results")
  expect_false(is.null(result$transport_results))
  expect_false(is.null(result$transport_results$transport_weights))
  expect_equal(length(result$transport_results$transport_weights), nrow(data))
  expect_false(is.null(result$cross_design_results))
  expect_true(all(c("rct", "obs", "combined", "heterogeneity", "publication_bias") %in% names(result$cross_design_results)))
})

test_that("cbamm_fast returns cross-design results for mixed designs", {
  data <- data.frame(
    study_id = paste0("Study_", 1:8),
    yi = c(0.20, 0.35, 0.10, 0.30, 0.22, 0.28, 0.18, 0.26),
    se = c(0.10, 0.12, 0.11, 0.10, 0.12, 0.11, 0.10, 0.12),
    study_type = factor(rep(c("RCT", "OBS"), each = 4), levels = c("RCT", "OBS")),
    age_mean = c(58, 60, 62, 59, 63, 64, 61, 65),
    female_pct = c(0.45, 0.50, 0.48, 0.46, 0.55, 0.58, 0.57, 0.60),
    bmi_mean = c(26, 27, 27, 26, 29, 30, 29, 31),
    charlson = c(1.4, 1.6, 1.5, 1.4, 2.0, 2.2, 2.1, 2.3)
  )

  target_population <- list(
    age_mean = 65,
    female_pct = 0.55,
    bmi_mean = 28,
    charlson = 1.8
  )

  config <- cbamm_config(
    methods = list(
      transport = TRUE,
      bayesian = FALSE,
      pet_peese = FALSE,
      cross_design = TRUE,
      conflict_detection = FALSE,
      missing_studies = FALSE,
      multiverse = FALSE
    ),
    output = list(verbose = FALSE, plots = FALSE),
    sensitivity = list(publication_bias_methods = character(0))
  )

  expect_silent(result <- suppressMessages(cbamm_fast(data, config = config, target_population = target_population)))
  expect_s3_class(result, "cbamm_results")
  expect_false(is.null(result$cross_design_results))
  expect_true(all(c("rct", "obs", "combined") %in% names(result$cross_design_results)))
})

test_that("cross_design_synthesis only labels bias-adjusted effect when applied", {
  rct <- data.frame(
    effect = c(0.45, 0.50, 0.55),
    variance = rep(0.02, 3)
  )
  obs <- data.frame(
    effect = c(0.60, 0.62, 0.58),
    variance = rep(0.03, 3)
  )

  result <- cross_design_synthesis(rct, obs, bias_offset = 0.1, use_pet_peese = TRUE)

  expect_false(result$obs$publication_bias$asymmetry_detected)
  expect_true(is.na(result$obs$bias_adjusted_effect))
  expect_equal(result$obs$effect, result$obs$raw_effect)
})

test_that("fast performance mode actually simplifies configuration", {
  config <- cbamm_config()
  optimized <- optimize_config_for_performance(
    config,
    performance_mode = "fast",
    data = data.frame(study_id = 1:12)
  )

  expect_equal(optimized$performance_mode, "fast")
  expect_false(optimized$methods$bayesian)
  expect_false(optimized$output$plots)
  expect_false(optimized$methods$pet_peese)
  expect_false(optimized$methods$missing_studies)
  expect_false(optimized$methods$multiverse)
  expect_equal(optimized$sensitivity$publication_bias_methods, character(0))
})

test_that("summary includes cross-design results when present", {
  data <- data.frame(
    study_id = paste0("Study_", 1:8),
    yi = c(0.20, 0.35, 0.10, 0.30, 0.22, 0.28, 0.18, 0.26),
    se = c(0.10, 0.12, 0.11, 0.10, 0.12, 0.11, 0.10, 0.12),
    study_type = factor(rep(c("RCT", "OBS"), each = 4), levels = c("RCT", "OBS")),
    age_mean = c(58, 60, 62, 59, 63, 64, 61, 65),
    female_pct = c(0.45, 0.50, 0.48, 0.46, 0.55, 0.58, 0.57, 0.60),
    bmi_mean = c(26, 27, 27, 26, 29, 30, 29, 31),
    charlson = c(1.4, 1.6, 1.5, 1.4, 2.0, 2.2, 2.1, 2.3)
  )

  target_population <- list(
    age_mean = 65,
    female_pct = 0.55,
    bmi_mean = 28,
    charlson = 1.8
  )

  config <- cbamm_config(
    methods = list(
      transport = TRUE,
      bayesian = FALSE,
      pet_peese = FALSE,
      cross_design = TRUE,
      conflict_detection = FALSE,
      missing_studies = FALSE,
      multiverse = FALSE
    ),
    output = list(verbose = FALSE, plots = FALSE),
    sensitivity = list(publication_bias_methods = character(0))
  )

  result <- suppressMessages(cbamm(data, config = config, target_population = target_population))
  summary_df <- summary(result)

  expect_true("Cross-design combined" %in% summary_df$Method)
  expect_output(print(result), "Cross-design effect")
})


test_that("cbamm_fast standardizes input aliases and keeps fast mode metadata", {
  data <- data.frame(
    TE = c(0.20, 0.35, 0.10, 0.30),
    seTE = c(0.10, 0.12, 0.11, 0.10)
  )

  config <- cbamm_config(
    methods = list(
      transport = FALSE,
      bayesian = FALSE,
      pet_peese = FALSE,
      cross_design = FALSE,
      conflict_detection = FALSE,
      missing_studies = FALSE,
      multiverse = FALSE
    ),
    output = list(verbose = FALSE, plots = FALSE),
    sensitivity = list(publication_bias_methods = character(0))
  )

  expect_silent(result <- cbamm_fast(data, config = config, method = "DL"))
  expect_s3_class(result, "cbamm_fast")
  expect_s3_class(result, "cbamm_results")
  expect_equal(result$diagnostics$performance_mode, "fast")
  expect_equal(result$config$estimators[[1]], "DL")
  expect_equal(result$data_summary$n_studies, 4)
})

test_that("cbamm_optimized maps legacy aliases onto canonical performance modes", {
  data <- data.frame(
    effect_size = c(0.20, 0.35, 0.10, 0.30),
    standard_error = c(0.10, 0.12, 0.11, 0.10),
    trial = paste0("Study_", 1:4)
  )

  config <- cbamm_config(
    methods = list(
      transport = FALSE,
      bayesian = FALSE,
      pet_peese = FALSE,
      cross_design = FALSE,
      conflict_detection = FALSE,
      missing_studies = FALSE,
      multiverse = FALSE
    ),
    output = list(verbose = FALSE, plots = FALSE),
    sensitivity = list(publication_bias_methods = character(0))
  )

  expect_silent(result <- cbamm_optimized(
    data,
    mode = "full",
    config = config,
    effect_col = "effect_size",
    se_col = "standard_error",
    study_col = "trial"
  ))
  expect_s3_class(result, "cbamm_optimized")
  expect_s3_class(result, "cbamm_results")
  expect_equal(result$diagnostics$performance_mode, "comprehensive")
  expect_equal(result$diagnostics$requested_mode, "full")
})

test_that("simulate_cbamm_data supports modern and legacy signatures", {
  modern <- simulate_cbamm_data(n_rct = 5, n_obs = 3, n_mr = 2, seed = 123)
  legacy <- simulate_cbamm_data(n_studies = 7, seed = 123)

  expect_equal(nrow(modern), 10)
  expect_true(all(c("study_id", "yi", "se", "study_type", "vi", "sei", "ni") %in% names(modern)))
  expect_equal(sort(unique(as.character(modern$study_type))), c("MR", "OBS", "RCT"))

  expect_equal(nrow(legacy), 7)
  expect_true(all(c("study_id", "yi", "se", "study_type", "vi", "sei", "ni") %in% names(legacy)))
  expect_true(all(as.character(legacy$study_type) %in% c("RCT", "OBS")))
})

test_that("detect_study_conflicts handles small datasets with uninformative predictors", {
  conflicted_data <- data.frame(
    study_id = paste0("Study_", 1:10),
    yi = c(rep(0.2, 5), rep(0.8, 5)),
    se = rep(0.1, 10),
    study_type = factor(rep("RCT", 10)),
    analysis_weights = rep(1, 10)
  )

  expect_silent(
    conflicts <- detect_study_conflicts(
      data = conflicted_data,
      threshold = 0.15,
      k_candidates = 2:4
    )
  )
  expect_false(is.null(conflicts))
  expect_equal(conflicts$K, 2)
  expect_true(conflicts$threshold_met)
})
