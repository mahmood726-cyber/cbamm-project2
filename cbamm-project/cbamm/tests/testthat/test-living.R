context("Living Systematic Review")

test_that("Living review framework works", {
  # Generate test data
  initial_data <- simulate_sequential_data(n_studies = 5)
  
  # Create living review
  result <- living_systematic_review(
    initial_data = initial_data,
    effect_size = "effect",
    variance = "variance",
    study_id = "study_id",
    method = "sequential",
    boundaries = "obrien-fleming"
  )
  
  expect_s3_class(result, "cbamm_living")
  expect_true(!is.null(result$sequential$cumulative_z))
  
  # Test updating
  new_data <- simulate_sequential_data(n_studies = 1, start_id = 6)
  updated <- update_living_review(result, new_data)
  
  expect_true(nrow(updated$data) > nrow(result$data))
})

