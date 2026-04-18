
# CBAMM Plotting Functions Test Script

cat("CBAMM Plotting Functions Test\n")
cat("============================\n\n")

# Load required libraries
library(ggplot2)
library(scales)

# Source the fixed functions
source("R/fixed_plotting_functions.R")

# Create comprehensive test data
create_test_data <- function() {
  set.seed(42)
  data.frame(
    study = c("Smith 2020", "Jones 2021", "Brown 2019", 
              "Davis 2022", "Wilson 2021", "Taylor 2020"),
    effect_size = c(0.25, 0.45, 0.15, 0.65, 0.35, 0.55),
    se = c(0.12, 0.15, 0.18, 0.14, 0.16, 0.13)
  )
}

# Run all tests
test_data <- create_test_data()

cat("Testing all fixed plotting functions:\n\n")

# Test each function
tryCatch({
  cat("1. Forest Plot Classic... ")
  forest_plot_classic(test_data)
  cat("SUCCESS\n")
}, error = function(e) cat("ERROR:", e$message, "\n"))

tryCatch({
  cat("2. Funnel Plot Classic... ")
  funnel_plot_classic(test_data)
  cat("SUCCESS\n")
}, error = function(e) cat("ERROR:", e$message, "\n"))

tryCatch({
  cat("3. Funnel Plot Contour... ")
  funnel_plot_contour(test_data)
  cat("SUCCESS\n")
}, error = function(e) cat("ERROR:", e$message, "\n"))

tryCatch({
  cat("4. Funnel Plot Trim Fill... ")
  funnel_plot_trim_fill(test_data)
  cat("SUCCESS\n")
}, error = function(e) cat("ERROR:", e$message, "\n"))

tryCatch({
  cat("5. P-Curve Plot... ")
  p_curve_plot(test_data)
  cat("SUCCESS\n")
}, error = function(e) cat("ERROR:", e$message, "\n"))

cat("\nAll tests completed!\n")

