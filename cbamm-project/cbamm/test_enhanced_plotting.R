
# Enhanced CBAMM Test Script with Advanced Plotting
# Tests all cumulative meta-analysis and plotting functionality

cat("Enhanced CBAMM Package Testing\n")
cat("=============================\n\n")

# Load package
devtools::load_all()

# Load example data
load("data/cumulative_example_data.rda")

# Run cumulative meta-analysis
result <- cumulative_meta_analysis(cumulative_example_data, order_by = "year")

cat("Basic Analysis Results:\n")
print(result)
cat("\n")

# Test enhanced plotting functions
cat("Testing Enhanced Plotting Functions:\n")
cat("-----------------------------------\n")

if (requireNamespace("ggplot2", quietly = TRUE)) {
  
  cat("1. Creating comprehensive dashboard...\n")
  create_cumulative_dashboard(result)
  
  cat("2. Creating evidence trajectory plot...\n")
  plot_evidence_trajectory(result, highlight_studies = c(10, 15, 20))
  
  cat("3. Creating stability assessment plot...\n")
  plot_stability_assessment(result)
  
  cat("4. Creating forest plot evolution...\n")
  plot_forest_evolution(result, steps_to_show = c(1, 5, 10, 15, 20, 24))
  
  cat("5. Testing save functionality...\n")
  create_cumulative_dashboard(result, save_plot = TRUE, filename = "test_dashboard.png")
  
  cat("\n✓ All enhanced plotting functions working correctly\n")
  
} else {
  cat("ggplot2 not available - install for plotting functionality\n")
}

cat("\nPackage Enhancement Complete!\n")
cat("============================\n\n")

cat("Available Functions:\n")
cat("- cumulative_meta_analysis(): Core analysis function\n")
cat("- plot(): Basic plotting (estimate, ci, precision, heterogeneity)\n")
cat("- create_cumulative_dashboard(): Comprehensive 4-panel dashboard\n")
cat("- plot_evidence_trajectory(): Advanced trajectory visualization\n")
cat("- plot_stability_assessment(): Stability metrics over time\n")
cat("- plot_forest_evolution(): Evolving forest plots\n")

