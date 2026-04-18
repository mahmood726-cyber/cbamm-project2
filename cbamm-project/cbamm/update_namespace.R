#!/usr/bin/env Rscript
#' Auto-generated script to update CBAMM package exports

# Read current NAMESPACE
namespace <- readLines('NAMESPACE')

# Functions to add
new_exports <- c(
  "network_meta_analysis",
  "validate_network_inputs",
  "validate_network_structure",
  "is_network_connected",
  "prepare_network_data",
  "run_frequentist_network_ma",
  "extract_nma_results",
  "run_basic_nma",
  "assess_network_inconsistency",
  "run_bayesian_network_ma",
  "print.cbamm_network",
  "summary.cbamm_network",
  "plot.cbamm_network",
  "plot_forest_network",
  "plot_network_diagram",
  "plot_network_enhanced",
  "plot_nma_forest",
  "plot_funnel_enhanced",
  "plot_influence_analysis",
  "plot_outlier_detection"
)

# Add export statements
for (func in new_exports) {
  export_line <- sprintf('export(%s)', func)
  if (!any(grepl(export_line, namespace, fixed = TRUE))) {
    namespace <- c(namespace, export_line)
    cat('Added:', export_line, '\n')
  }
}

# Write updated NAMESPACE
writeLines(namespace, 'NAMESPACE')
cat('\nNAMESPACE updated successfully!\n')
cat('Run devtools::document() to update documentation\n')
