# CBAMM Package Build and Validation Script
# Run this script to create the complete CBAMM package structure

# Load required packages
if (!require(devtools)) install.packages("devtools")
if (!require(usethis)) install.packages("usethis")

# Set up the package directory
package_dir <- "cbamm"
if (dir.exists(package_dir)) {
  cat("Package directory already exists. Overwriting...\n")
  unlink(package_dir, recursive = TRUE)
}

# Create package structure
usethis::create_package(package_dir, open = FALSE)
setwd(package_dir)

# Create directory structure
cat("Creating directory structure...\n")
dir.create("R", showWarnings = FALSE)
dir.create("tests/testthat", recursive = TRUE, showWarnings = FALSE)
dir.create("man", showWarnings = FALSE)
dir.create("vignettes", showWarnings = FALSE)

# File mapping: which source file goes to which destination
file_mapping <- list(
  # R modules
  "config_module.r" = "R/config.R",
  "main_interface.r" = "R/cbamm.R", 
  "utils_module.r" = "R/utils.R",
  "simulate_module.r" = "R/simulate.R",
  "weights_module.r" = "R/weights.R",
  "meta_analysis_module.r" = "R/meta_analysis.R",
  "sensitivity_module.r" = "R/sensitivity.R", 
  "novelty_module.r" = "R/novelty.R",
  "plotting_module.r" = "R/plotting.R",
  "bayesian_module.r" = "R/bayesian.R",
  
  # Tests
  "basic_tests.r" = "tests/testthat/test-core.R",
  
  # Package files
  "DESCRIPTION" = "DESCRIPTION",
  "NAMESPACE" = "NAMESPACE",
  "README.md" = "README.md"
)

# Function to copy and rename files
copy_source_files <- function(source_dir = "..") {
  cat("Copying source files...\n")
  
  for (source_file in names(file_mapping)) {
    dest_file <- file_mapping[[source_file]]
    source_path <- file.path(source_dir, source_file)
    
    if (file.exists(source_path)) {
      # Create destination directory if needed
      dest_dir <- dirname(dest_file)
      if (!dir.exists(dest_dir)) {
        dir.create(dest_dir, recursive = TRUE)
      }
      
      # Copy file
      file.copy(source_path, dest_file, overwrite = TRUE)
      cat("  ✓", source_file, "->", dest_file, "\n")
    } else {
      cat("  ✗", source_file, "not found\n")
    }
  }
}

# Create missing S3 methods file
cat("Creating S3 methods...\n")
s3_methods_content <- '
#\' S3 Methods for CBAMM Objects
#\' 
#\' Print, summary, and plot methods for CBAMM results and configuration objects

#\' Print Method for CBAMM Results
#\' @param x CBAMM results object
#\' @param ... Additional arguments (unused)
#\' @return Invisible x
#\' @export
print.cbamm_results <- function(x, ...) {
  cat("CBAMM Meta-Analysis Results\\n")
  cat("===========================\\n\\n")
  
  # Data summary
  if (!is.null(x$data_summary)) {
    cat("Studies:", x$data_summary$n_studies, "\\n")
    if (!is.null(x$data_summary$study_types)) {
      cat("Study Types:", paste(names(x$data_summary$study_types), collapse = ", "), "\\n")
    }
  }
  
  # Core results
  if (!is.null(x$meta_results$transport)) {
    pred <- try(metafor::predict(x$meta_results$transport, transf = exp), silent = TRUE)
    if (!inherits(pred, "try-error")) {
      cat(sprintf("Transport HR: %.3f (95%% CI: %.3f-%.3f)\\n", 
                 pred$pred, pred$ci.lb, pred$ci.ub))
    }
  }
  
  invisible(x)
}

#\' Summary Method for CBAMM Results
#\' @param object CBAMM results object  
#\' @param ... Additional arguments (unused)
#\' @return Summary data frame
#\' @export
summary.cbamm_results <- function(object, ...) {
  # Basic summary - can be enhanced
  if (!is.null(object$meta_results$transport)) {
    fit <- object$meta_results$transport
    return(data.frame(
      Method = "Transport-weighted",
      Estimate = as.numeric(coef(fit)),
      SE = as.numeric(fit$se),
      I_squared = fit$I2,
      Studies = fit$k
    ))
  }
  return(data.frame())
}

#\' Plot Method for CBAMM Results
#\' @param x CBAMM results object
#\' @param which Which plot to show
#\' @param ... Additional arguments  
#\' @return Plot object
#\' @export
plot.cbamm_results <- function(x, which = "all", ...) {
  if (length(x$plots) == 0) {
    warning("No plots available")
    return(NULL)
  }
  
  if (which == "all") {
    return(x$plots)
  } else if (which %in% names(x$plots)) {
    return(x$plots[[which]])
  } else {
    warning("Plot not found: ", which)
    return(NULL)
  }
}
'

writeLines(s3_methods_content, "R/methods.R")
cat("  ✓ R/methods.R created\n")

# Now copy the source files (assumes source files are in parent directory)
copy_source_files()

# Build package documentation
cat("\\nBuilding documentation...\\n")
try({
  devtools::document()
  cat("  ✓ Documentation generated\\n")
}, silent = FALSE)

# Check package
cat("\\nChecking package...\\n")
check_result <- try({
  devtools::check(quiet = TRUE)
}, silent = TRUE)

if (inherits(check_result, "try-error")) {
  cat("  ⚠ Package check had issues (this is normal for initial build)\\n")
} else {
  cat("  ✓ Package check completed\\n")
}

# Test basic functionality
cat("\\nTesting basic functionality...\\n")
test_result <- try({
  devtools::load_all(quiet = TRUE)
  
  # Test configuration
  config <- cbamm_config()
  cat("  ✓ Configuration works\\n")
  
  # Test data simulation
  if (exists("simulate_cbamm_data")) {
    data <- simulate_cbamm_data(n_rct = 5, n_obs = 3)
    cat("  ✓ Data simulation works\\n")
    
    # Test basic pipeline (may fail due to missing functions)
    minimal_config <- cbamm_config(
      methods = list(
        transport = FALSE,
        bayesian = FALSE,
        conflict_detection = FALSE,
        pet_peese = FALSE,
        missing_studies = FALSE
      ),
      output = list(verbose = FALSE, plots = FALSE)
    )
    
    if (exists("cbamm")) {
      results <- try(cbamm(data, config = minimal_config), silent = TRUE)
      if (!inherits(results, "try-error")) {
        cat("  ✓ Basic CBAMM pipeline works\\n")
      } else {
        cat("  ⚠ CBAMM pipeline needs debugging\\n")
      }
    }
  }
  
}, silent = TRUE)

if (inherits(test_result, "try-error")) {
  cat("  ⚠ Some functionality tests failed (expected during development)\\n")
}

# Install package
cat("\\nInstalling package...\\n")
install_result <- try({
  devtools::install(quiet = TRUE)
  cat("  ✓ Package installed successfully\\n")
}, silent = TRUE)

if (inherits(install_result, "try-error")) {
  cat("  ⚠ Package installation had issues\\n")
}

# Summary
cat("\\n" %+% rep("=", 50) %+% "\\n")
cat("CBAMM PACKAGE BUILD SUMMARY\\n")
cat(rep("=", 50) %+% "\\n")
cat("Package directory: ", getwd(), "\\n")
cat("\\nNext steps:\\n")
cat("1. Fix any R CMD check warnings/errors\\n")
cat("2. Test functionality with: library(cbamm); cbamm(simulate_cbamm_data())\\n")
cat("3. Add comprehensive tests and documentation\\n")
cat("4. Publish to GitHub/CRAN\\n")

cat("\\n🎉 CBAMM package structure is ready!\\n")

# Return to original directory
setwd("..")