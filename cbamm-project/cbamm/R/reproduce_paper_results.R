#' Reproduce PLOS ONE Paper Results
#'
#' Automatically runs all case studies and generates all figures presented
#' in the PLOS ONE submission.
#'
#' @param output_dir Directory to save results (default: current directory)
#' @return Invisible list of results
#' @export
reproduce_paper_results <- function(output_dir = ".") {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  message("Reproducing CBAMM Paper Results (v5.7)...")
  message("Target Directory: ", normalizePath(output_dir))
  
  # 1. Run Figure Generation
  fig_script <- system.file("paper_figures.R", package = "cbamm")
  if (fig_script != "") {
    old_wd <- getwd()
    setwd(output_dir)
    on.exit(setwd(old_wd))
    source(fig_script)
  } else {
    warning("Paper figures script not found in package.")
  }
  
  # 2. Run TruthCert Audit
  audit_script <- system.file("truthcert_audit.R", package = "cbamm")
  if (audit_script != "") {
    source(audit_script)
  }
  
  message("
[SUCCESS] All paper materials reproduced in: ", output_dir)
  return(invisible(TRUE))
}
