
#' Launch CBAMM Dashboard
#' 
#' Launches an interactive Shiny dashboard for meta-analysis
#' 
#' @export
#' @examples
#' \dontrun{
#' launch_cbamm_dashboard()
#' }
launch_cbamm_dashboard <- function() {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required. Please install it.")
  }
  if (!requireNamespace("shinydashboard", quietly = TRUE)) {
    stop("Package 'shinydashboard' is required. Please install it.")
  }
  
  app_dir <- system.file("shiny", package = "cbamm")
  if (app_dir == "") {
    stop("Could not find Shiny app directory. Try reinstalling the package.")
  }
  
  shiny::runApp(app_dir, launch.browser = TRUE)
}

