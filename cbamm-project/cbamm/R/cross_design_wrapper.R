
#' Cross-Design Synthesis Wrapper
#' @export
cross_design_synthesis_unified <- function(data, design_var = "design", ...) {
  if (!is.null(data[[design_var]])) {
    # Split data by design
    rct_data <- data[data[[design_var]] %in% c("RCT", "rct"), ]
    obs_data <- data[data[[design_var]] %in% c("Obs", "Observational", "OBS", "obs"), ]
    
    # Call original function
    return(cross_design_synthesis(rct_data = rct_data, obs_data = obs_data, ...))
  } else {
    # Assume all data is RCT if no design variable
    return(cross_design_synthesis(rct_data = data, obs_data = data.frame(), ...))
  }
}

