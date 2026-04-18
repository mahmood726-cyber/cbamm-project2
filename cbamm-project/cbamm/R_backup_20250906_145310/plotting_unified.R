
# CBAMM Unified Plotting Functions
# Using established packages: metafor and meta

cbamm_forest <- function(data) {
  require(metafor)
  ma <- rma(yi = data$effect_size, sei = data$se, slab = data$study)
  forest(ma, main = "Forest Plot", xlab = "Effect Size")
  return(ma)
}

cbamm_funnel <- function(data, type = "standard") {
  require(metafor)
  ma <- rma(yi = data$effect_size, sei = data$se)
  
  if (type == "contour") {
    funnel(ma, level = c(90, 95, 99), 
           shade = c("white", "gray55", "gray75"),
           main = "Contour-Enhanced Funnel Plot")
  } else if (type == "trimfill") {
    tf <- trimfill(ma)
    funnel(tf, main = "Trim-and-Fill Funnel Plot")
    return(tf)
  } else {
    funnel(ma, main = "Funnel Plot")
  }
  return(ma)
}

cbamm_diagnostic <- function(data) {
  require(metafor)
  ma <- rma(yi = data$effect_size, sei = data$se, slab = data$study)
  par(mfrow = c(2, 2))
  forest(ma, main = "Forest", cex = 0.8)
  funnel(ma, main = "Funnel")
  baujat(ma, main = "Baujat")
  radial(ma, main = "Radial")
  par(mfrow = c(1, 1))
  return(ma)
}

