
#' Example Network Meta-Analysis Dataset
#'
#' A simulated dataset for demonstrating network meta-analysis functionality.
#' Contains effect sizes and standard errors for multiple treatments across studies.
#'
#' @format A data frame with variables:
#' \describe{
#'   \item{study_id}{Study identifier}
#'   \item{treatment}{Treatment name (A, B, C, D)}
#'   \item{effect}{Effect size (log odds ratio)}
#'   \item{se}{Standard error of effect size}
#'   \item{sample_size}{Sample size per arm}
#' }
#' @examples
#' data(network_example_data)
#' head(network_example_data)
"network_example_data"

# Generate example network data
set.seed(12345)

# Create a connected network with 4 treatments
studies <- paste0("Study_", 1:15)
treatments <- c("Placebo", "Treatment_A", "Treatment_B", "Treatment_C")

# Generate study-treatment combinations ensuring connectivity
network_data_list <- list()

# Studies with Placebo vs Treatment_A
for(i in 1:5) {
  network_data_list[[length(network_data_list) + 1]] <- data.frame(
    study_id = studies[i],
    treatment = "Placebo",
    effect = rnorm(1, 0, 0.1),
    se = runif(1, 0.15, 0.25),
    sample_size = sample(50:200, 1)
  )
  
  network_data_list[[length(network_data_list) + 1]] <- data.frame(
    study_id = studies[i],
    treatment = "Treatment_A", 
    effect = rnorm(1, 0.3, 0.15),
    se = runif(1, 0.15, 0.25),
    sample_size = sample(50:200, 1)
  )
}

# Studies with Placebo vs Treatment_B
for(i in 6:10) {
  network_data_list[[length(network_data_list) + 1]] <- data.frame(
    study_id = studies[i],
    treatment = "Placebo",
    effect = rnorm(1, 0, 0.1),
    se = runif(1, 0.15, 0.25),
    sample_size = sample(50:200, 1)
  )
  
  network_data_list[[length(network_data_list) + 1]] <- data.frame(
    study_id = studies[i],
    treatment = "Treatment_B",
    effect = rnorm(1, 0.25, 0.15),
    se = runif(1, 0.15, 0.25),
    sample_size = sample(50:200, 1)
  )
}

# Studies with Treatment_A vs Treatment_C  
for(i in 11:13) {
  network_data_list[[length(network_data_list) + 1]] <- data.frame(
    study_id = studies[i],
    treatment = "Treatment_A",
    effect = rnorm(1, 0.3, 0.15),
    se = runif(1, 0.15, 0.25),
    sample_size = sample(50:200, 1)
  )
  
  network_data_list[[length(network_data_list) + 1]] <- data.frame(
    study_id = studies[i],
    treatment = "Treatment_C",
    effect = rnorm(1, 0.4, 0.15),
    se = runif(1, 0.15, 0.25),
    sample_size = sample(50:200, 1)
  )
}

# Multi-arm studies
for(i in 14:15) {
  for(trt in c("Placebo", "Treatment_A", "Treatment_B")) {
    effect_val <- switch(trt,
                        "Placebo" = rnorm(1, 0, 0.1),
                        "Treatment_A" = rnorm(1, 0.3, 0.15),
                        "Treatment_B" = rnorm(1, 0.25, 0.15))
    
    network_data_list[[length(network_data_list) + 1]] <- data.frame(
      study_id = studies[i],
      treatment = trt,
      effect = effect_val,
      se = runif(1, 0.15, 0.25),
      sample_size = sample(50:200, 1)
    )
  }
}

network_example_data <- do.call(rbind, network_data_list)
rownames(network_example_data) <- NULL

# Save to data directory
if (!dir.exists("data")) {
  dir.create("data")
}

save(network_example_data, file = "data/network_example_data.rda")

