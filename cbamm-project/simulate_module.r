#' Simulate CBAMM Data for Testing and Examples
#'
#' @param n_rct Number of RCT studies
#' @param n_obs Number of observational studies  
#' @param n_mr Number of Mendelian randomization studies
#' @param true_log_hr True log hazard ratio
#' @param seed Random seed for reproducibility
#' 
#' @return Data frame with simulated study data
#' @export
#'
#' @examples
#' # Default simulation
#' data <- simulate_cbamm_data()
#' 
#' # Custom simulation
#' data <- simulate_cbamm_data(n_rct = 10, n_obs = 15, n_mr = 5)
simulate_cbamm_data <- function(n_rct = 18, n_obs = 18, n_mr = 8,
                                true_log_hr = log(0.90), seed = 42) {
  set.seed(seed)
  
  N <- n_rct + n_obs + n_mr
  study_type <- factor(c(rep("RCT", n_rct), rep("OBS", n_obs), rep("MR", n_mr)),
                       levels = c("RCT", "OBS", "MR"))
  study_id <- sprintf("Study_%02d", seq_len(N))
  
  se <- numeric(N)
  se[study_type == "RCT"] <- runif(sum(study_type == "RCT"), 0.05, 0.15)
  se[study_type == "OBS"] <- runif(sum(study_type == "OBS"), 0.06, 0.18)
  se[study_type == "MR"]  <- runif(sum(study_type == "MR"),  0.10, 0.25)
  
  true_effect <- numeric(N)
  true_effect[study_type == "RCT"] <- rnorm(sum(study_type == "RCT"), true_log_hr,         0.05)
  true_effect[study_type == "OBS"] <- rnorm(sum(study_type == "OBS"), true_log_hr - 0.10,  0.08)
  true_effect[study_type == "MR"]  <- rnorm(sum(study_type == "MR"),  true_log_hr - 0.05,  0.12)
  
  yi <- rnorm(N, true_effect, se)
  
  grade <- character(N)
  grade[study_type == "RCT"] <- sample(c("High", "Moderate"), sum(study_type == "RCT"), TRUE, prob = c(0.7, 0.3))
  grade[study_type == "OBS"] <- sample(c("Moderate", "Low", "Very low"), sum(study_type == "OBS"), TRUE, prob = c(0.3, 0.5, 0.2))
  grade[study_type == "MR"]  <- sample(c("Low", "Very low"), sum(study_type == "MR"),  TRUE, prob = c(0.6, 0.4))
  grade <- factor(grade, levels = c("High", "Moderate", "Low", "Very low"))
  
  neg_ctrl <- numeric(N)
  neg_ctrl[study_type == "OBS"] <- rnorm(sum(study_type == "OBS"), -0.08, 0.04)
  neg_ctrl[study_type != "OBS"] <- rnorm(sum(study_type != "OBS"),  0.00, 0.01)
  
  year <- sample(2010:2025, N, replace = TRUE)
  age_mean   <- round(rnorm(N, 68, 6), 1)
  female_pct <- pmax(0.15, pmin(0.85, rnorm(N, 0.42, 0.10)))
  bmi_mean   <- round(rnorm(N, 28.5, 2.5), 1)
  charlson   <- round(pmax(0, rnorm(N, 1.6, 0.6)), 1)
  
  tte_compliant <- integer(N)
  tte_compliant[study_type == "OBS"] <- rbinom(sum(study_type == "OBS"), 1, 0.65)
  tte_compliant[study_type != "OBS"] <- 1L
  
  data.frame(
    study_id = study_id,
    study_type = study_type,
    yi = yi,
    se = se,
    grade = grade,
    neg_ctrl = neg_ctrl,
    year = year,
    age_mean = age_mean,
    female_pct = female_pct,
    bmi_mean = bmi_mean,
    charlson = charlson,
    tte_compliant = tte_compliant,
    stringsAsFactors = FALSE
  )
}