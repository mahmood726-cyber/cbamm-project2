simulate_cbamm_data <- function(n_rct = NULL,
                                n_obs = NULL,
                                n_mr = NULL,
                                true_log_hr = log(0.90),
                                seed = 42,
                                n_studies = NULL,
                                n_min = 20,
                                n_max = 200,
                                true_effect = NULL,
                                tau = 0.2) {
  set.seed(seed)

  using_legacy_signature <- !is.null(n_studies) &&
    is.null(n_rct) && is.null(n_obs) && is.null(n_mr)

  if (using_legacy_signature) {
    n_rct <- floor(n_studies / 2)
    n_obs <- n_studies - n_rct
    n_mr <- 0L
  } else {
    if (is.null(n_rct)) n_rct <- 18L
    if (is.null(n_obs)) n_obs <- 18L
    if (is.null(n_mr)) n_mr <- 8L
  }

  if (!is.null(true_effect)) {
    true_log_hr <- true_effect
  }

  n_rct <- as.integer(n_rct)
  n_obs <- as.integer(n_obs)
  n_mr <- as.integer(n_mr)

  N <- n_rct + n_obs + n_mr
  study_type <- factor(
    c(rep("RCT", n_rct), rep("OBS", n_obs), rep("MR", n_mr)),
    levels = c("RCT", "OBS", "MR")
  )
  study_id <- sprintf("Study_%02d", seq_len(N))

  ni <- round(runif(N, n_min, n_max))
  se <- numeric(N)
  se[study_type == "RCT"] <- runif(sum(study_type == "RCT"), 0.05, 0.15)
  se[study_type == "OBS"] <- runif(sum(study_type == "OBS"), 0.06, 0.18)
  se[study_type == "MR"] <- runif(sum(study_type == "MR"), 0.10, 0.25)

  heterogeneity <- if (using_legacy_signature) tau else 0.05
  true_effect_vec <- numeric(N)
  true_effect_vec[study_type == "RCT"] <- rnorm(sum(study_type == "RCT"), true_log_hr, heterogeneity)
  true_effect_vec[study_type == "OBS"] <- rnorm(sum(study_type == "OBS"), true_log_hr - 0.10, max(heterogeneity, 0.08))
  true_effect_vec[study_type == "MR"] <- rnorm(sum(study_type == "MR"), true_log_hr - 0.05, max(heterogeneity, 0.12))

  yi <- rnorm(N, true_effect_vec, se)
  vi <- se^2

  grade <- character(N)
  grade[study_type == "RCT"] <- sample(c("High", "Moderate"), sum(study_type == "RCT"), TRUE, prob = c(0.7, 0.3))
  grade[study_type == "OBS"] <- sample(c("Moderate", "Low", "Very low"), sum(study_type == "OBS"), TRUE, prob = c(0.3, 0.5, 0.2))
  grade[study_type == "MR"] <- sample(c("Low", "Very low"), sum(study_type == "MR"), TRUE, prob = c(0.6, 0.4))
  grade <- factor(grade, levels = c("High", "Moderate", "Low", "Very low"))

  neg_ctrl <- numeric(N)
  neg_ctrl[study_type == "OBS"] <- rnorm(sum(study_type == "OBS"), -0.08, 0.04)
  neg_ctrl[study_type != "OBS"] <- rnorm(sum(study_type != "OBS"), 0.00, 0.01)

  tte_compliant <- integer(N)
  tte_compliant[study_type == "OBS"] <- rbinom(sum(study_type == "OBS"), 1, 0.65)
  tte_compliant[study_type != "OBS"] <- 1L

  data.frame(
    study_id = study_id,
    study_type = study_type,
    yi = yi,
    vi = vi,
    sei = se,
    se = se,
    ni = ni,
    n = ni,
    grade = grade,
    neg_ctrl = neg_ctrl,
    year = sample(2010:2025, N, replace = TRUE),
    region = factor(sample(c("US", "EU", "Asia"), N, replace = TRUE)),
    age_mean = round(rnorm(N, 68, 6), 1),
    female_pct = pmax(0.15, pmin(0.85, rnorm(N, 0.42, 0.10))),
    bmi_mean = round(rnorm(N, 28.5, 2.5), 1),
    charlson = round(pmax(0, rnorm(N, 1.6, 0.6)), 1),
    tte_compliant = tte_compliant,
    stringsAsFactors = FALSE
  )
}
