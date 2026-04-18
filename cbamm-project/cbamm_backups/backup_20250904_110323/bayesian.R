#' Run Bayesian Analysis with brms and JAGS Fallback
#'
#' @param data Study data
#' @param config CBAMM configuration
#'
#' @return Bayesian analysis results
#' @keywords internal
run_bayesian_analysis <- function(data, config) {
  if (!config$methods$bayesian) {
    return(NULL)
  }
  
  # Prepare weights
  w_ext <- if ("analysis_weights_grade" %in% names(data)) {
    dplyr::coalesce(data$analysis_weights_grade, data$analysis_weights)
  } else {
    data$analysis_weights
  }
  
  # Try brms first if available
  brms_result <- .try_brms_analysis(data, w_ext, config)
  if (!is.null(brms_result)) {
    return(brms_result)
  }
  
  # Fall back to JAGS
  if (config$output$verbose) message("brms failed, trying JAGS fallback...")
  jags_result <- .try_jags_analysis(data, w_ext, config)
  
  return(jags_result)
}

#' Try brms Analysis
#' @keywords internal
.try_brms_analysis <- function(data, w_ext, config) {
  if (!check_package("brms", quietly = TRUE) || 
      !check_package("posterior", quietly = TRUE) || 
      !check_package("loo", quietly = TRUE)) {
    return(NULL)
  }
  
  with_fallback(
    pkg = "brms",
    fun_call = {
      # Prepare data
      bayes_data <- data %>% 
        dplyr::mutate(
          study_type_fct = factor(study_type, levels = c("RCT","OBS","MR")),
          grade_numeric  = as.numeric(grade),
          w = w_ext
        )
      
      # Set up priors
      priors1 <- if (!is.null(config$bayesian$priors)) {
        config$bayesian$priors 
      } else {
        c(
          brms::prior(normal(0, 0.3), class = Intercept),
          brms::prior(normal(0, 0.2), class = b),
          brms::prior(student_t(3, 0, 0.3), class = sd)
        )
      }
      
      priors2 <- if (!is.null(config$bayesian$priors_alt)) {
        config$bayesian$priors_alt
      } else {
        c(
          brms::prior(normal(0, 0.5), class = Intercept),
          brms::prior(normal(0, 0.3), class = b),
          brms::prior(student_t(3, 0, 0.5), class = sd)
        )
      }
      
      # Fit models
      m1 <- try(brms::brm(
        yi | se(se) ~ 1 + study_type_fct + grade_numeric + (1|study_id),
        data = bayes_data, 
        family = gaussian(),
        prior = priors1,
        chains = config$bayesian$chains, 
        iter = config$bayesian$iter, 
        warmup = config$bayesian$warmup, 
        cores = min(config$bayesian$cores, config$n_cores %||% 1),
        control = list(adapt_delta = 0.98), 
        seed = 1001,
        weights = w_ext,
        save_pars = brms::save_pars(all = TRUE),
        refresh = 0
      ), silent = TRUE)
      
      m2 <- try(brms::brm(
        yi | se(se) ~ 1 + study_type_fct + grade_numeric + (1|study_id),
        data = bayes_data, 
        family = gaussian(),
        prior = priors2,
        chains = config$bayesian$chains, 
        iter = config$bayesian$iter, 
        warmup = config$bayesian$warmup, 
        cores = min(config$bayesian$cores, config$n_cores %||% 1),
        control = list(adapt_delta = 0.98), 
        seed = 1002,
        weights = w_ext,
        save_pars = brms::save_pars(all = TRUE),
        refresh = 0
      ), silent = TRUE)
      
      # Process successful fits
      fits <- list()
      if (!inherits(m1, "try-error")) {
        fits$cons <- m1
        .check_brms_convergence(m1, "model 1")
      }
      if (!inherits(m2, "try-error")) {
        fits$mod <- m2
        .check_brms_convergence(m2, "model 2")
      }
      
      if (length(fits) == 0) return(NULL)
      
      # Model stacking
      stacking_result <- .perform_model_stacking(fits)
      
      return(list(
        engine = "brms",
        weights = stacking_result$weights,
        draws = stacking_result$mixed_draws,
        summary = stacking_result$summary,
        models = fits
      ))
    },
    fallback_fun = function() NULL,
    fallback_msg = "brms analysis failed"
  )
}

#' Check brms Convergence
#' @keywords internal
.check_brms_convergence <- function(model, label) {
  rhats <- try(brms::rhat(model), silent = TRUE)
  if (!inherits(rhats, "try-error") && any(rhats > 1.05, na.rm = TRUE)) {
    warning("High R-hat values (>1.05) detected in brms ", label, ". Check convergence.")
  }
}

#' Perform Model Stacking
#' @keywords internal
.perform_model_stacking <- function(fits) {
  # Compute LOO for each model
  llist <- lapply(fits, function(f) {
    loo_result <- try(loo::loo(f, pointwise = TRUE), silent = TRUE)
    if (inherits(loo_result, "try-error")) return(NULL)
    loo_result
  })
  
  # Remove failed LOO computations
  llist <- llist[!sapply(llist, is.null)]
  
  # Compute stacking weights
  if (length(llist) > 1 && all(sapply(llist, function(x) inherits(x, "loo")))) {
    sw <- loo::loo_model_weights(llist, method = "stacking")
  } else {
    sw <- rep(1/length(fits), length(fits))
    names(sw) <- names(fits)
  }
  
  # Extract and mix draws
  draws <- lapply(fits, function(f) posterior::as_draws_df(f)[,"b_Intercept"])
  M <- min(sapply(draws, length))
  set.seed(2024)
  
  mix <- c()
  for (i in seq_along(draws)) {
    take <- max(1, round(M * sw[i]))
    mix <- c(mix, draws[[i]][sample(seq_len(length(draws[[i]])), size = take, replace = TRUE)])
  }
  
  hr <- exp(mix)
  
  if (length(sw) > 1) {
    message("Stacking weights: ", paste(sprintf("%s=%.2f", names(sw), sw), collapse = ", "))
  }
  message(sprintf("Bayesian (brms) baseline HR: %.3f (95%% CrI %.3f–%.3f)",
              median(hr), quantile(hr, 0.025), quantile(hr, 0.975)))
  
  list(
    weights = sw,
    mixed_draws = mix,
    summary = list(
      median_hr = median(hr), 
      cri = c(quantile(hr, c(0.025, 0.975)))
    )
  )
}

#' Try JAGS Analysis
#' @keywords internal
.try_jags_analysis <- function(data, w_ext, config) {
  if (!check_package("rjags", quietly = TRUE)) {
    return(NULL)
  }
  
  with_fallback(
    pkg = "rjags",
    fun_call = {
      # Prepare data for JAGS
      dd <- data %>% 
        dplyr::mutate(
          is_obs = as.integer(study_type == "OBS"),
          is_mr  = as.integer(study_type == "MR"),
          grade_num = as.numeric(grade),
          study_idx = as.integer(factor(study_id)),
          w = w_ext
        )
      
      jdat <- list(
        N = nrow(dd),
        S = max(dd$study_idx),
        yi = dd$yi,
        se = dd$se,
        w  = as.numeric(dd$w),
        is_obs = dd$is_obs,
        is_mr  = dd$is_mr,
        grade  = dd$grade_num,
        study = dd$study_idx
      )
      
      # JAGS model
      model_string <- "
      model{
        for(i in 1:N){
          prec[i] <- w[i] / pow(se[i], 2)
          yi[i] ~ dnorm(eta[i], prec[i])
          eta[i] <- mu + b_obs * is_obs[i] + b_mr * is_mr[i] + b_grade * grade[i] + u[study[i]]
        }
        for(s in 1:S){
          u[s] ~ dnorm(0, tau_u)
        }
        mu ~ dnorm(0, pow(0.5, -2))
        b_obs ~ dnorm(0, pow(0.3, -2))
        b_mr  ~ dnorm(0, pow(0.3, -2))
        b_grade ~ dnorm(0, pow(0.3, -2))
        sigma_u ~ dunif(0, 1)
        tau_u <- pow(sigma_u, -2)
      }"
      
      # Fit JAGS model
      jm <- rjags::jags.model(
        textConnection(model_string), 
        data = jdat, 
        n.chains = config$bayesian$chains, 
        quiet = TRUE
      )
      
      # Use S3 update method (not rjags::update)
      update(jm, config$bayesian$warmup, progress.bar = "none")
      
      sm <- rjags::coda.samples(
        jm, 
        c("mu", "b_obs", "b_mr", "b_grade", "sigma_u"),
        n.iter = max(1, config$bayesian$iter - config$bayesian$warmup),
        thin = 2, 
        progress.bar = "none"
      )
      
      # Convert to matrix (handle S3 method dispatch)
      draws <- try(as.matrix(sm), silent = TRUE)
      if (inherits(draws, "try-error")) {
        draws <- try(coda::as.matrix.mcmc.list(sm), silent = TRUE)
      }
      
      if (inherits(draws, "try-error")) {
        warning("Could not convert JAGS samples to matrix")
        return(list(
          engine = "jags", 
          draws_mu = NULL,
          summary = NULL
        ))
      }
      
      hr <- exp(draws[,"mu"])
      message(sprintf("Bayesian (JAGS) baseline HR: %.3f (95%% CrI %.3f–%.3f)",
                  median(hr), quantile(hr, 0.025), quantile(hr, 0.975)))
      
      list(
        engine = "jags", 
        draws_mu = draws[,"mu"],
        summary = list(
          median_hr = median(hr), 
          cri = c(quantile(hr, c(0.025, 0.975)))
        )
      )
    },
    fallback_fun = function() NULL,
    fallback_msg = "JAGS analysis failed"
  )
}