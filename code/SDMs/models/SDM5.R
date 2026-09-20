################################################################################
##                                                                            ##
##                           SDM.5 FITTED WITH brms                           ##
##                            (SDMs/models/SDM5.R)                            ##
##                                                                            ##
##       Addressing Missing Covariates in Species Distribution Models:        ##
##  Inferential Impacts and Mitigation via Joint Species Distribution Models  ##
##                                                                            ##
##                  Arthur F. Rossignol & Frédéric Gosselin                   ##
##                                                                            ##
##                                    2026                                    ##
##                                                                            ##
################################################################################

fit_SDM.5 <- function(Y,
                      X1, 
                      model_type,
                      use_OLRE,
                      seed) {
  
  obs_id  <- factor(seq_along(Y))
  data_df <- data.frame(Y = Y, X1 = X1, obs_id = obs_id)
  
  get_family_brms <- function(model_type) {
    switch(model_type,
           "gaussian"          = gaussian(),
           "poisson"           = poisson(link = "log"),
           "bernoulli_probit"  = bernoulli(link = "probit"),
           "bernoulli_logit"   = bernoulli(link = "logit"),
           "bernoulli_cloglog" = bernoulli(link = "cloglog"))
  }
  
  family <- get_family_brms(model_type)
  
  if (use_OLRE) {
    formula_brms <- bf(Y ~ X1 + (1 | obs_id))
  } else {
    formula_brms <- bf(Y ~ X1)
  }
  
  prior_spec <- c(prior(normal(0, 10), class = Intercept),
                  prior(normal(0, 10), class = b))
  
  if (use_OLRE) {
    prior_spec <- c(prior_spec, prior(cauchy(0, 1), class = sd))
  }
  
  fit <- brm(formula = formula_brms,
             data    = data_df,
             family  = family,
             prior   = prior_spec,
             chains  = 3,
             iter    = 40000,
             warmup  = 20000,
             thin    = 5,
             cores   = 3,
             control = list(adapt_delta = 0.95, max_treedepth = 12),
             seed    = seed)
  
  post_summary <- posterior_summary(fit, pars = c("b_Intercept", "b_X1"))
  
  result <- list(alpha_est = post_summary["b_Intercept", "Estimate"],
                 alpha_SE  = post_summary["b_Intercept", "Est.Error"],
                 beta_est  = post_summary["b_X1", "Estimate"],
                 beta_SE   = post_summary["b_X1", "Est.Error"],
                 converged = NA)
  
  rhat_vals <- rhat(fit)
  result$converged <- !any(rhat_vals > 1.1, na.rm = TRUE)
  
  if (use_OLRE) {
    re_summary <- VarCorr(fit, summary = TRUE)
    sd_site <- re_summary$obs_id$sd[1, "Estimate"]
    result$sigma2_est <- sd_site^2
  } else if (model_type == "gaussian") {
    sigma_est <- posterior_summary(fit, pars = "sigma")["sigma", "Estimate"]
    result$sigma2_est <- sigma_est^2
  }
  
  return(result)
}

################################################################################
