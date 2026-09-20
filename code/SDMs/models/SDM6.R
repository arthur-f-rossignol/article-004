################################################################################
##                                                                            ##
##                           SDM.6 FITTED WITH INLA                           ##
##                            (SDMs/models/SDM6.R)                            ##
##                                                                            ##
##       Addressing Missing Covariates in Species Distribution Models:        ##
##  Inferential Impacts and Mitigation via Joint Species Distribution Models  ##
##                                                                            ##
##                  Arthur F. Rossignol & Frédéric Gosselin                   ##
##                                                                            ##
##                                    2026                                    ##
##                                                                            ##
################################################################################

fit_SDM.6 <- function(Y, 
                      X1, 
                      model_type, 
                      use_OLRE) {
  
  obs_id  <- factor(seq_along(Y))
  data_df <- data.frame(Y = Y, X1 = X1, obs_id = obs_id)
  
  get_family_inla <- function(model_type) {
    switch(model_type,
           "gaussian"          = "gaussian",
           "poisson"           = "poisson",
           "bernoulli_probit"  = "binomial",
           "bernoulli_logit"   = "binomial",
           "bernoulli_cloglog" = "binomial",
           NULL)
  }
  
  get_link_inla <- function(model_type) {
    switch(model_type,
           "gaussian"          = "identity",
           "poisson"           = "log",
           "bernoulli_probit"  = "probit",
           "bernoulli_logit"   = "logit",
           "bernoulli_cloglog" = "cloglog",
           "identity")
  }
  
  family <- get_family_inla(model_type)
  link_spec   <- get_link_inla(model_type)

  if (use_OLRE) {
    formula_inla <- Y ~ X1 + f(obs_id, model = "iid")
  } else {
    formula_inla <- Y ~ X1
  }

  # INLA binomial requires explicit Ntrials
  Ntrials <- if (family == "binomial") rep(1L, length(Y)) else NULL

  fit <- inla(formula           = formula_inla,
              family            = family,
              data              = data_df,
              Ntrials           = Ntrials,
              control.family    = list(link = link_spec),
              control.compute   = list(config = TRUE),
              control.predictor = list(compute = TRUE),
              verbose           = FALSE)
  
  fixed_effects <- fit$summary.fixed
  
  result <- list(alpha_est = fixed_effects["(Intercept)", "mean"],
                 alpha_SE  = fixed_effects["(Intercept)", "sd"],
                 beta_est  = fixed_effects["X1", "mean"],
                 beta_SE   = fixed_effects["X1", "sd"],
                 converged = TRUE)
  
  if (!is.null(fit$mode$mode.status) && fit$mode$mode.status != 0) {
    result$converged <- FALSE
  }
  
  if (use_OLRE && !is.null(fit$summary.hyperpar)) {
    hyper_summary     <- fit$summary.hyperpar
    precision_row     <- grep("Precision for obs_id", rownames(hyper_summary))
    precision_mean    <- hyper_summary[precision_row, "mean"]
    result$sigma2_est <- 1 / precision_mean
  } else if (model_type == "gaussian" && !is.null(fit$summary.hyperpar)) {
    hyper_summary     <- fit$summary.hyperpar
    precision_row     <- grep("Precision for the Gaussian observations",
                              rownames(hyper_summary))
    precision_mean    <- hyper_summary[precision_row, "mean"]
    result$sigma2_est <- 1 / precision_mean
  }
  
  return(result)
}

################################################################################
