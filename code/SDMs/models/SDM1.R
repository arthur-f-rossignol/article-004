################################################################################
##                                                                            ##
##                      SDM.1 FITTED WITH stats AND lme4                      ##
##                            (SDMs/models/SDM1.R)                            ##
##                                                                            ##
##       Addressing Missing Covariates in Species Distribution Models:        ##
##  Inferential Impacts and Mitigation via Joint Species Distribution Models  ##
##                                                                            ##
##                  Arthur F. Rossignol & Frédéric Gosselin                   ##
##                                                                            ##
##                                    2026                                    ##
##                                                                            ##
################################################################################

fit_SDM.1 <- function(Y,
                      X1,
                      model_type, 
                      use_OLRE) {
  
  obs_id  <- factor(seq_along(Y))
  data_df <- data.frame(Y = Y, X1 = X1, obs_id = obs_id)
  
  get_family_stats <- function(model_type) {
    switch(model_type,
           "gaussian"          = gaussian(),
           "poisson"           = poisson(link = "log"),
           "bernoulli_probit"  = binomial(link = "probit"),
           "bernoulli_logit"   = binomial(link = "logit"),
           "bernoulli_cloglog" = binomial(link = "cloglog"))
  }

  family <- get_family_stats(model_type)
  
  if (use_OLRE) {
    if (model_type == "gaussian") {
      fit <- lmer(Y ~ X1 + (1 | obs_id), data = data_df)
    } else {
      fit <- glmer(Y ~ X1 + (1 | obs_id), family = family, data = data_df)
    }
  } else {
    if (model_type == "gaussian") {
      fit <- lm(Y ~ X1, data = data_df)
    } else {
      fit <- glm(Y ~ X1, family = family, data = data_df)
    }
  }
  
  summary_fit <- summary(fit)
  
  if (use_OLRE) {
    result <- list(alpha_est  = summary_fit$coefficients[1, 1],
                   alpha_SE   = summary_fit$coefficients[1, 2],
                   beta_est   = summary_fit$coefficients[2, 1],
                   beta_SE    = summary_fit$coefficients[2, 2],
                   sigma2_est = summary_fit$varcor$obs_id[1],
                   converged  = (fit@optinfo$conv$opt == 0))
  } else {
    coef_summary <- summary(fit)$coefficients
    result <- list(alpha_est = coef_summary["(Intercept)", "Estimate"],
                   alpha_SE  = coef_summary["(Intercept)", "Std. Error"],
                   beta_est  = coef_summary["X1", "Estimate"],
                   beta_SE   = coef_summary["X1", "Std. Error"],
                   converged = if (model_type == "gaussian") TRUE else fit$converged)
    if (model_type == "gaussian") {
      result$sigma2_est <- summary_fit$sigma^2
    }
  }
  
  return(result)
}

################################################################################
