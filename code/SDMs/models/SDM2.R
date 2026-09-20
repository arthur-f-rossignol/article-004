################################################################################
##                                                                            ##
##                         SDM.2 FITTED WITH glmmTMB                          ##
##                            (SDMs/models/SDM2.R)                            ##
##                                                                            ##
##       Addressing Missing Covariates in Species Distribution Models:        ##
##  Inferential Impacts and Mitigation via Joint Species Distribution Models  ##
##                                                                            ##
##                  Arthur F. Rossignol & Frédéric Gosselin                   ##
##                                                                            ##
##                                    2026                                    ##
##                                                                            ##
################################################################################

fit_SDM.2 <- function(Y, 
                      X1,
                      model_type, 
                      use_OLRE) {
  
  obs_id  <- factor(seq_along(Y))
  data_df <- data.frame(Y = Y, X1 = X1, obs_id = obs_id)
  
  get_family_glmmTMB <- function(model_type) {
    switch(model_type,
           "gaussian"          = gaussian(),
           "poisson"           = poisson(link = "log"),
           "bernoulli_probit"  = binomial(link = "probit"),
           "bernoulli_logit"   = binomial(link = "logit"),
           "bernoulli_cloglog" = binomial(link = "cloglog"))
  }

  family <- get_family_glmmTMB(model_type)
  
  if (use_OLRE) {
    fit <- glmmTMB(Y ~ X1 + (1 | obs_id), family = family, data = data_df)
  } else {
    fit <- glmmTMB(Y ~ X1, family = family, data = data_df)
  }
  
  fit_summary <- summary(fit)
  result <- list(alpha_est = fit_summary$coefficients$cond[1, 1],
                 alpha_SE  = fit_summary$coefficients$cond[1, 2],
                 beta_est  = fit_summary$coefficients$cond[2, 1],
                 beta_SE   = fit_summary$coefficients$cond[2, 2])
  
  if (use_OLRE) {
    result$sigma2_est <- fit_summary$varcor$cond$obs_id[1]
  } else if (model_type == "gaussian") {
    result$sigma2_est <- sigma(fit)^2
  }
  
  result$converged <- fit$sdr$pdHess
  
  return(result)
}

################################################################################
