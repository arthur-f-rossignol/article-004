################################################################################
##                                                                            ##
##                        SDM.9 FITTED WITH robustbase                        ##
##                            (SDMs/models/SDM9.R)                            ##
##                                                                            ##
##       Addressing Missing Covariates in Species Distribution Models:        ##
##  Inferential Impacts and Mitigation via Joint Species Distribution Models  ##
##                                                                            ##
##                  Arthur F. Rossignol & Frédéric Gosselin                   ##
##                                                                            ##
##                                    2026                                    ##
##                                                                            ##
################################################################################

fit_SDM.9 <- function(Y, 
                      X1,
                      model_type) {
  
  data_df <- data.frame(Y = Y, X1 = X1)

  get_family_robustbase <- function(model_type) {
    switch(model_type,
           "gaussian"          = gaussian(),
           "poisson"           = poisson(link = "log"),
           "bernoulli_probit"  = binomial(link = "probit"),
           "bernoulli_logit"   = binomial(link = "logit"),
           "bernoulli_cloglog" = binomial(link = "cloglog"))
  }

  if (model_type == "gaussian") {
    fit <- lmrob(Y ~ X1, data = data_df)
  } else {
    family <- get_family_robustbase(model_type)
    fit <- glmrob(Y ~ X1, family = family, data = data_df, method = "Mqle")
  }
  
  coef_summary <- summary(fit)$coefficients
  
  result <- list(alpha_est = coef_summary["(Intercept)", "Estimate"],
                 alpha_SE  = coef_summary["(Intercept)", "Std. Error"],
                 beta_est  = coef_summary["X1", "Estimate"],
                 beta_SE   = coef_summary["X1", "Std. Error"],
                 converged = fit$converged)
  
  if (model_type == "gaussian") {
    result$sigma2_est <- fit$scale^2
  }
  
  return(result)
}

################################################################################
