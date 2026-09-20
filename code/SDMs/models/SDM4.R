################################################################################
##                                                                            ##
##                       SDM.4 FITTED WITH GLMMadaptive                       ##
##                            (SDMs/models/SDM4.R)                            ##
##                                                                            ##
##       Addressing Missing Covariates in Species Distribution Models:        ##
##  Inferential Impacts and Mitigation via Joint Species Distribution Models  ##
##                                                                            ##
##                  Arthur F. Rossignol & Frédéric Gosselin                   ##
##                                                                            ##
##                                    2026                                    ##
##                                                                            ##
################################################################################

fit_SDM.4 <- function(Y, 
                      X1, 
                      model_type) {
  
  obs_id  <- factor(seq_along(Y))
  data_df <- data.frame(Y = Y, X1 = X1, obs_id = obs_id)
  
  get_family_GLMMadaptive <- function(model_type) {
    switch(model_type,
           "gaussian"         = gaussian(),
           "poisson"          = poisson(),
           "bernoulli_probit" = binomial(link = "probit"),
           "bernoulli_logit"  = binomial(link = "logit"))
  }
  
  family <- get_family_GLMMadaptive(model_type)
  
  fit <- mixed_model(fixed  = Y ~ X1,
                     random = ~ 1 | obs_id,
                     data   = data_df,
                     family = family)
  
  fixed_effects <- fixef(fit)
  se_fixed      <- sqrt(diag(vcov(fit)))
  
  result <- list(alpha_est = fixed_effects[1],
                 alpha_SE  = se_fixed[1],
                 beta_est  = fixed_effects[2],
                 beta_SE   = se_fixed[2],
                 converged = fit$converged)
  
  result$sigma2_est <- fit$D[1, 1]
  
  return(result)
}

################################################################################
