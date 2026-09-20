################################################################################
##                                                                            ##
##                          SDM.3 FITTED WITH spaMM                           ##
##                            (SDMs/models/SDM3.R)                            ##
##                                                                            ##
##       Addressing Missing Covariates in Species Distribution Models:        ##
##  Inferential Impacts and Mitigation via Joint Species Distribution Models  ##
##                                                                            ##
##                  Arthur F. Rossignol & Frédéric Gosselin                   ##
##                                                                            ##
##                                    2026                                    ##
##                                                                            ##
################################################################################

fit_SDM.3 <- function(Y, 
                      X1, 
                      model_type,
                      use_OLRE) {
  
  obs_id  <- factor(seq_along(Y))
  data_df <- data.frame(Y = Y, X1 = X1, obs_id = obs_id)
  
  get_family_spaMM <- function(model_type) {
    switch(model_type,
           "gaussian"          = gaussian(),
           "poisson"           = poisson(),
           "bernoulli_probit"  = binomial(link = "probit"),
           "bernoulli_logit"   = binomial(link = "logit"),
           "bernoulli_cloglog" = binomial(link = "cloglog"))
  }
  
  family <- get_family_spaMM(model_type)
  
  if (use_OLRE) {
    fit <- fitme(Y ~ X1 + (1 | obs_id), family = family, data = data_df)
    
    fixed_coefs <- fixef(fit)
    fit_summary <- summary(fit, verbose = FALSE)
    se_values   <- fit_summary$beta_table[, "Cond. SE"]
    
    result <- list(alpha_est = fixed_coefs[1],
                   alpha_SE  = se_values[1],
                   beta_est  = fixed_coefs[2],
                   beta_SE   = se_values[2],
                   converged = NA)
    
    lambda <- VarCorr(fit)
    if (!is.null(lambda) && length(lambda) > 0) {
      if ("obs_id" %in% names(lambda)) {
        result$sigma2_est <- as.numeric(lambda[["obs_id"]])
      } else if (length(lambda) > 0) {
        result$sigma2_est <- as.numeric(lambda[[1]])
      }
    }
    
  } else {
    fit <- fitme(Y ~ X1, family = family, data = data_df)
    
    fixed_coefs <- fixef(fit)
    fit_summary <- summary(fit, verbose = FALSE)
    se_values   <- fit_summary$beta_table[, "Cond. SE"]
    
    result <- list(alpha_est = fixed_coefs[1],
                   alpha_SE  = se_values[1],
                   beta_est  = fixed_coefs[2],
                   beta_SE   = se_values[2],
                   converged = NA)
    
    if (model_type == "gaussian") {
      result$sigma2_est <- as.numeric(fit$phi)[1]
    }
  }
  
  return(result)
}

################################################################################
