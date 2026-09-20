################################################################################
##                                                                            ##
##                          SDM.10 FITTED WITH mgcv                           ##
##                           (SDMs/models/SDM10.R)                            ##
##                                                                            ##
##       Addressing Missing Covariates in Species Distribution Models:        ##
##  Inferential Impacts and Mitigation via Joint Species Distribution Models  ##
##                                                                            ##
##                  Arthur F. Rossignol & Frédéric Gosselin                   ##
##                                                                            ##
##                                    2026                                    ##
##                                                                            ##
################################################################################

fit_SDM.10 <- function(Y, 
                       X1,
                       model_type,
                       use_OLRE) {
  
  obs_id  <- factor(seq_along(Y))
  data_df <- data.frame(Y = Y, X1 = X1, obs_id = obs_id)
  
  get_family_mgcv <- function(model_type) {
    switch(model_type,
           "gaussian"         = gaussian(),
           "poisson"          = poisson(),
           "bernoulli_probit" = binomial(link = "probit"),
           "bernoulli_logit"  = binomial(link = "logit"))
  }
  
  family         <- get_family_mgcv(model_type)
  smooth_formula <- Y ~ s(X1, bs = "tp", k = 10)
  
  compute_AME <- function(fit_gam, newdata) {
    d_plus     <- newdata
    d_plus$X1  <- d_plus$X1 + 1e-5
    d_minus    <- newdata
    d_minus$X1 <- d_minus$X1 - 1e-5
    Xp_plus    <- predict(fit_gam, newdata = d_plus,  type = "lpmatrix")
    Xp_minus   <- predict(fit_gam, newdata = d_minus, type = "lpmatrix")
    Xd         <- (Xp_plus - Xp_minus) / (2 * 1e-5)
    grad       <- colMeans(Xd)
    AME_est    <- sum(grad * coef(fit_gam))
    AME_SE     <- sqrt(as.numeric(t(grad) %*% vcov(fit_gam) %*% grad))
    list(est = AME_est, SE = AME_SE)
  }
  
  if (use_OLRE) {
    fit <- gamm(formula = smooth_formula,
                random  = list(obs_id = ~ 1),
                family  = family,
                data    = data_df)
    
    gam_summary <- summary(fit$gam)
    AME         <- compute_AME(fit$gam, data_df)
    
    result <- list(alpha_est = unname(coef(fit$gam)[1]),
                   alpha_SE  = gam_summary$se[1],
                   beta_est  = AME$est,
                   beta_SE   = AME$SE,
                   edf_X1    = unname(gam_summary$edf[1]),
                   converged = NA)
    
    vc                <- VarCorr(fit$lme)
    result$sigma2_est <- as.numeric(vc["obs_id", "Variance"])
    
  } else {
    fit <- gam(formula = smooth_formula,
               family  = family,
               data    = data_df)
    
    gam_summary <- summary(fit)
    AME         <- compute_AME(fit, data_df)
    
    result <- list(alpha_est = unname(coef(fit)[1]),
                   alpha_SE  = gam_summary$se[1],
                   beta_est  = AME$est,
                   beta_SE   = AME$SE,
                   edf_X1    = unname(gam_summary$edf[1]),
                   converged = fit$converged)
    
    if (model_type == "gaussian") {
      result$sigma2_est <- fit$sig2
    }
  }
  
  return(result)
}

################################################################################
