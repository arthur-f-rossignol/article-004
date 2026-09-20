################################################################################
##                                                                            ##
##                          SDM.7 FITTED WITH nimble                          ##
##                            (SDMs/models/SDM7.R)                            ##
##                                                                            ##
##       Addressing Missing Covariates in Species Distribution Models:        ##
##  Inferential Impacts and Mitigation via Joint Species Distribution Models  ##
##                                                                            ##
##                  Arthur F. Rossignol & Frédéric Gosselin                   ##
##                                                                            ##
##                                    2026                                    ##
##                                                                            ##
################################################################################

fit_SDM.7 <- function(Y, 
                      X1, 
                      model_type,
                      use_OLRE) {
  
  n_obs  <- length(Y)
  obs_id <- seq_len(n_obs)
  
  if (model_type == "gaussian") {
    nimble_code <- nimbleCode({
      alpha ~ dnorm(0, sd = 10)
      beta ~ dnorm(0, sd = 10)
      tau ~ dgamma(0.001, 0.001)
      sigma <- 1 / sqrt(tau)
      sigma2 <- sigma^2
      for(i in 1:n_obs) {
        mu[i] <- alpha + beta * X1[i]
        Y[i] ~ dnorm(mu[i], sd = sigma)
      }
    })
  } else if (model_type %in% c("bernoulli_probit")) {
    if (use_OLRE) {
      nimble_code <- nimbleCode({
        alpha ~ dnorm(0, sd = 10)
        beta ~ dnorm(0, sd = 10)
        tau ~ dgamma(0.001, 0.001)
        sigma <- 1 / sqrt(tau)
        sigma2 <- sigma^2
        for(j in 1:n_obs) {
          site_effect[j] ~ dnorm(0, sd = sigma)
        }
        for(i in 1:n_obs) {
          probit(p[i]) <- alpha + beta * X1[i] + site_effect[obs_id[i]]
          Y[i] ~ dbern(p[i])
        }
      })
    } else {
      nimble_code <- nimbleCode({
        alpha ~ dnorm(0, sd = 10)
        beta ~ dnorm(0, sd = 10)
        for(i in 1:n_obs) {
          probit(p[i]) <- alpha + beta * X1[i]
          Y[i] ~ dbern(p[i])
        }
      })
    }
  } else if (model_type %in% c("bernoulli_logit")) {
    if (use_OLRE) {
      nimble_code <- nimbleCode({
        alpha ~ dnorm(0, sd = 10)
        beta ~ dnorm(0, sd = 10)
        tau ~ dgamma(0.001, 0.001)
        sigma <- 1 / sqrt(tau)
        sigma2 <- sigma^2
        for(j in 1:n_obs) {
          site_effect[j] ~ dnorm(0, sd = sigma)
        }
        for(i in 1:n_obs) {
          logit(p[i]) <- alpha + beta * X1[i] + site_effect[obs_id[i]]
          Y[i] ~ dbern(p[i])
        }
      })
    } else {
      nimble_code <- nimbleCode({
        alpha ~ dnorm(0, sd = 10)
        beta ~ dnorm(0, sd = 10)
        for(i in 1:n_obs) {
          logit(p[i]) <- alpha + beta * X1[i]
          Y[i] ~ dbern(p[i])
        }
      })
    }
  } else if (model_type %in% c("bernoulli_cloglog")) {
    if (use_OLRE) {
      nimble_code <- nimbleCode({
        alpha ~ dnorm(0, sd = 10)
        beta ~ dnorm(0, sd = 10)
        tau ~ dgamma(0.001, 0.001)
        sigma <- 1 / sqrt(tau)
        sigma2 <- sigma^2
        for(j in 1:n_obs) {
          site_effect[j] ~ dnorm(0, sd = sigma)
        }
        for(i in 1:n_obs) {
          p[i] <- 1 - exp(-exp(alpha + beta * X1[i] + site_effect[obs_id[i]]))
          Y[i] ~ dbern(p[i])
        }
      })
    } else {
      nimble_code <- nimbleCode({
        alpha ~ dnorm(0, sd = 10)
        beta ~ dnorm(0, sd = 10)
        for(i in 1:n_obs) {
          p[i] <- 1 - exp(-exp(alpha + beta * X1[i]))
          Y[i] ~ dbern(p[i])
        }
      })
    }
  } else if (model_type %in% c("poisson")) {
    if (use_OLRE) {
      nimble_code <- nimbleCode({
        alpha ~ dnorm(0, sd = 10)
        beta ~ dnorm(0, sd = 10)
        tau ~ dgamma(0.001, 0.001)
        sigma <- 1 / sqrt(tau)
        sigma2 <- sigma^2
        for(j in 1:n_obs) {
          site_effect[j] ~ dnorm(0, sd = sigma)
        }
        for(i in 1:n_obs) {
          log(lambda[i]) <- alpha + beta * X1[i] + site_effect[obs_id[i]]
          Y[i] ~ dpois(lambda[i])
        }
      })
    } else {
      nimble_code <- nimbleCode({
        alpha ~ dnorm(0, sd = 10)
        beta ~ dnorm(0, sd = 10)
        for(i in 1:n_obs) {
          log(lambda[i]) <- alpha + beta * X1[i]
          Y[i] ~ dpois(lambda[i])
        }
      })
    }
  }
  
  nimble_data <- list(Y      = Y,
                      X1     = X1,
                      obs_id = obs_id)
  
  nimble_constants <- list(n_obs = n_obs)
  
  nimble_inits <- function(chain_id) {
    inits <- list(alpha = rnorm(1, 0, 1),
                  beta  = rnorm(1, 0, 1))
    if (model_type == "gaussian") {
      inits$tau <- 1 / runif(1, 0.1, 2)^2
    }
    if (use_OLRE && model_type != "gaussian") {
      inits$tau         <- 1 / runif(1, 0.1, 2)^2
      inits$site_effect <- rnorm(n_obs, 0, 0.1)
    }
    return(inits)
  }
  
  nimble_params <- c("alpha", "beta")
  
  if (use_OLRE || model_type == "gaussian") {
    nimble_params <- c(nimble_params, "sigma")
  }
  
  Nchains <- 3
  
  nimble_out <- runMCMC_btadjust(code         = nimble_code,
                                 constants    = nimble_constants,
                                 data         = nimble_data,
                                 inits        = lapply(1:Nchains, nimble_inits),
                                 params       = nimble_params,
                                 niter.min    = 10000,
                                 niter.max    = Inf,
                                 nburnin.min  = 10000,
                                 nburnin.max  = Inf,
                                 thin.min     = 1,
                                 thin.max     = Inf,
                                 Nchains      = Nchains,
                                 conv.max     = 1.05,
                                 neff.min     = 5000,
                                 control      = list(time.max                   = 24 * 3600,
                                                     round.thinmult             = TRUE,
                                                     print.diagnostics          = TRUE,
                                                     Ncycles.target             = 2,
                                                     check.convergence.firstrun = TRUE,
                                                     convtype                   = 'Gelman'),
                                 control.MCMC = list(parallelize = TRUE))
  
  samples <- do.call(rbind, nimble_out)
  attrs   <- attributes(nimble_out)
  
  result <- list(alpha_est = mean(samples[, "alpha"]),
                 alpha_SE  = sd(samples[, "alpha"]),
                 beta_est  = mean(samples[, "beta"]),
                 beta_SE   = sd(samples[, "beta"]),
                 converged = attrs$final.params$converged)
  
  if (use_OLRE || model_type == "gaussian") {
    result$sigma2_est <- mean(samples[, "sigma"]^2)
  }
  
  return(result)
}

################################################################################
