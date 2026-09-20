################################################################################
##                                                                            ##
##                         SDM.8 FITTED WITH runjags                          ##
##                            (SDMs/models/SDM8.R)                            ##
##                                                                            ##
##       Addressing Missing Covariates in Species Distribution Models:        ##
##  Inferential Impacts and Mitigation via Joint Species Distribution Models  ##
##                                                                            ##
##                  Arthur F. Rossignol & Frédéric Gosselin                   ##
##                                                                            ##
##                                    2026                                    ##
##                                                                            ##
################################################################################

fit_SDM.8 <- function(Y, 
                      X1, 
                      model_type, 
                      use_OLRE) {
  
  n_obs  <- length(Y)
  obs_id <- seq_len(n_obs)
  
  jags_models <- list()
  
  jags_models$gaussian_no_OLRE <- "
    model {
      alpha ~ dnorm(0, 0.01)
      beta ~ dnorm(0, 0.01)
      tau ~ dgamma(0.001, 0.001)
      sigma <- 1 / sqrt(tau)
      sigma2 <- sigma^2
      for (i in 1:n_obs) {
        mu[i] <- alpha + beta * X1[i]
        Y[i] ~ dnorm(mu[i], tau)
      }
    }"
  
  jags_models$bernoulli_probit_no_OLRE <- "
    model {
      alpha ~ dnorm(0, 0.01)
      beta ~ dnorm(0, 0.01)
      for (i in 1:n_obs) {
        z[i] <- alpha + beta * X1[i]
        p[i] <- phi(z[i])
        Y[i] ~ dbern(p[i])
      }
    }"
  
  jags_models$bernoulli_probit_OLRE <- "
    model {
      alpha ~ dnorm(0, 0.01)
      beta ~ dnorm(0, 0.01)
      tau ~ dgamma(0.001, 0.001)
      sigma <- 1 / sqrt(tau)
      sigma2 <- sigma^2
      for (j in 1:n_obs) {
        obs_effect[j] ~ dnorm(0, tau)
      }
      for (i in 1:n_obs) {
        z[i] <- alpha + beta * X1[i] + obs_effect[obs_id[i]]
        p[i] <- phi(z[i])
        Y[i] ~ dbern(p[i])
      }
    }"
  
  jags_models$bernoulli_logit_no_OLRE <- "
    model {
      alpha ~ dnorm(0, 0.01)
      beta ~ dnorm(0, 0.01)
      for (i in 1:n_obs) {
        logit(p[i]) <- alpha + beta * X1[i]
        Y[i] ~ dbern(p[i])
      }
    }"
  
  jags_models$bernoulli_logit_OLRE <- "
    model {
      alpha ~ dnorm(0, 0.01)
      beta ~ dnorm(0, 0.01)
      tau ~ dgamma(0.001, 0.001)
      sigma <- 1 / sqrt(tau)
      sigma2 <- sigma^2
      for (j in 1:n_obs) {
        obs_effect[j] ~ dnorm(0, tau)
      }
      for (i in 1:n_obs) {
        logit(p[i]) <- alpha + beta * X1[i] + obs_effect[obs_id[i]]
        Y[i] ~ dbern(p[i])
      }
    }"
  
  jags_models$bernoulli_cloglog_no_OLRE <- "
    model {
      alpha ~ dnorm(0, 0.01)
      beta ~ dnorm(0, 0.01)
      for (i in 1:n_obs) {
        p[i] <- 1 - exp(-exp(alpha + beta * X1[i]))
        Y[i] ~ dbern(p[i])
      }
    }"
  
  jags_models$bernoulli_cloglog_OLRE <- "
    model {
      alpha ~ dnorm(0, 0.01)
      beta ~ dnorm(0, 0.01)
      tau ~ dgamma(0.001, 0.001)
      sigma <- 1 / sqrt(tau)
      sigma2 <- sigma^2
      for (j in 1:n_obs) {
        obs_effect[j] ~ dnorm(0, tau)
      }
      for (i in 1:n_obs) {
        p[i] <- 1 - exp(-exp(alpha + beta * X1[i] + obs_effect[obs_id[i]]))
        Y[i] ~ dbern(p[i])
      }
    }"
  
  jags_models$poisson_no_OLRE <- "
    model {
      alpha ~ dnorm(0, 0.01)
      beta  ~ dnorm(0, 0.01)
      for (i in 1:n_obs) {
        log(lambda[i]) <- alpha + beta * X1[i]
        Y[i] ~ dpois(lambda[i])
      }
    }"
  
  jags_models$poisson_OLRE <- "
    model {
      alpha ~ dnorm(0, 0.01)
      beta  ~ dnorm(0, 0.01)
      tau ~ dgamma(0.001, 0.001)
      sigma <- 1 / sqrt(tau)
      sigma2 <- sigma^2
      for (j in 1:n_obs) {
        obs_effect[j] ~ dnorm(0, tau)
      }
      for (i in 1:n_obs) {
        log(lambda[i]) <- alpha + beta * X1[i] + obs_effect[obs_id[i]]
        Y[i] ~ dpois(lambda[i])
      }
    }"
  
  model_key         <- paste0(model_type, ifelse(use_OLRE, "_OLRE", "_no_OLRE"))
  jags_model_string <- jags_models[[model_key]]
  
  jags_model <- tempfile(fileext = ".txt")
  writeLines(jags_model_string, jags_model)
  
  jags_data <- list(Y      = as.numeric(Y),
                    X1     = as.numeric(X1),
                    obs_id = as.integer(obs_id),
                    n_obs  = as.integer(n_obs))
  
  jags_params <- c("alpha", "beta")
  
  if (use_OLRE || model_type == "gaussian") {
    jags_params <- c(jags_params, "sigma2")
  }
  
  make_inits <- function(chain_id) {
    inits <- list(alpha = rnorm(1, 0, 1),
                  beta  = rnorm(1, 0, 1))
    if (model_type == "gaussian") {
      sd0       <- runif(1, 0.5, 2)
      inits$tau <- 1 / (sd0^2)
    }
    if (use_OLRE) {
      sd_OLRE           <- runif(1, 0.1, 1.0)
      inits$tau         <- 1 / (sd_OLRE^2)
      inits$obs_effect  <- rnorm(n_obs, 0, sd = sd_OLRE * 0.5)
    }
    return(inits)
  }
  
  Nchains <- 3
  
  jags_out <- runMCMC_btadjust(MCMC_language = "Jags",
                               code          = jags_model,
                               data          = jags_data,
                               inits         = lapply(1:Nchains, make_inits),
                               params        = jags_params,
                               niter.min     = 10000,
                               niter.max     = Inf,
                               nburnin.min   = 10000,
                               nburnin.max   = Inf,
                               thin.min      = 1,
                               thin.max      = Inf,
                               Nchains       = Nchains,
                               conv.max      = 1.05,
                               neff.min      = 5000,
                               control       = list(time.max                   = 24 * 3600,
                                                    round.thinmult             = TRUE,
                                                    print.diagnostics          = TRUE,
                                                    Ncycles.target             = 2,
                                                    check.convergence.firstrun = TRUE,
                                                    convtype                   = 'Gelman'),
                               control.MCMC  = list(parallelize = TRUE))
  
  combined_samples <- do.call(rbind, jags_out)
  attrs            <- attributes(jags_out)
  
  result <- list(alpha_est = mean(combined_samples[, "alpha"]),
                 alpha_SE  = sd(combined_samples[, "alpha"]),
                 beta_est  = mean(combined_samples[, "beta"]),
                 beta_SE   = sd(combined_samples[, "beta"]),
                 converged = attrs$final.params$converged)
  
  if ("sigma2" %in% colnames(combined_samples)) {
    result$sigma2_est = mean(combined_samples[, "sigma2"])
  }
  
  unlink(jags_model)
  
  return(result)
}

################################################################################
