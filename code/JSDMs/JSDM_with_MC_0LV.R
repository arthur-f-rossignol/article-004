################################################################################
##                                                                            ##
##          JSDMs WITH ONE MISSING COVARIATE AND NO LATENT VARIABLE           ##
##                                                                            ##
##       Addressing Missing Covariates in Species Distribution Models:        ##
##  Inferential Impacts and Mitigation via Joint Species Distribution Models  ##
##                                                                            ##
##                  Arthur F. Rossignol & Frédéric Gosselin                   ##
##                                                                            ##
##                                    2026                                    ##
##                                                                            ##
################################################################################

## PACKAGES ####################################################################

library(gllvm)
library(jSDM)
library(Hmsc)
library(nimble)
library(runjags)
library(coda)
library(parallel)
library(runMCMCbtadjust)

## MAIN PARAMETERS #############################################################

DEFAULT_PARAMS <- list(n      = 500,      # number of sites
                       p      = 10,       # number of species
                       q      = 2,        # number of covariates (excluding intercept)
                       num_lv = 0,        # number of latent variables 
                       mcmc   = list(niter_min      = 25000,
                                     nburnin_min    = 25000,
                                     thin_min       = 10,
                                     nchains        = 3,
                                     conv_max       = 1.05,
                                     neff_min       = 1000,
                                     time_max_hours = 12))

## SCENARIO PARAMETERIZATION ###################################################

build_scenarios <- function(p) {
  list(S.1 = list(name       = "S.1",
                  true_alpha = rep(1, p),
                  true_beta  = rep(1, p),
                  true_gamma = rep(1, p)),

       S.2 = list(name       = "S.2",
                  true_alpha = rep(0, p),
                  true_beta  = rep(0.5, p),
                  true_gamma = rep(0.5, p)),

       S.3 = list(name       = "S.3",
                  true_alpha = rep(-1.5, p),
                  true_beta  = rep(0.5, p),
                  true_gamma = rep(0.5, p)),

       S.4 = list(name       = "S.4",
                  true_alpha = seq(-1.5, 1.5, length.out = p),
                  true_beta  = rep(0.5, p),
                  true_gamma = rep(0.5, p)),

       S.5 = list(name       = "S.5",
                  true_alpha = rep(0, p),
                  true_beta  = rep(1.5, p),
                  true_gamma = rep(0.2, p)),

       S.6 = list(name       = "S.6",
                  true_alpha = seq(-1, 1, length.out = p),
                  true_beta  = rep(1, p),
                  true_gamma = rep(-0.5, p)))
}

## DATA GENERATING PROCESS #####################################################

simulate_data <- function(seed, 
                          scenario, 
                          params = NULL) {
  
  set.seed(seed)
  
  params <- DEFAULT_PARAMS
  
  n <- params$n
  p <- params$p

  X1 <- rnorm(n, mean = 0, sd = 1)
  X2 <- rnorm(n, mean = 0, sd = 1)
  X  <- cbind(intercept = 1, X1 = X1, X2 = X2)

  expand_coef <- function(coef, p) {
    if (length(coef) == 1) {
      rep(coef, p) else coef
    }
  }

  alphas <- expand_coef(scenario$true_alpha, p)
  betas  <- expand_coef(scenario$true_beta, p)
  gammas <- expand_coef(scenario$true_gamma, p)

  B <- rbind(alphas, betas, gammas)
  rownames(B) <- c("alpha", "beta", "gamma")
  colnames(B) <- paste0("sp", 1:p)

  probabilities <- pnorm(as.vector(X %*% B))
  
  Y <- matrix(rbinom(n * p, size = 1, prob = probabilities),
              nrow = n,
              ncol = p,
              dimnames = list(NULL, paste0("sp", 1:p)))

  prevalence <- colMeans(Y)

  list(Y        = Y,
       X1       = X1,
       X2       = X2,
       X        = X,
       B        = B,
       scenario = scenario,
       params   = params)
}

## RESULT TEMPLATE INITIALIZATION ##############################################

initialize_results <- function(p, 
                               num_lv) {
  
  list(alpha_estimates       = rep(NA_real_, p),
       alpha_standard_errors = rep(NA_real_, p),
       beta_estimates        = rep(NA_real_, p),
       beta_standard_errors  = rep(NA_real_, p),
       num_latent_variables  = num_lv,
       latent_variables      = NULL,
       loadings              = NULL,
       convergence           = FALSE,
       ess_alpha             = NA_real_,
       ess_beta              = NA_real_,
       psrf_alpha            = NA_real_,
       psrf_beta             = NA_real_,
       psrf_max              = NA_real_,
       computation_time      = NA_real_)
}

## JSDM.1 ######################################################################

fit_JSDM.1 <- function(Y, 
                       X1,  
                       num_lv,
                       seed) {
  
  set.seed(seed)
  
  p <- ncol(Y)
  
  result <- initialize_results(p, num_lv)
  
  start_time <- Sys.time()

  fit <- gllvm::gllvm(y             = Y,
                      X             = data.frame(X1 = X1),
                      family        = binomial(link = "probit"),
                      method        = "LA",
                      num.lv        = num_lv,
                      formula       = ~ 1 + X1,
                      seed          = seed,
                      trace         = FALSE,
                      control       = list(optimizer = "nlminb", 
                                           max.iter  = 1e6),
                      control.start = list(starting.val = "res", 
                                           n.init       = 50, 
                                           n.init.max   = 1e4))
  
  result$computation_time <- as.numeric(difftime(Sys.time(), 
                                                 start_time, 
                                                 units = "secs"))
  
  result$alpha_estimates       <- as.vector(fit$params$beta0)
  result$beta_estimates        <- as.vector(fit$params$Xcoef)
  result$alpha_standard_errors <- as.vector(fit$sd$beta0)
  result$beta_standard_errors  <- as.vector(fit$sd$Xcoef)

  result$convergence <- fit$convergence

  result
}

## JSDM.2 ######################################################################

fit_JSDM.2 <- function(Y, 
                       X1, 
                       num_lv, 
                       seed) {
  
  set.seed(seed)
  
  p <- ncol(Y)
  n <- nrow(Y)
  
  result <- initialize_results(p, num_lv)
 
  start_time <- Sys.time()
  
  fit <- jSDM::jSDM_binomial_probit(burnin        = 100000,
                                    mcmc          = 100000,
                                    thin          = 5,
                                    presence_data = Y,
                                    site_formula  = ~ X1,
                                    site_data     = data.frame(X1 = X1),
                                    n_latent      = num_lv,
                                    site_effect   = "none",
                                    beta_start    = 0,
                                    lambda_start  = 0,
                                    W_start       = 0,
                                    alpha_start   = 0,
                                    V_alpha       = 1,
                                    mu_beta       = 0,
                                    V_beta        = 1,
                                    mu_lambda     = 0,
                                    V_lambda      = 1,
                                    shape_Valpha  = 0.5,
                                    rate_Valpha   = 0.0005,
                                    seed          = seed,
                                    verbose       = 1)
  
  result$computation_time <- as.numeric(difftime(Sys.time(), 
                                                 start_time, 
                                                 units = "secs"))
  
  n_covars <- 2  
  
  mean_beta <- matrix(NA, nrow = p, ncol = n_covars)
  se_beta   <- matrix(NA, nrow = p, ncol = n_covars)

  for (j in 1:p) {
    sp_samples <- as.matrix(fit$mcmc.sp[[j]])

    for (cv in 1:n_covars) {
      mean_beta[j, cv] <- mean(sp_samples[, cv])
      se_beta[j, cv]   <- sd(sp_samples[, cv])
    }
  }

  result$alpha_estimates       <- mean_beta[, 1]
  result$alpha_standard_errors <- se_beta[, 1]
  result$beta_estimates        <- mean_beta[, 2]
  result$beta_standard_errors  <- se_beta[, 2]

  ess_alpha <- numeric(p)
  ess_beta  <- numeric(p)
  for (j in 1:p) {
    sp_mcmc <- fit$mcmc.sp[[j]]
    ess_all <- coda::effectiveSize(sp_mcmc)
    if (length(ess_all) >= 1) {
      ess_alpha[j] <- ess_all[1]
    }
    if (length(ess_all) >= 2) {
      ess_beta[j]  <- ess_all[2]
    }
  }
  result$ess_alpha <- mean(ess_alpha, na.rm = TRUE)
  result$ess_beta  <- mean(ess_beta, na.rm = TRUE)

  geweke_vals <- numeric(p)
  for (j in 1:p) {
    gw <- coda::geweke.diag(fit$mcmc.sp[[j]])$z
    if (is.numeric(gw) && length(gw) >= 2) {
      geweke_vals[j] <- max(abs(gw[1:2]), na.rm = TRUE)
    }
  }
  max_geweke <- max(geweke_vals, na.rm = TRUE)
  result$convergence <- is.finite(max_geweke) && max_geweke < 2.0
  result$psrf_max    <- NA
  result$geweke_max  <- max_geweke

  if (!is.null(fit$mcmc.Deviance)) {
    deviance_samples     <- as.matrix(fit$mcmc.Deviance)
    result$deviance_mean <- mean(deviance_samples)
    result$deviance_sd   <- sd(deviance_samples)
  }

  result
}

## JSDM.3 ######################################################################

fit_JSDM.3 <- function(Y,
                       X1, 
                       num_lv,
                       seed) {
  
  set.seed(seed)
  
  p <- ncol(Y)
  n <- nrow(Y)
  
  result <- initialize_results(p, num_lv = num_lv)
  
  model <- Hmsc::Hmsc(Y        = Y,
                      XData    = data.frame(X1 = X1),
                      XFormula = ~ 1 + X1,
                      distr    = "probit")
  
  start_time <- Sys.time()
  
  fit <- Hmsc::sampleMcmc(model,
                          thin      = 5,
                          samples   = 50000,
                          transient = 50000,
                          nChains   = 3,
                          nParallel = 3,
                          verbose   = 0,
                          initPar   = "fixed effects")
  
  result$computation_time <- as.numeric(difftime(Sys.time(), 
                                                 start_time, 
                                                 units = "secs"))
  
  postBeta <- Hmsc::getPostEstimate(fit, parName = "Beta")
  
  result$alpha_estimates <- postBeta$mean[1, ]
  result$beta_estimates  <- postBeta$mean[2, ]

  mpost <- Hmsc::convertToCodaObject(fit)
  
  beta_samples  <- do.call(rbind, mpost$Beta)
  intercept_idx <- seq(1, p * 2, by = 2)
  x1_idx        <- seq(2, p * 2, by = 2)
  
  result$alpha_standard_errors <- apply(beta_samples[, intercept_idx, drop = FALSE], 2, sd)
  result$beta_standard_errors  <- apply(beta_samples[, x1_idx, drop = FALSE], 2, sd)

  ess_values       <- effectiveSize(mpost$Beta)
  result$ess_alpha <- mean(ess_values[intercept_idx], na.rm = TRUE)
  result$ess_beta  <- mean(ess_values[x1_idx], na.rm = TRUE)

  if (length(mpost$Beta) > 1) {
    psrf_values        <- gelman.diag(mpost$Beta, multivariate = FALSE)$psrf[, "Point est."]
    result$psrf_max    <- max(psrf_values, na.rm = TRUE)
    result$convergence <- is.finite(result$psrf_max) && result$psrf_max < 1.05
  }

  result
}

## JSDM.4 ######################################################################

fit_JSDM.4 <- function(Y,
                       X1, 
                       num_lv,
                       seed) {

  set.seed(seed)
  
  Y     <- as.matrix(Y)
  X_mat <- cbind(1, as.matrix(X1))

  n_obs     <- nrow(Y)
  p_species <- ncol(Y)
  n_covars  <- ncol(X_mat)

  result <- initialize_results(p_species, num_lv)

  modelCode <- nimble::nimbleCode({

    for (i in 1:nvar) {
      meanenvsppar[i] ~ dnorm(0, 1)
    }
    V[1:nvar, 1:nvar] ~ dinvwish(S = V0[1:nvar, 1:nvar], df = f0)

    for (j in 1:nspecies) {
      envsppar[1:nvar, j] ~ dmnorm(meanenvsppar[1:nvar], cov = V[1:nvar, 1:nvar])
    }

    if (nlvbigger0) {
      sigmalambdapar0 ~ dgamma(a[1], b[1])
      sigmalambdapar[1] <- sigmalambdapar0

      if (nlvbigger1) {
        for (i in 1:(nlv - 1)) {
          probslambdapar[i] ~ dgamma(a[2], b[2])
          sigmalambdapar[i + 1] <- sigmalambdapar[i] * probslambdapar[i]
        }
      }

      for (i in 1:nlv) {
        meanlambdapar[i] <- 0

        phi[i, i] ~ dgamma(nu / 2, nu / 2)
        lambdapar_raw_d[i] ~ dnorm(0, sigmalambdapar[i] * phi[i, i])
        lambdapar[i, i] <- meanlambdapar[i] + lambdapar_raw_d[i]
        lambdapar_abs_d[i] <- abs(lambdapar_raw_d[i])

        for (j in 1:nsites) {
          upar[j, i] ~ dnorm(0, sd = 1)
          upar_abs[j, i] <- abs(upar[j, i])
        }
      }

      for (i in 1:sum_idx_s) {
        phib[i] ~ dgamma(nu / 2, nu / 2)
        lambdapar_raw_s_np[i] ~ dnorm(0, sigmalambdapar[indices_idx_s[i, 2]] * phib[i])
      }

      for (i in 1:n_idx_s) {
        lambdapar[indices_idx_s[i, 1], indices_idx_s[i, 2]] <-
          meanlambdapar[indices_idx_s[i, 2]] + lambdapar_raw_s_np[i]
        lambdapar_abs_s_np[i] <- abs(lambdapar[indices_idx_s[i, 1], indices_idx_s[i, 2]])
      }

      if (runn_idx_u) {
        for (i in 1:n_idx_u) {
          phibb[i] ~ dgamma(nu / 2, nu / 2)
          lambdapar_raw_u_np[i] ~ dnorm(0, sigmalambdapar[indices_idx_u[i, 2]] * phibb[i])
          lambdapar[indices_idx_u[i, 1], indices_idx_u[i, 2]] <- lambdapar_raw_u_np[i]
          lambdapar_abs_u_np[i] <- abs(lambdapar_raw_u_np[i])
        }
      }
    }

    for (i in 1:nobs) {
      if (nvarbigger1) {
        CLel.nlvar[i] <- sum(envsppar[1:nvar, beta0[i]] * env[i, 1:nvar])
      } else {
        CLel.nlvar[i] <- envsppar[1, beta0[i]] * env[i, 1]
      }

      if (nlvbigger1) {
        CLel.nlv[i] <- sum(lambdapar[beta0lambda[i], 1:nlv] * upar[alpha[i], 1:nlv])
      } else {
        if (nlvbigger0) {
          CLel.nlv[i] <- lambdapar[beta0lambda[i], 1] * upar[alpha[i], 1]
        } else {
          CLel.nlv[i] <- 0.0
        }
      }

      CL[i]  <- CLel.nlvar[i]
      CL2[i] <- CLel.nlv[i]
      CL3[i] ~ dnorm(CL[i] + CL2[i], 1)
      CL4[i] <- step(CL3[i])
      Y[i]   ~ dbern(CL4[i])
    }
  })

  Y_long      <- as.vector(t(Y))
  env_long    <- X_mat[rep(1:n_obs, each = p_species), , drop = FALSE]
  beta0       <- rep(1:p_species, times = n_obs)
  beta0lambda <- beta0
  alpha       <- rep(1:n_obs, each = p_species)
  nobs        <- length(Y_long)

  loading_info <- list(indices_idx_s = matrix(0, nrow = 0, ncol = 2),
                       n_idx_s       = 0,
                       sum_idx_s     = 0,
                       indices_idx_u = matrix(0, nrow = 0, ncol = 2),
                       n_idx_u       = 0,
                       runn_idx_u    = FALSE)

  modelConsts <- list(nobs        = nobs,
                      nsites      = n_obs,
                      nspecies    = p_species,
                      nvar        = n_covars,
                      nlv         = num_lv,
                      nlvbigger0  = (num_lv > 0),
                      nlvbigger1  = (num_lv > 1),
                      nvarbigger1 = (n_covars > 1),
                      beta0       = beta0,
                      beta0lambda = beta0lambda,
                      alpha       = alpha,
                      a           = c(50, 50),
                      b           = c(1, 1),
                      nu          = 3,
                      V0          = diag(n_covars),
                      f0          = n_covars + 1)

  modelData <- list(Y   = Y_long, 
                    env = env_long)

  make_inits <- function() {
    inits <- list(meanenvsppar    = rnorm(n_covars, 0, 0.1),
                  V               = diag(n_covars) * 0.5,
                  envsppar        = matrix(rnorm(n_covars * p_species, 0, 0.2),
                                           n_covars, p_species),
                  CL3             = rnorm(nobs, 0, 0.5))

    inits
  }

  inits_list <- replicate(3, make_inits(), simplify = FALSE)

  params <- c("meanenvsppar", "envsppar", "V")

  mcmc_params <- DEFAULT_PARAMS$mcmc
  
  start_time <- Sys.time()
  
  out <- runMCMCbtadjust::runMCMC_btadjust(
    MCMC_language = "Nimble",
    code          = modelCode,
    constants     = modelConsts,
    data          = modelData,
    inits         = inits_list,
    params        = params,
    niter.min     = mcmc_params$niter_min,
    niter.max     = Inf,
    nburnin.min   = mcmc_params$nburnin_min,
    nburnin.max   = Inf,
    thin.min      = mcmc_params$thin_min,
    thin.max      = Inf,
    Nchains       = length(inits_list),
    conv.max      = mcmc_params$conv_max,
    neff.min      = mcmc_params$neff_min,
    control = list(time.max                   = 3600 * mcmc_params$time_max_hours,
                   round.thinmult             = TRUE,
                   print.diagnostics          = TRUE,
                   Ncycles.target             = 2,
                   check.convergence.firstrun = TRUE,
                   convtype                   = "Gelman",
                   seed                       = seed),
    control.MCMC = list(parallelize = TRUE))
  
  result$computation_time <- as.numeric(difftime(Sys.time(), 
                                                 start_time, 
                                                 units = "secs"))
  
  combined <- do.call(rbind, out)
  attrs    <- attributes(out)

  envsp_cols <- grep("^envsppar\\[", colnames(combined))
  if (length(envsp_cols) > 0) {
    envsp_arr <- array(NA, dim = c(nrow(combined), n_covars, p_species))
    for (idx in envsp_cols) {
      nm <- colnames(combined)[idx]
      ij <- as.numeric(strsplit(gsub("envsppar\\[|\\]", "", nm), ",")[[1]])
      envsp_arr[, ij[1], ij[2]] <- combined[, idx]
    }

    for (sp in 1:p_species) {
      result$alpha_estimates[sp]       <- mean(envsp_arr[, 1, sp], na.rm = TRUE)
      result$alpha_standard_errors[sp] <- sd(envsp_arr[, 1, sp], na.rm = TRUE)
      if (n_covars >= 2) {
        result$beta_estimates[sp]       <- mean(envsp_arr[, 2, sp], na.rm = TRUE)
        result$beta_standard_errors[sp] <- sd(envsp_arr[, 2, sp], na.rm = TRUE)
      }
    }

    result$ess_alpha <- mean(sapply(1:p_species, function(sp) {
      effectiveSize(as.mcmc(envsp_arr[, 1, sp]))
    }), na.rm = TRUE)

    if (n_covars >= 2) {
      result$ess_beta <- mean(sapply(1:p_species, function(sp) {
        effectiveSize(as.mcmc(envsp_arr[, 2, sp]))
      }), na.rm = TRUE)
    }
  }

  result$convergence <- attrs$final.params$converged

  result
}

## JSDM.5 ######################################################################

fit_JSDM.5 <- function(Y, 
                       X1, 
                       num_lv, 
                       seed) {
  
  set.seed(seed)
  
  Y     <- as.matrix(Y)
  X_mat <- as.matrix(X1)
  if (is.vector(X_mat)) {
    X_mat <- matrix(X_mat, ncol = 1)
  }

  n_obs     <- nrow(Y)
  p_species <- ncol(Y)
  n_covars  <- ncol(X_mat)

  result <- initialize_results(p_species, num_lv)
  
  model_string <- "
    model {
      for(i in 1:n) {
        for(j in 1:p) {
          probit(p_y[i, j]) <- beta0[j] + inprod(beta[j, ], X[i, ])
          y[i, j] ~ dbern(p_y[i, j])
        }
      }

      for(j in 1:p) {
        beta0[j] ~ dnorm(0, 0.1)
      }

      for(j in 1:p) {
        for(m in 1:n_covars) {
          beta[j,m] ~ dnorm(0, 0.1)
        }
      }
    }"

  jags_data <- list(y        = Y, 
                    X        = X_mat,
                    n        = n_obs, 
                    p        = p_species,
                    n_covars = n_covars)

  make_inits <- function() {
    inits <- list(beta0 = rnorm(p_species, 0, 0.5),
                  beta  = matrix(rnorm(p_species * n_covars, 0, 0.3), 
                                 p_species, n_covars))
    inits
  }
  
  inits_list <- replicate(3, make_inits(), simplify = FALSE)

  monitor_params <- c("beta0", "beta")

  mcmc_params <- DEFAULT_PARAMS$mcmc

  random_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  model_file <- file.path(
    tempdir(),
    sprintf("jsdm5_%s_%04d.txt", random_id, sample.int(9999, 1)))
  writeLines(model_string, model_file)

  start_time <- Sys.time()
  
  out <- runMCMCbtadjust::runMCMC_btadjust(
    MCMC_language = "Jags",
    code          = model_file,
    data          = jags_data,
    inits         = inits_list,
    params        = monitor_params,
    params.conv   = c("beta0", "beta"),
    niter.min     = mcmc_params$niter_min,
    niter.max     = Inf,
    nburnin.min   = mcmc_params$nburnin_min,
    nburnin.max   = Inf,
    thin.min      = mcmc_params$thin_min,
    thin.max      = Inf,
    Nchains       = length(inits_list),
    conv.max      = mcmc_params$conv_max,
    neff.min      = mcmc_params$neff_min,
    control = list(time.max                   = 3600 * mcmc_params$time_max_hours,
                   round.thinmult             = TRUE,
                   print.diagnostics          = FALSE,
                   Ncycles.target             = 2,
                   check.convergence.firstrun = FALSE,
                   convtype                   = "Gelman",
                   seed                       = seed),
    control.MCMC = list(parallelize = TRUE))
  
  result$computation_time <- as.numeric(difftime(Sys.time(), 
                                                 start_time, 
                                                 units = "secs"))
  
  combined <- do.call(rbind, out)
  attrs    <- attributes(out)

  beta0_cols <- grep("^beta0\\[", colnames(combined))
  beta0_mat  <- matrix(NA, nrow = nrow(combined), ncol = p_species)
  for (idx in beta0_cols) {
    j <- as.numeric(gsub("beta0\\[|\\]", "", colnames(combined)[idx]))
    beta0_mat[, j] <- combined[, idx]
  }

  beta_cols <- grep("^beta\\[", colnames(combined))
  beta_arr  <- array(NA, dim = c(nrow(combined), p_species, n_covars))
  for (idx in beta_cols) {
    nm <- colnames(combined)[idx]
    ij <- as.numeric(strsplit(gsub("beta\\[|\\]", "", nm), ",")[[1]])
    beta_arr[, ij[1], ij[2]] <- combined[, idx]
  }

  if (num_lv > 0) {
    lambda_cols <- grep("^lambda\\[", colnames(combined))
    lambda_arr  <- array(NA, dim = c(nrow(combined), p_species, num_lv))
    for (idx in lambda_cols) {
      nm <- colnames(combined)[idx]
      ij <- as.numeric(strsplit(gsub("lambda\\[|\\]", "", nm), ",")[[1]])
      lambda_arr[, ij[1], ij[2]] <- combined[, idx]
    }
    result$loadings <- apply(lambda_arr, c(2, 3), mean, na.rm = TRUE)

    W_cols <- grep("^W\\[", colnames(combined))
    if (length(W_cols) > 0) {
      W_arr <- array(NA, dim = c(nrow(combined), n_obs, num_lv))
      for (idx in W_cols) {
        nm <- colnames(combined)[idx]
        ik <- as.numeric(strsplit(gsub("W\\[|\\]", "", nm), ",")[[1]])
        W_arr[, ik[1], ik[2]] <- combined[, idx]
      }
      result$latent_variables <- apply(W_arr, c(2, 3), mean)
    }
  }

  for (sp in 1:p_species) {
    result$alpha_estimates[sp]       <- mean(beta0_mat[, sp], na.rm = TRUE)
    result$alpha_standard_errors[sp] <- sd(beta0_mat[, sp], na.rm = TRUE)
    result$beta_estimates[sp]        <- mean(beta_arr[, sp, 1], na.rm = TRUE)
    result$beta_standard_errors[sp]  <- sd(beta_arr[, sp, 1], na.rm = TRUE)
  }

  ess_alpha        <- sapply(1:p_species, function(sp) effectiveSize(as.mcmc(beta0_mat[, sp])))
  ess_beta         <- sapply(1:p_species, function(sp) effectiveSize(as.mcmc(beta_arr[, sp, 1])))
  result$ess_alpha <- mean(ess_alpha, na.rm = TRUE)
  result$ess_beta  <- mean(ess_beta, na.rm = TRUE)

  result$convergence <- attrs$final.params$converged

  unlink(model_file)

  result
}

## JSDM.6 ######################################################################

fit_JSDM.6 <- function(Y, 
                       X1, 
                       num_lv, 
                       seed) {

  set.seed(seed)
  
  Y     <- as.matrix(Y)
  X_mat <- as.matrix(X1)
  if (is.vector(X_mat)) {
    X_mat <- matrix(X_mat, ncol = 1)
  }

  n_obs     <- nrow(Y)
  p_species <- ncol(Y)
  n_covars  <- ncol(X_mat)

  result <- initialize_results(p_species, num_lv)

  model_string <- "
    model {
      for(i in 1:n) {
        for(j in 1:p) {
          probit(p_y[i, j]) <- beta0[j] + inprod(beta[j, ], X[i, ])
          y[i, j] ~ dbern(p_y[i, j])
        }
      }
      for(j in 1:p) {
        beta0[j] ~ dnorm(0, 0.1)
        for(m in 1:n_covars) {
          beta[j,m] ~ dnorm(0, 0.1)
        }
      }
    }"

  jags_data <- list(y        = Y, 
                    X        = X_mat,
                    n        = n_obs, 
                    p        = p_species, 
                    n_covars = n_covars)

  make_inits <- function() {
    inits <- list(beta0 = rnorm(p_species, 0, 0.5),
                  beta  = matrix(rnorm(p_species * n_covars, 0, 0.3), 
                                 p_species, n_covars))
    inits
  }
  inits_list <- replicate(3, make_inits(), simplify = FALSE)

  monitor_params <- c("beta0", "beta")
  
  mcmc_params <- DEFAULT_PARAMS$mcmc

  random_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  model_file <- file.path(tempdir(),
                          sprintf("jsdm6_%s_%04d.txt", random_id, sample.int(9999, 1)))
  writeLines(model_string, model_file)
  
  start_time <- Sys.time()
  
  out <- runMCMCbtadjust::runMCMC_btadjust(
    MCMC_language = "Jags",
    code          = model_file,
    data          = jags_data,
    inits         = inits_list,
    params        = monitor_params,
    params.conv   = c("beta0", "beta"),
    niter.min     = mcmc_params$niter_min,
    niter.max     = Inf,
    nburnin.min   = mcmc_params$nburnin_min,
    nburnin.max   = Inf,
    thin.min      = mcmc_params$thin_min,
    thin.max      = Inf,
    Nchains       = length(inits_list),
    conv.max      = mcmc_params$conv_max,
    neff.min      = mcmc_params$neff_min,
    control       = list(time.max                   = 3600 * mcmc_params$time_max_hours,
                         round.thinmult             = TRUE,
                         print.diagnostics          = FALSE,
                         Ncycles.target             = 2,
                         check.convergence.firstrun = FALSE,
                         convtype                   = "Gelman",
                         seed                       = seed),
    control.MCMC  = list(parallelize = TRUE))
  
  result$computation_time <- as.numeric(difftime(Sys.time(), 
                                                 start_time, 
                                                 units = "secs"))
  
  combined <- do.call(rbind, out)
  attrs    <- attributes(out)

  beta0_cols <- grep("^beta0\\[", colnames(combined))
  beta0_mat  <- matrix(NA, nrow = nrow(combined), ncol = p_species)
  for (idx in beta0_cols) {
    j <- as.numeric(gsub("beta0\\[|\\]", "", colnames(combined)[idx]))
    beta0_mat[, j] <- combined[, idx]
  }

  beta_cols <- grep("^beta\\[", colnames(combined))
  beta_arr  <- array(NA, dim = c(nrow(combined), p_species, n_covars))
  for (idx in beta_cols) {
    nm <- colnames(combined)[idx]
    ij <- as.numeric(strsplit(gsub("beta\\[|\\]", "", nm), ",")[[1]])
    beta_arr[, ij[1], ij[2]] <- combined[, idx]
  }

  for (sp in 1:p_species) {
    result$alpha_estimates[sp]       <- mean(beta0_mat[, sp], na.rm = TRUE)
    result$alpha_standard_errors[sp] <- sd(beta0_mat[, sp], na.rm = TRUE)
    result$beta_estimates[sp]        <- mean(beta_arr[, sp, 1], na.rm = TRUE)
    result$beta_standard_errors[sp]  <- sd(beta_arr[, sp, 1], na.rm = TRUE)
  }

  ess_alpha        <- sapply(1:p_species, function(sp) effectiveSize(as.mcmc(beta0_mat[, sp])))
  ess_beta         <- sapply(1:p_species, function(sp) effectiveSize(as.mcmc(beta_arr[, sp, 1])))
  result$ess_alpha <- mean(ess_alpha, na.rm = TRUE)
  result$ess_beta  <- mean(ess_beta, na.rm = TRUE)

  result$convergence <- attrs$final.params$converged

  unlink(model_file)

  result
}

## DATA ANALYSIS ###############################################################

run_JSDM <- function(seed, 
                     output_dir, 
                     scenario_id,
                     methods = 1:6, 
                     params = NULL) {

  if (is.null(params)) {
    params <- DEFAULT_PARAMS
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  scenarios <- build_scenarios(params$p)
  scenario <- scenarios[[scenario_id]]

  data <- simulate_data(seed     = seed, 
                        scenario = scenario, 
                        params   = params)

  results <- list(seed        = seed,
                  timestamp   = Sys.time(),
                  scenario_id = scenario_id,
                  scenario    = scenario,
                  parameters  = list(n      = params$n,
                                     p      = params$p,
                                     q      = params$q,
                                     num_lv = params$num_lv,
                                     B      = data$B))

  total_start <- Sys.time()

  if (1 %in% methods) {
    results$JSDM.1 <- fit_JSDM.1(Y      = data$Y,
                                 X1     = data$X1,
                                 num_lv = params$num_lv,
                                 seed   = seed)
  }

  if (2 %in% methods) {
    results$JSDM.2 <- fit_JSDM.2(Y      = data$Y,
                                 X1     = data$X1,
                                 num_lv = params$num_lv,
                                 seed   = seed)
  }

  if (3 %in% methods) {
    results$JSDM.3 <- fit_JSDM.3(Y      = data$Y, 
                                 X1     = data$X1,
                                 num_lv = params$num_lv,
                                 seed   = seed)
  }

  if (4 %in% methods) {
    results$JSDM.4 <- fit_JSDM.4(Y      = data$Y,
                                 X1.    = data$X1,
                                 num_lv = params$num_lv,
                                 seed   = seed)
  }

  if (5 %in% methods) {
    results$JSDM.5 <- fit_JSDM.5(Y      = data$Y,
                                 X1     = data$X1,
                                 num_lv = params$num_lv,
                                 seed   = seed)
  }

  if (6 %in% methods) {
    results$JSDM.6 <- fit_JSDM.6(Y      = data$Y,
                                 X1     = data$X1,
                                 num_lv = params$num_lv,
                                 seed   = seed)
  }

  results$total_computation_time <- as.numeric(difftime(Sys.time(), 
                                                        total_start,
                                                        units = "secs"))

  if (length(methods) == 1) {
    output_file <- file.path(output_dir,
                             sprintf("JSDM%d_S%d_seed%06d.RData", methods,
                                     scenario_id, seed))
  } else {
    method_str <- paste(methods, collapse = "_")
    output_file <- file.path(output_dir, 
                             sprintf("JSDM_S%d_seed%06d_methods%s.RData", 
                                     scenario_id, seed, method_str))
  }
  
  save(results, file = output_file)

  for (m in methods) {
    method_name <- paste0("JSDM.", m)
    if (!is.null(results[[method_name]])) {
      r <- results[[method_name]]
    }
  }
}

## RUN #########################################################################

args <- commandArgs(trailingOnly = TRUE)

seed_value  <- as.integer(args[1])
output_dir  <- args[2]
scenario_id <- if (length(args) >= 3) as.integer(args[3]) else 1L
methods     <- if (length(args) >= 4) as.integer(strsplit(args[4], ",")[[1]]) else 1:6

run_JSDM(seed        = seed_value,
         output_dir  = output_dir,
         scenario_id = scenario_id,
         methods     = methods)

################################################################################
