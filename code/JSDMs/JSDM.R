################################################################################
##                                                                            ##
##          JSDMs WITH OPTIONAL MISSING COVARIATE AND LATENT VARIABLES        ##
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

## ARGUMENTS FROM SLURM ########################################################

ARGS <- commandArgs(trailingOnly = TRUE)

seed <- as.integer(ARGS[1])
results_dir <- ARGS[2]

if (length(ARGS) >= 3) {
  scenario_id <- as.integer(ARGS[3]) 
} else {
  scenario_id <- 1
}

if (length(ARGS) >= 4) {
  methods <- as.integer(strsplit(ARGS[4], ",")[[1]]) 
} else {
  methods <- 1:6
}

if (length(ARGS) >= 5) {
  n_lv <- as.integer(ARGS[5])
}

if (length(ARGS) >= 6) {
  missing_covariate <- as.logical(ARGS[6])
}

## MAIN PARAMETERS #############################################################

DEFAULT_PARAMS <- list(n_obs             = 1000,     # number of sites
                       n_sp              = 10,       # number of species
                       n_covars          = 2,        # number of covariates (intercept excluded)
                       n_lv              = 1,        # number of latent variables
                       missing_covariate = TRUE,     # X2 is missing (unobserved)
                       mcmc              = list(n_iter_min     = 25000,
                                                n_burnin_min   = 25000,
                                                thin_min       = 10,
                                                n_chains       = 3,
                                                conv_max       = 1.05,
                                                n_eff_min      = 1000,
                                                time_max_hours = 24))

## SCENARIO PARAMETERIZATION ###################################################

build_scenarios <- function(n_sp) {
  list(S.1 = list(name       = "S.1",
                  true_alpha = rep(1, n_sp),
                  true_beta  = rep(1, n_sp),
                  true_gamma = rep(1, n_sp)),

       S.2 = list(name       = "S.2",
                  true_alpha = rep(0, n_sp),
                  true_beta  = rep(0.5, n_sp),
                  true_gamma = rep(0.5, n_sp)),

       S.3 = list(name       = "S.3",
                  true_alpha = rep(-1.5, n_sp),
                  true_beta  = rep(0.5, n_sp),
                  true_gamma = rep(0.5, n_sp)),

       S.4 = list(name       = "S.4",
                  true_alpha = seq(-1.5, 1.5, length.out = n_sp),
                  true_beta  = rep(0.5, n_sp),
                  true_gamma = rep(0.5, n_sp)),

       S.5 = list(name       = "S.5",
                  true_alpha = rep(0, n_sp),
                  true_beta  = rep(1.5, n_sp),
                  true_gamma = rep(0.2, n_sp)),

       S.6 = list(name       = "S.6",
                  true_alpha = seq(-1, 1, length.out = n_sp),
                  true_beta  = rep(1, n_sp),
                  true_gamma = rep(-0.5, n_sp)))
}

## DATA GENERATING PROCESS #####################################################

simulate_data <- function(seed,
                          scenario,
                          params = NULL) {

  set.seed(seed)

  params <- DEFAULT_PARAMS

  n_obs <- params$n_obs
  n_sp  <- params$n_sp

  X1 <- rnorm(n_obs, mean = 0, sd = 1)
  X2 <- rnorm(n_obs, mean = 0, sd = 1)
  X  <- cbind(intercept = 1, X1 = X1, X2 = X2)

  expand_coef <- function(coef, n_sp) {
    if (length(coef) == 1) {
      rep(coef, n_sp)
    } else {
      coef
    }
  }

  alphas <- expand_coef(scenario$true_alpha, n_sp)
  betas  <- expand_coef(scenario$true_beta, n_sp)
  gammas <- expand_coef(scenario$true_gamma, n_sp)

  B <- rbind(alphas, betas, gammas)
  rownames(B) <- c("alpha", "beta", "gamma")
  colnames(B) <- paste0("sp", 1:n_sp)

  probabilities <- pnorm(as.vector(X %*% B))

  Y <- matrix(rbinom(n_obs * n_sp, size = 1, prob = probabilities),
              nrow = n_obs,
              ncol = n_sp,
              dimnames = list(NULL, paste0("sp", 1:n_sp)))

  prevalence <- colMeans(Y)

  list(Y        = Y,
       X1       = X1,
       X2       = X2,
       X        = X,
       B        = B,
       scenario = scenario,
       params   = params)
}

## COVARIATE SELECTION #########################################################

select_covariates <- function(data,
                              missing_covariate) {

  if (missing_covariate) {
    list(X1 = data$X1)
  } else {
    list(X1 = data$X1,
         X2 = data$X2)
  }
}

## RESULT TEMPLATE INITIALIZATION ##############################################

results_initialization <- function(n_sp,
                                   n_lv,
                                   fit_gamma) {

  result <- list(alpha_estimates         = rep(NA, n_sp),
                 alpha_standard_errors   = rep(NA, n_sp),
                 beta_estimates          = rep(NA, n_sp),
                 beta_standard_errors    = rep(NA, n_sp),
                 latent_variables        = NULL,
                 loadings                = NULL,
                 convergence             = NA,
                 ess_alpha               = NA,
                 ess_beta                = NA,
                 psrf_alpha              = NA,
                 psrf_beta               = NA,
                 psrf_max                = NA,
                 computation_time        = NA)

  if (fit_gamma) {
    result$gamma_estimates       <- rep(NA, n_sp)
    result$gamma_standard_errors <- rep(NA, n_sp)
    result$ess_gamma             <- NA
  }

  return(result)
}

## JSDM.1 ######################################################################

fit_JSDM.1 <- function(Y,
                       covariates,
                       n_lv,
                       seed) {

  set.seed(seed)

  n_sp <- ncol(Y)

  fit_gamma <- (length(covariates) >= 2)

  result <- results_initialization(n_sp, n_lv, fit_gamma = fit_gamma)

  X_data  <- as.data.frame(covariates)
  formula <- as.formula(paste("~ 1 +", paste(names(covariates), collapse = " + ")))

  start_time <- Sys.time()

  fit <- gllvm::gllvm(y             = Y,
                      X             = X_data,
                      family        = binomial(link = "probit"),
                      method        = "LA",
                      num.lv        = n_lv,
                      formula       = formula,
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

  Xcoef_estimates       <- fit$params$Xcoef
  Xcoef_standard_errors <- fit$sd$Xcoef

  result$alpha_estimates       <- as.vector(fit$params$beta0)
  result$alpha_standard_errors <- as.vector(fit$sd$beta0)
  result$beta_estimates        <- as.vector(Xcoef_estimates[, 1])
  result$beta_standard_errors  <- as.vector(Xcoef_standard_errors[, 1])

  if (fit_gamma) {
    result$gamma_estimates       <- as.vector(Xcoef_estimates[, 2])
    result$gamma_standard_errors <- as.vector(Xcoef_standard_errors[, 2])
  }

  if (fit$num.lv > 0) {
    result$latent_variables <- fit$lvs
    result$loadings         <- fit$params$theta
    result$sigma_lv         <- fit$params$sigma.lv
  }
  
  result$convergence <- fit$convergence

  return(result)
}

## JSDM.2 ######################################################################

fit_JSDM.2 <- function(Y,
                       covariates,
                       n_lv,
                       seed) {

  set.seed(seed)

  n_sp  <- ncol(Y)
  n_obs <- nrow(Y)

  fit_gamma <- (length(covariates) >= 2)

  result <- results_initialization(n_sp, n_lv, fit_gamma = fit_gamma)

  X_data       <- as.data.frame(covariates)
  site_formula <- as.formula(paste("~", paste(names(covariates), collapse = " + ")))

  start_time <- Sys.time()

  fit <- jSDM::jSDM_binomial_probit(burnin        = 100000,
                                    mcmc          = 100000,
                                    thin          = 5,
                                    presence_data = Y,
                                    site_formula  = site_formula,
                                    site_data     = X_data,
                                    n_latent      = n_lv,
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

  n_covars <- length(covariates) + 1

  mean_beta <- matrix(NA, nrow = n_sp, ncol = n_covars)
  se_beta   <- matrix(NA, nrow = n_sp, ncol = n_covars)

  if (n_lv > 0) {
    mean_lambda <- matrix(NA, nrow = n_sp, ncol = n_lv)
  }

  for (j in 1:n_sp) {
    sp_samples <- as.matrix(fit$mcmc.sp[[j]])

    for (covar_id in 1:n_covars) {
      mean_beta[j, covar_id] <- mean(sp_samples[, covar_id])
      se_beta[j, covar_id]   <- sd(sp_samples[, covar_id])
    }

    if (n_lv > 0) {
      for (l in 1:n_lv) {
        mean_lambda[j, l] <- mean(sp_samples[, n_covars + l])
      }
    }
  }

  result$alpha_estimates       <- mean_beta[, 1]
  result$alpha_standard_errors <- se_beta[, 1]
  result$beta_estimates        <- mean_beta[, 2]
  result$beta_standard_errors  <- se_beta[, 2]

  if (fit_gamma) {
    result$gamma_estimates       <- mean_beta[, 3]
    result$gamma_standard_errors <- se_beta[, 3]
  }

  if (n_lv > 0) {
    result$loadings <- mean_lambda
  }

  if (n_lv > 0 && !is.null(fit$mcmc.latent) && length(fit$mcmc.latent) > 0) {
    W_mean <- matrix(NA, nrow = n_obs, ncol = n_lv)
    for (l in 1:n_lv) {
      lv_name <- paste0("lv_", l)
      if (lv_name %in% names(fit$mcmc.latent)) {
        lv_samples <- as.matrix(fit$mcmc.latent[[lv_name]])
        W_mean[, l] <- colMeans(lv_samples)
      }
    }
    result$latent_variables     <- W_mean
  }

  ess_alpha <- numeric(n_sp)
  ess_beta  <- numeric(n_sp)
  if (fit_gamma) {
    ess_gamma <- numeric(n_sp)
  }
  for (j in 1:n_sp) {
    sp_mcmc <- fit$mcmc.sp[[j]]
    ess_all <- coda::effectiveSize(sp_mcmc)
    if (length(ess_all) >= 1) {
      ess_alpha[j] <- ess_all[1]
    }
    if (length(ess_all) >= 2) {
      ess_beta[j]  <- ess_all[2]
    }
    if (fit_gamma && length(ess_all) >= 3) {
      ess_gamma[j] <- ess_all[3]
    }
  }
  result$ess_alpha <- mean(ess_alpha, na.rm = TRUE)
  result$ess_beta  <- mean(ess_beta, na.rm = TRUE)
  if (fit_gamma) {
    result$ess_gamma <- mean(ess_gamma, na.rm = TRUE)
  }

  geweke_vals <- numeric(n_sp)
  for (j in 1:n_sp) {
    geweke_z <- coda::geweke.diag(fit$mcmc.sp[[j]])$z
    if (is.numeric(geweke_z) && length(geweke_z) >= n_covars) {
      geweke_vals[j] <- max(abs(geweke_z[1:n_covars]), na.rm = TRUE)
    }
  }
  max_geweke <- max(geweke_vals, na.rm = TRUE)
  result$convergence   <- is.finite(max_geweke) && max_geweke < 2.0
  result$psrf_beta_max <- NA
  result$geweke_max    <- max_geweke

  if (!is.null(fit$mcmc.Deviance)) {
    deviance_samples     <- as.matrix(fit$mcmc.Deviance)
    result$deviance_mean <- mean(deviance_samples)
    result$deviance_sd   <- sd(deviance_samples)
  }

  if (!is.null(fit$probit_theta_latent)) {
    result$probit_theta <- fit$probit_theta_latent
  }
  if (!is.null(fit$theta_latent)) {
    result$theta <- fit$theta_latent
  }

  return(result)
}

## JSDM.3 ######################################################################

fit_JSDM.3 <- function(Y,
                       covariates,
                       n_lv,
                       seed) {

  set.seed(seed)

  n_sp  <- ncol(Y)
  n_obs <- nrow(Y)

  fit_gamma <- length(covariates) >= 2

  result <- results_initialization(n_sp, n_lv = n_lv, fit_gamma = fit_gamma)

  X_data   <- as.data.frame(covariates)
  XFormula <- as.formula(paste("~ 1 +", paste(names(covariates), collapse = " + ")))

  if (n_lv > 0) {
    studyDesign        <- data.frame(site = as.factor(1:n_obs))
    random_level       <- Hmsc::HmscRandomLevel(units = studyDesign$site)
    random_level$nfMin <- n_lv
    random_level$nfMax <- n_lv

    model <- Hmsc::Hmsc(Y           = Y,
                        XData       = X_data,
                        XFormula    = XFormula,
                        studyDesign = studyDesign,
                        ranLevels   = list(site = random_level),
                        distr       = "probit")
  } else {
    model <- Hmsc::Hmsc(Y        = Y,
                        XData    = X_data,
                        XFormula = XFormula,
                        distr    = "probit")
  }

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
  if (fit_gamma) {
    result$gamma_estimates <- postBeta$mean[3, ]
  }

  mpost <- Hmsc::convertToCodaObject(fit)

  beta_samples <- do.call(rbind, mpost$Beta)

  n_fixed <- length(covariates) + 1   # intercept + fitted covariates

  alpha_id <- seq(1, n_sp * n_fixed, by = n_fixed)
  beta_id  <- seq(2, n_sp * n_fixed, by = n_fixed)

  result$alpha_standard_errors <- apply(beta_samples[, alpha_id, drop = FALSE], 2, sd)
  result$beta_standard_errors  <- apply(beta_samples[, beta_id, drop = FALSE], 2, sd)

  if (fit_gamma) {
    gamma_id                    <- seq(3, n_sp * n_fixed, by = n_fixed)
    result$gamma_standard_errors <- apply(beta_samples[, gamma_id, drop = FALSE], 2, sd)
  }

  if (n_lv > 0) {
    postEta    <- Hmsc::getPostEstimate(fit, parName = "Eta")
    postLambda <- Hmsc::getPostEstimate(fit, parName = "Lambda")

    if (!is.null(postEta)) {
      result$latent_variables     <- postEta$mean
    }

    if (!is.null(postLambda)) {
      result$loadings <- postLambda$mean
    }
  }

  ess_values       <- effectiveSize(mpost$Beta)
  result$ess_alpha <- mean(ess_values[alpha_id], na.rm = TRUE)
  result$ess_beta  <- mean(ess_values[beta_id], na.rm = TRUE)
  if (fit_gamma) {
    result$ess_gamma <- mean(ess_values[gamma_id], na.rm = TRUE)
  }

  if (length(mpost$Beta) > 1) {
    psrf_values          <- gelman.diag(mpost$Beta, multivariate = FALSE)$psrf[, "Point est."]
    result$psrf_beta_max <- max(psrf_values, na.rm = TRUE)
    result$psrf_max      <- result$psrf_beta_max
    result$convergence   <- is.finite(result$psrf_beta_max) && result$psrf_beta_max < 1.05
  }

  return(result)
}

## JSDM.4 ######################################################################

fit_JSDM.4 <- function(Y,
                       covariates,
                       n_lv,
                       seed) {

  set.seed(seed)

  Y     <- as.matrix(Y)
  X_mat <- cbind(1, do.call(cbind, lapply(covariates, as.matrix)))

  n_obs    <- nrow(Y)
  n_sp     <- ncol(Y)
  n_covars <- ncol(X_mat)

  fit_gamma <- (length(covariates) >= 2)

  result <- results_initialization(n_sp, n_lv, fit_gamma = fit_gamma)

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

      for (i in 1:sum_id_s) {
        phib[i] ~ dgamma(nu / 2, nu / 2)
        lambdapar_raw_s_np[i] ~ dnorm(0, sigmalambdapar[indices_id_s[i, 2]] * phib[i])
      }

      for (i in 1:n_id_s) {
        lambdapar[indices_id_s[i, 1], indices_id_s[i, 2]] <-
          meanlambdapar[indices_id_s[i, 2]] + lambdapar_raw_s_np[i]
        lambdapar_abs_s_np[i] <- abs(lambdapar[indices_id_s[i, 1], indices_id_s[i, 2]])
      }

      if (runn_id_u) {
        for (i in 1:n_id_u) {
          phibb[i] ~ dgamma(nu / 2, nu / 2)
          lambdapar_raw_u_np[i] ~ dnorm(0, sigmalambdapar[indices_id_u[i, 2]] * phibb[i])
          lambdapar[indices_id_u[i, 1], indices_id_u[i, 2]] <- lambdapar_raw_u_np[i]
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
  env_long    <- X_mat[rep(1:n_obs, each = n_sp), , drop = FALSE]
  beta0       <- rep(1:n_sp, times = n_obs)
  beta0lambda <- beta0
  alpha       <- rep(1:n_obs, each = n_sp)
  nobs        <- length(Y_long)

  if (n_lv == 0) {
    loading_info <- list(indices_id_s = matrix(0, nrow = 0, ncol = 2),
                         n_id_s       = 0,
                         sum_id_s     = 0,
                         indices_id_u = matrix(0, nrow = 0, ncol = 2),
                         n_id_u       = 0,
                         runn_id_u    = FALSE)
  } else {
    indices_s <- matrix(0, nrow = 0, ncol = 2)
    for (j in 1:n_sp) {
      for (k in 1:n_lv) {
        if (j > k && j <= n_lv) {
          indices_s <- rbind(indices_s, c(j, k))
        }
      }
    }

    indices_u <- matrix(0, nrow = 0, ncol = 2)
    for (j in 1:n_sp) {
      for (k in 1:n_lv) {
        if (j < k) {
          indices_u <- rbind(indices_u, c(j, k))
        }
      }
    }

    if (n_sp > n_lv) {
      for (j in (n_lv + 1):n_sp) {
        for (k in 1:n_lv) {
          indices_s <- rbind(indices_s, c(j, k))
        }
      }
    }

    loading_info <- list(indices_id_s = indices_s,
                         n_id_s       = nrow(indices_s),
                         sum_id_s     = nrow(indices_s),
                         indices_id_u = indices_u,
                         n_id_u       = nrow(indices_u),
                         runn_id_u    = nrow(indices_u) > 0)
  }

  modelConsts <- list(nobs        = nobs,
                      nsites      = n_obs,
                      nspecies    = n_sp,
                      nvar        = n_covars,
                      nlv         = n_lv,
                      nlvbigger0  = (n_lv > 0),
                      nlvbigger1  = (n_lv > 1),
                      nvarbigger1 = (n_covars > 1),
                      beta0       = beta0,
                      beta0lambda = beta0lambda,
                      alpha       = alpha,
                      a           = c(50, 50),
                      b           = c(1, 1),
                      nu          = 3,
                      V0          = diag(n_covars),
                      f0          = n_covars + 1)

  if (n_lv > 0) {
    modelConsts$indices_id_s <- loading_info$indices_id_s
    modelConsts$n_id_s       <- loading_info$n_id_s
    modelConsts$sum_id_s     <- loading_info$sum_id_s
    modelConsts$indices_id_u <- loading_info$indices_id_u
    modelConsts$n_id_u       <- loading_info$n_id_u
    modelConsts$runn_id_u    <- loading_info$runn_id_u
  }

  modelData <- list(Y   = Y_long,
                    env = env_long)

  make_inits <- function() {
    inits <- list(meanenvsppar    = rnorm(n_covars, 0, 0.1),
                  V               = diag(n_covars) * 0.5,
                  envsppar        = matrix(rnorm(n_covars * n_sp, 0, 0.2),
                                           n_covars, n_sp),
                  CL3             = rnorm(nobs, 0, 0.5))

    if (n_lv > 0) {
      inits$sigmalambdapar0 <- 1.0
      inits$upar            <- matrix(rnorm(n_obs * n_lv, 0, 0.3), n_obs, n_lv)
      inits$lambdapar_raw_d <- rnorm(n_lv, 0, 0.3)
      inits$phi             <- matrix(1, n_lv, n_lv)
      diag(inits$phi)       <- rep(1.5, n_lv)

      inits$lambdapar <- matrix(0, n_sp, n_lv)
      for (i in 1:min(n_lv, n_sp)) {
        inits$lambdapar[i, i] <- inits$lambdapar_raw_d[i]
      }

      if (loading_info$sum_id_s > 0) {
        inits$lambdapar_raw_s_np <- rnorm(loading_info$sum_id_s, 0, 0.2)
        inits$phib               <- rep(1.5, loading_info$sum_id_s)
        for (i in 1:loading_info$n_id_s) {
          inits$lambdapar[loading_info$indices_id_s[i, 1],
                          loading_info$indices_id_s[i, 2]] <- inits$lambdapar_raw_s_np[i]
        }
      }

      if (n_lv > 1) {
        inits$probslambdapar <- rep(0.9, n_lv - 1)
        inits$sigmalambdapar <- rep(1, n_lv)
        inits$sigmalambdapar[1] <- inits$sigmalambdapar0
        for (i in 2:n_lv) {
          inits$sigmalambdapar[i] <- inits$sigmalambdapar[i - 1] * inits$probslambdapar[i - 1]
        }
      } else {
        inits$sigmalambdapar <- inits$sigmalambdapar0
      }

      if (loading_info$runn_id_u) {
        inits$lambdapar_raw_u_np <- rnorm(loading_info$n_id_u, 0, 0.2)
        inits$phibb              <- rep(1.5, loading_info$n_id_u)
        for (i in 1:loading_info$n_id_u) {
          inits$lambdapar[loading_info$indices_id_u[i, 1],
                          loading_info$indices_id_u[i, 2]] <- inits$lambdapar_raw_u_np[i]
        }
      }
    }

    return(inits)
  }

  mcmc_params <- DEFAULT_PARAMS$mcmc
  
  inits_list <- replicate(mcmc_params$n_chains, make_inits(), simplify = FALSE)

  params <- c("meanenvsppar", 
              "envsppar", 
              "V")
  if (n_lv > 0) {
    params <- c(params, 
                "upar", 
                "lambdapar", 
                "sigmalambdapar")
  }

  start_time <- Sys.time()

  out <- runMCMCbtadjust::runMCMC_btadjust(MCMC_language = "Nimble",
                                           code          = modelCode,
                                           constants     = modelConsts,
                                           data          = modelData,
                                           inits         = inits_list,
                                           params        = params,
                                           niter.min     = mcmc_params$n_iter_min,
                                           niter.max     = Inf,
                                           nburnin.min   = mcmc_params$n_burnin_min,
                                           nburnin.max   = Inf,
                                           thin.min      = mcmc_params$thin_min,
                                           thin.max      = Inf,
                                           Nchains       = mcmc_params$n_chains,
                                           conv.max      = mcmc_params$conv_max,
                                           neff.min      = mcmc_params$n_eff_min,
                                           control       = list(time.max                   = 3600 * mcmc_params$time_max_hours,
                                                                round.thinmult             = TRUE,
                                                                print.diagnostics          = TRUE,
                                                                Ncycles.target             = 2,
                                                                check.convergence.firstrun = TRUE,
                                                                convtype                   = "Gelman",
                                                                seed                       = seed),
                                           control.MCMC  = list(parallelize = TRUE))

  result$computation_time <- as.numeric(difftime(Sys.time(),
                                                 start_time,
                                                 units = "secs"))

  combined <- do.call(rbind, out)
  attrs    <- attributes(out)

  neff_values <- attrs$final.diags$neff
  neff_names  <- rownames(neff_values)

  envsp_cols <- grep("^envsppar\\[", colnames(combined))
  if (length(envsp_cols) > 0) {
    envsp_arr <- array(NA, dim = c(nrow(combined), n_covars, n_sp))
    for (id in envsp_cols) {
      col_name      <- colnames(combined)[id]
      param_indices <- as.numeric(strsplit(gsub("envsppar\\[|\\]", "", col_name), ",")[[1]])
      envsp_arr[, param_indices[1], param_indices[2]] <- combined[, id]
    }

    for (species in 1:n_sp) {
      result$alpha_estimates[species]       <- mean(envsp_arr[, 1, species], na.rm = TRUE)
      result$alpha_standard_errors[species] <- sd(envsp_arr[, 1, species], na.rm = TRUE)
      if (n_covars >= 2) {
        result$beta_estimates[species]       <- mean(envsp_arr[, 2, species], na.rm = TRUE)
        result$beta_standard_errors[species] <- sd(envsp_arr[, 2, species], na.rm = TRUE)
      }
      if (fit_gamma && n_covars >= 3) {
        result$gamma_estimates[species]       <- mean(envsp_arr[, 3, species], na.rm = TRUE)
        result$gamma_standard_errors[species] <- sd(envsp_arr[, 3, species], na.rm = TRUE)
      }
    }

    result$ess_alpha <- mean(neff_values[grep("^envsppar\\[1,", neff_names)], na.rm = TRUE)

    if (n_covars >= 2) {
      result$ess_beta <- mean(neff_values[grep("^envsppar\\[2,", neff_names)], na.rm = TRUE)
    }

    if (fit_gamma && n_covars >= 3) {
      result$ess_gamma <- mean(neff_values[grep("^envsppar\\[3,", neff_names)], na.rm = TRUE)
    }
  }

  if (n_lv > 0) {
    u_cols <- grep("^upar\\[", colnames(combined))
    if (length(u_cols) > 0) {
      u_arr <- array(NA, dim = c(nrow(combined), n_obs, n_lv))
      for (id in u_cols) {
        col_name      <- colnames(combined)[id]
        param_indices <- as.numeric(strsplit(gsub("upar\\[|\\]", "", col_name), ",")[[1]])
        u_arr[, param_indices[1], param_indices[2]] <- combined[, id]
      }
      result$latent_variables <- apply(u_arr, c(2, 3), mean)
    }

    lambda_cols <- grep("^lambdapar\\[", colnames(combined))
    if (length(lambda_cols) > 0) {
      lambda_arr <- array(NA, dim = c(nrow(combined), n_sp, n_lv))
      for (id in lambda_cols) {
        col_name      <- colnames(combined)[id]
        param_indices <- as.numeric(strsplit(gsub("lambdapar\\[|\\]", "", col_name), ",")[[1]])
        lambda_arr[, param_indices[1], param_indices[2]] <- combined[, id]
      }
      result$loadings <- apply(lambda_arr, c(2, 3), mean)
    }
  }

  result$convergence <- attrs$final.params$converged

  return(result)
}

## JSDM.5 ######################################################################

fit_JSDM.5 <- function(Y,
                       covariates,
                       n_lv,
                       seed) {

  set.seed(seed)

  Y     <- as.matrix(Y)
  X_mat <- do.call(cbind, lapply(covariates, as.matrix))
  if (is.vector(X_mat)) {
    X_mat <- matrix(X_mat, ncol = 1)
  }

  n_obs    <- nrow(Y)
  n_sp     <- ncol(Y)
  n_covars <- ncol(X_mat)

  fit_gamma <- (length(covariates) >= 2)

  result <- results_initialization(n_sp, n_lv, fit_gamma = fit_gamma)

  if (n_lv > 0) {
    model_str <- "
      model {
        for(i in 1:n_obs) {
          for(j in 1:n_sp) {
            eta[i, j] <- inprod(lambda[j, ], W[i, ]) + inprod(beta[j, ], X[i, ])
            probit(p_y[i, j]) <- beta0[j] + eta[i, j]
            y[i, j] ~ dbern(p_y[i, j])
          }
        }
        for(i in 1:n_obs) {
          for(k in 1:n_lv) {
            W[i,k] ~ dnorm(0, 1)
          }
        }
        for(j in 1:n_sp) {
          beta0[j] ~ dnorm(0, 0.1)
        }
        for(i in 1:(n_lv - 1)) {
          for(j in (i + 1):n_lv) {
            lambda[i, j] <- 0
          }
        }
        for(i in 1:n_lv) {
          lambda[i, i] ~ dnorm(0, 0.1) T(0,)
        }
        for(i in 2:n_lv) {
          for(j in 1:(i - 1)) {
            lambda[i, j] ~ dnorm(0, 0.1)
          }
        }
        for(i in (n_lv + 1):n_sp) {
          for(j in 1:n_lv) {
            lambda[i, j] ~ dnorm(0, 0.1)
          }
        }
        for(j in 1:n_sp) {
          for(m in 1:n_covars) {
            beta[j, m] ~ dnorm(0, 0.1)
          }
        }
      }"
  } else {
    model_str <- "
      model {
        for(i in 1:n_obs) {
          for(j in 1:n_sp) {
            probit(p_y[i, j]) <- beta0[j] + inprod(beta[j, ], X[i, ])
            y[i, j] ~ dbern(p_y[i, j])
          }
        }
        for(j in 1:n_sp) {
          beta0[j] ~ dnorm(0, 0.1)
        }
        for(j in 1:n_sp) {
          for(m in 1:n_covars) {
            beta[j, m] ~ dnorm(0, 0.1)
          }
        }
      }"
  }

  data_list <- list(y        = Y,
                    X        = X_mat,
                    n_obs    = n_obs,
                    n_sp     = n_sp,
                    n_covars = n_covars)
  if (n_lv > 0) {
    data_list$n_lv <- n_lv
  }

  make_inits <- function() {
    inits <- list(beta0 = rnorm(n_sp, 0, 0.5),
                  beta  = matrix(rnorm(n_sp * n_covars, 0, 0.3),
                                 n_sp, n_covars))
    if (n_lv > 0) {
      inits$W <- matrix(rnorm(n_obs * n_lv, 0, 0.5), n_obs, n_lv)
    }
    return(inits)
  }
  
  mcmc_params <- DEFAULT_PARAMS$mcmc
  
  inits_list <- replicate(mcmc_params$n_chains, make_inits(), simplify = FALSE)

  monitor_params <- c("beta0",
                      "beta")
  if (n_lv > 0) {
    monitor_params <- c(monitor_params,
                        "lambda", 
                        "W")
  }

  random_id <- format(Sys.time(), "%Y%m%d_%H%M%S")

  model_file <- file.path(tempdir(),
                          sprintf("jsdm5_%s_%04d.txt", random_id, sample.int(9999, 1)))
  writeLines(model_str, model_file)

  start_time <- Sys.time()

  out <- runMCMCbtadjust::runMCMC_btadjust(MCMC_language = "Jags",
                                           code          = model_file,
                                           data          = data_list,
                                           inits         = inits_list,
                                           params        = monitor_params,
                                           params.conv   = c("beta0", "beta"),
                                           niter.min     = mcmc_params$n_iter_min,
                                           niter.max     = Inf,
                                           nburnin.min   = mcmc_params$n_burnin_min,
                                           nburnin.max   = Inf,
                                           thin.min      = mcmc_params$thin_min,
                                           thin.max      = Inf,
                                           Nchains       = mcmc_params$n_chains,
                                           conv.max      = mcmc_params$conv_max,
                                           neff.min      = mcmc_params$n_eff_min,
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

  neff_values <- attrs$final.diags$neff
  neff_names  <- rownames(neff_values)

  beta0_cols <- grep("^beta0\\[", colnames(combined))
  beta0_mat  <- matrix(NA, nrow = nrow(combined), ncol = n_sp)
  for (id in beta0_cols) {
    j <- as.numeric(gsub("beta0\\[|\\]", "", colnames(combined)[id]))
    beta0_mat[, j] <- combined[, id]
  }

  beta_cols <- grep("^beta\\[", colnames(combined))
  beta_arr  <- array(NA, dim = c(nrow(combined), n_sp, n_covars))
  for (id in beta_cols) {
    col_name      <- colnames(combined)[id]
    param_indices <- as.numeric(strsplit(gsub("beta\\[|\\]", "", col_name), ",")[[1]])
    beta_arr[, param_indices[1], param_indices[2]] <- combined[, id]
  }

  if (n_lv > 0) {
    lambda_cols <- grep("^lambda\\[", colnames(combined))
    lambda_arr  <- array(NA, dim = c(nrow(combined), n_sp, n_lv))
    for (id in lambda_cols) {
      col_name      <- colnames(combined)[id]
      param_indices <- as.numeric(strsplit(gsub("lambda\\[|\\]", "", col_name), ",")[[1]])
      lambda_arr[, param_indices[1], param_indices[2]] <- combined[, id]
    }
    result$loadings <- apply(lambda_arr, c(2, 3), mean, na.rm = TRUE)

    W_cols <- grep("^W\\[", colnames(combined))
    if (length(W_cols) > 0) {
      W_arr <- array(NA, dim = c(nrow(combined), n_obs, n_lv))
      for (id in W_cols) {
        col_name      <- colnames(combined)[id]
        param_indices <- as.numeric(strsplit(gsub("W\\[|\\]", "", col_name), ",")[[1]])
        W_arr[, param_indices[1], param_indices[2]] <- combined[, id]
      }
      result$latent_variables <- apply(W_arr, c(2, 3), mean)
    }
  }

  for (species in 1:n_sp) {
    result$alpha_estimates[species]       <- mean(beta0_mat[, species], na.rm = TRUE)
    result$alpha_standard_errors[species] <- sd(beta0_mat[, species], na.rm = TRUE)
    result$beta_estimates[species]        <- mean(beta_arr[, species, 1], na.rm = TRUE)
    result$beta_standard_errors[species]  <- sd(beta_arr[, species, 1], na.rm = TRUE)
    if (fit_gamma) {
      result$gamma_estimates[species]       <- mean(beta_arr[, species, 2], na.rm = TRUE)
      result$gamma_standard_errors[species] <- sd(beta_arr[, species, 2], na.rm = TRUE)
    }
  }

  result$ess_alpha <- mean(neff_values[grep("^beta0\\[", neff_names)], na.rm = TRUE)
  result$ess_beta  <- mean(neff_values[grep("^beta\\[[0-9]+,1\\]", neff_names)], na.rm = TRUE)
  if (fit_gamma) {
    result$ess_gamma <- mean(neff_values[grep("^beta\\[[0-9]+,2\\]", neff_names)], na.rm = TRUE)
  }

  result$convergence <- attrs$final.params$converged

  unlink(model_file)

  return(result)
}

## JSDM.6 ######################################################################

fit_JSDM.6 <- function(Y,
                       covariates,
                       n_lv,
                       seed) {

  set.seed(seed)

  Y     <- as.matrix(Y)
  X_mat <- do.call(cbind, lapply(covariates, as.matrix))
  if (is.vector(X_mat)) {
    X_mat <- matrix(X_mat, ncol = 1)
  }

  n_obs    <- nrow(Y)
  n_sp     <- ncol(Y)
  n_covars <- ncol(X_mat)

  fit_gamma <- (length(covariates) >= 2)

  result <- results_initialization(n_sp, n_lv, fit_gamma = fit_gamma)
  
  if (n_lv > 0) {
    result$sigma_W_estimate <- rep(NA, n_lv)
    result$sigma_W_se       <- rep(NA, n_lv)
  }

  if (n_lv > 0) {
    model_str <- "
      model {
        for(i in 1:n_obs) {
          for(j in 1:n_sp) {
            eta[i, j] <- inprod(lambda[j, ], W[i, ]) + inprod(beta[j, ], X[i, ])
            probit(p_y[i, j]) <- beta0[j] + eta[i, j]
            y[i, j] ~ dbern(p_y[i, j])
          }
        }
        for(k in 1:n_lv) {
          tau_W[k] ~ dgamma(1, 1)
          sigma_W[k] <- pow(tau_W[k], -0.5)
        }
        for(i in 1:n_obs) {
          for(k in 1:n_lv) {
            W[i,k] ~ dnorm(0, tau_W[k])
          }
        }
        for(j in 1:n_sp) {
          beta0[j] ~ dnorm(0, 0.1)
          for(m in 1:n_covars) {
            beta[j, m] ~ dnorm(0, 0.1)
          }
        }
        for(i in 1:(n_lv - 1)) {
          for(j in (i + 1):n_lv) {
            lambda[i, j] <- 0
          }
        }
        for(i in 1:n_lv) {
          lambda[i, i] <- 1
        }
        for(i in 2:n_lv) {
          for(j in 1:(i - 1)) {
            lambda[i, j] ~ dnorm(0, 0.1)
          }
        }
        for(i in (n_lv + 1):n_sp) {
          for(j in 1:n_lv) {
            lambda[i, j] ~ dnorm(0, 0.1)
          }
        }
      }"
  } else {
    model_str <- "
      model {
        for(i in 1:n_obs) {
          for(j in 1:n_sp) {
            probit(p_y[i, j]) <- beta0[j] + inprod(beta[j, ], X[i, ])
            y[i, j] ~ dbern(p_y[i, j])
          }
        }
        for(j in 1:n_sp) {
          beta0[j] ~ dnorm(0, 0.1)
          for(m in 1:n_covars) {
            beta[j, m] ~ dnorm(0, 0.1)
          }
        }
      }"
  }

  data_list <- list(y        = Y,
                    X        = X_mat,
                    n_obs    = n_obs,
                    n_sp     = n_sp,
                    n_covars = n_covars)
  
  if (n_lv > 0) {
    data_list$n_lv <- n_lv
  }

  make_inits <- function() {
    inits <- list(beta0 = rnorm(n_sp, 0, 0.5),
                  beta  = matrix(rnorm(n_sp * n_covars, 0, 0.3),
                                 n_sp, n_covars))
    if (n_lv > 0) {
      inits$W     <- matrix(rnorm(n_obs * n_lv, 0, 0.5), n_obs, n_lv)
      inits$tau_W <- rgamma(n_lv, shape = 1, rate = 1)
    }
    return(inits)
  }
  
  mcmc_params <- DEFAULT_PARAMS$mcmc
  
  inits_list <- replicate(mcmc_params$n_chains, make_inits(), simplify = FALSE)

  monitor_params <- c("beta0", 
                      "beta")
  if (n_lv > 0) {
    free_lambda <- c()
    if (n_lv >= 2) {
      for (i in 2:n_lv) {
        for (j in 1:(i - 1)) {
          free_lambda <- c(free_lambda, sprintf("lambda[%d,%d]", i, j))
        }
      }
    }
    if (n_sp > n_lv) {
      for (i in (n_lv + 1):n_sp) {
        for (j in 1:n_lv) {
          free_lambda <- c(free_lambda, sprintf("lambda[%d,%d]", i, j))
        }
      }
    }
    monitor_params <- c(monitor_params,
                        "W", 
                        free_lambda, 
                        "sigma_W")
  }

  random_id <- format(Sys.time(), "%Y%m%d_%H%M%S")

  model_file <- file.path(tempdir(),
                          sprintf("jsdm6_%s_%04d.txt", random_id, sample.int(9999, 1)))
  writeLines(model_str, model_file)

  start_time <- Sys.time()

  out <- runMCMCbtadjust::runMCMC_btadjust(MCMC_language = "Jags",
                                           code          = model_file,
                                           data          = data_list,
                                           inits         = inits_list,
                                           params        = monitor_params,
                                           params.conv   = c("beta0", "beta"),
                                           niter.min     = mcmc_params$n_iter_min,
                                           niter.max     = Inf,
                                           nburnin.min   = mcmc_params$n_burnin_min,
                                           nburnin.max   = Inf,
                                           thin.min      = mcmc_params$thin_min,
                                           thin.max      = Inf,
                                           Nchains       = mcmc_params$n_chains,
                                           conv.max      = mcmc_params$conv_max,
                                           neff.min      = mcmc_params$n_eff_min,
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

  neff_values <- attrs$final.diags$neff
  neff_names  <- rownames(neff_values)

  beta0_cols <- grep("^beta0\\[", colnames(combined))
  beta0_mat  <- matrix(NA, nrow = nrow(combined), ncol = n_sp)
  for (id in beta0_cols) {
    j <- as.numeric(gsub("beta0\\[|\\]", "", colnames(combined)[id]))
    beta0_mat[, j] <- combined[, id]
  }

  beta_cols <- grep("^beta\\[", colnames(combined))
  beta_arr  <- array(NA, dim = c(nrow(combined), n_sp, n_covars))
  for (id in beta_cols) {
    col_name      <- colnames(combined)[id]
    param_indices <- as.numeric(strsplit(gsub("beta\\[|\\]", "", col_name), ",")[[1]])
    beta_arr[, param_indices[1], param_indices[2]] <- combined[, id]
  }

  if (n_lv > 0) {
    W_cols <- grep("^W\\[", colnames(combined))
    if (length(W_cols) > 0) {
      W_arr <- array(NA, dim = c(nrow(combined), n_obs, n_lv))
      for (id in W_cols) {
        col_name      <- colnames(combined)[id]
        param_indices <- as.numeric(strsplit(gsub("W\\[|\\]", "", col_name), ",")[[1]])
        W_arr[, param_indices[1], param_indices[2]] <- combined[, id]
      }
      result$latent_variables <- apply(W_arr, c(2, 3), mean)
    }

    s_cols <- grep("^sigma_W\\[", colnames(combined))
    if (length(s_cols) > 0) {
      s_mat <- matrix(NA, nrow = nrow(combined), ncol = n_lv)
      for (id in s_cols) {
        k <- as.numeric(gsub("sigma_W\\[|\\]", "", colnames(combined)[id]))
        s_mat[, k] <- combined[, id]
      }
      result$sigma_W_estimate <- apply(s_mat, 2, mean, na.rm = TRUE)
      result$sigma_W_se       <- apply(s_mat, 2, sd, na.rm = TRUE)
    }

    lambda_cols <- grep("^lambda\\[", colnames(combined))
    lambda_mean <- matrix(NA, nrow = n_sp, ncol = n_lv)
    if (length(lambda_cols) > 0) {
      for (id in lambda_cols) {
        col_name      <- colnames(combined)[id]
        param_indices <- as.numeric(strsplit(gsub("lambda\\[|\\]", "", col_name), ",")[[1]])
        lambda_mean[param_indices[1], param_indices[2]] <- mean(combined[, id], na.rm = TRUE)
      }
    }
    for (k in 1:min(n_lv, n_sp)) {
      lambda_mean[k, k] <- 1.0
    }
    if (n_lv >= 2) {
      for (i in 1:(n_lv - 1)) {
        for (j in (i + 1):n_lv) {
          lambda_mean[i, j] <- 0.0
        }
      }
    }
    result$loadings <- lambda_mean
  }

  for (species in 1:n_sp) {
    result$alpha_estimates[species]       <- mean(beta0_mat[, species], na.rm = TRUE)
    result$alpha_standard_errors[species] <- sd(beta0_mat[, species], na.rm = TRUE)
    result$beta_estimates[species]        <- mean(beta_arr[, species, 1], na.rm = TRUE)
    result$beta_standard_errors[species]  <- sd(beta_arr[, species, 1], na.rm = TRUE)
    if (fit_gamma) {
      result$gamma_estimates[species]       <- mean(beta_arr[, species, 2], na.rm = TRUE)
      result$gamma_standard_errors[species] <- sd(beta_arr[, species, 2], na.rm = TRUE)
    }
  }

  result$ess_alpha <- mean(neff_values[grep("^beta0\\[", neff_names)], na.rm = TRUE)
  result$ess_beta  <- mean(neff_values[grep("^beta\\[[0-9]+,1\\]", neff_names)], na.rm = TRUE)
  if (fit_gamma) {
    result$ess_gamma <- mean(neff_values[grep("^beta\\[[0-9]+,2\\]", neff_names)], na.rm = TRUE)
  }

  result$convergence <- attrs$final.params$converged

  unlink(model_file)

  return(result)
}

## DATA ANALYSIS ###############################################################

run_JSDM <- function(seed,
                     results_dir,
                     scenario_id,
                     methods = 1:6,
                     params = NULL) {

  if (is.null(params)) {
    params <- DEFAULT_PARAMS
  }

  scenarios <- build_scenarios(params$n_sp)
  scenario  <- scenarios[[scenario_id]]

  data <- simulate_data(seed     = seed,
                        scenario = scenario,
                        params   = params)

  covariates <- select_covariates(data, params$missing_covariate)

  results <- list(seed              = seed,
                  timestamp         = Sys.time(),
                  scenario_id       = scenario_id,
                  scenario          = scenario,
                  parameters        = list(n_obs             = params$n_obs,
                                           n_sp              = params$n_sp,
                                           n_covars          = params$n_covars,
                                           n_lv            = params$n_lv,
                                           missing_covariate = params$missing_covariate,
                                           B                 = data$B))

  total_start <- Sys.time()

  if (1 %in% methods) {
    results$JSDM.1 <- fit_JSDM.1(Y          = data$Y,
                                 covariates = covariates,
                                 n_lv     = params$n_lv,
                                 seed       = seed)
  }

  if (2 %in% methods) {
    results$JSDM.2 <- fit_JSDM.2(Y          = data$Y,
                                 covariates = covariates,
                                 n_lv     = params$n_lv,
                                 seed       = seed)
  }

  if (3 %in% methods) {
    results$JSDM.3 <- fit_JSDM.3(Y          = data$Y,
                                 covariates = covariates,
                                 n_lv     = params$n_lv,
                                 seed       = seed)
  }

  if (4 %in% methods) {
    results$JSDM.4 <- fit_JSDM.4(Y          = data$Y,
                                 covariates = covariates,
                                 n_lv     = params$n_lv,
                                 seed       = seed)
  }

  if (5 %in% methods) {
    results$JSDM.5 <- fit_JSDM.5(Y          = data$Y,
                                 covariates = covariates,
                                 n_lv     = params$n_lv,
                                 seed       = seed)
  }

  if (6 %in% methods) {
    results$JSDM.6 <- fit_JSDM.6(Y          = data$Y,
                                 covariates = covariates,
                                 n_lv     = params$n_lv,
                                 seed       = seed)
  }

  results$total_computation_time <- as.numeric(difftime(Sys.time(),
                                                        total_start,
                                                        units = "secs"))

  mc_tag <- if (params$missing_covariate) "MC" else "noMC"
  lv_tag <- sprintf("%dLV", params$n_lv)

  if (length(methods) == 1) {
    output_file <- file.path(results_dir,
                             sprintf("JSDM%d_S%d_%s_%s_seed%06d.RData", methods,
                                     scenario_id, mc_tag, lv_tag, seed))
  } else {
    method_str <- paste(methods, collapse = "_")
    output_file <- file.path(results_dir,
                             sprintf("JSDM_S%d_%s_%s_seed%06d_methods%s.RData",
                                     scenario_id, mc_tag, lv_tag, seed, method_str))
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

params <- DEFAULT_PARAMS

if (!is.null(n_lv)) {
  params$n_lv <- n_lv
}

if (!is.null(missing_covariate)) {
  params$missing_covariate <- missing_covariate
}

run_JSDM(seed        = seed,
         results_dir  = results_dir,
         scenario_id = scenario_id,
         methods     = methods,
         params      = params)

################################################################################
