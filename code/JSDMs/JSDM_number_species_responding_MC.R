################################################################################
##                                                                            ##
##   JSDMs WITH VARYING NUMBER OF SPECIES RESPONDING TO A MISSING COVARIATE   ##
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

args        <- commandArgs(trailingOnly = TRUE)
task_id     <- as.integer(args[1])
results_dir <- args[2]

## PARAMETERS ##################################################################

DEFAULT_PARAMS <- list(n_obs      = 1000,
                       n_sp       = 10,
                       n_lv       = 1,
                       true_alpha = rep(1, 10),
                       true_beta  = rep(1, 10),
                       true_gamma = rep(1, 10),
                       mcmc       = list(n_iter_min   = 25000,
                                         n_burnin_min = 25000,
                                         thin_min     = 10,
                                         n_chains     = 3,
                                         conv_max     = 1.05,
                                         neff_min     = 1000,
                                         time_max     = 24))

## DATA GENERATING PROCESS #####################################################

simulate_data <- function(replicate_id,
                          k_responding,
                          params) {
  
  seed <- replicate_id
  set.seed(seed)
  
  n_obs <- params$n_obs
  n_sp  <- params$n_sp
  
  X1 <- rnorm(n, mean = 0, sd = 1)
  X2 <- rnorm(n, mean = 0, sd = 1)
  X  <- cbind(intercept = 1, X1 = X1, X2 = X2)
  
  gamma_active <- params$true_gamma
  if (k_responding < n_sp) {
    gamma_active[(k_responding + 1):n_sp] <- 0
  }
  
  B <- rbind(params$true_alpha,
             params$true_beta,
             gamma_active)
  rownames(B) <- c("alpha", "beta", "gamma")
  colnames(B) <- paste0("sp", 1:n_sp)
  
  M <- X %*% B
  Y <- matrix(NA_integer_, 
              nrow = n_obs, 
              ncol = n_sp,
              dimnames = list(NULL, paste0("sp", 1:n_sp)))
  for (j in 1:n_sp) {
    Y[, j] <- rbinom(n_obs, size = 1, prob = pnorm(M[, j]))
  }
  
  list(Y            = Y,
       X1           = X1,
       X2           = X2,
       X            = X,
       B            = B,
       params       = params,
       k_responding = k_responding,
       seed         = seed)
}

## RESULT TEMPLATE INITIALIZATION ##############################################

initialize_results <- function(n_sp,
                               n_lv) {
  
  list(
    # fixed effects
    alpha_estimates       = rep(NA_real_, n_sp),
    alpha_standard_errors = rep(NA_real_, n_sp),
    beta_estimates        = rep(NA_real_, n_sp),
    beta_standard_errors  = rep(NA_real_, n_sp),
    gamma_estimates       = rep(NA_real_, n_sp),
    gamma_standard_errors = rep(NA_real_, n_sp),
    # latent structure
    n_lv             = n_lv,
    latent_variables = NULL,
    loadings         = NULL,
    sigma_lv         = NA_real_,                   # JSDM.1 only
    sigma_W_estimate = rep(NA_real_, n_lv),        # JSDM.6 only
    sigma_W_se       = rep(NA_real_, n_lv),        # JSDM.6 only
    # convergence diagnostics
    convergence      = FALSE,
    ess_alpha        = NA_real_,
    ess_beta         = NA_real_,
    ess_gamma        = NA_real_,
    psrf_max         = NA_real_,             # Gelman-Rubin (JSDM.3–6)
    geweke_max       = NA_real_,             # Geweke (JSDM.2)
    # computation time
    computation_time = NA_real_)
}

## JSDM.1 ######################################################################

fit_JSDM.1 <- function(Y,
                       X1,
                       X2,
                       n_lv,
                       seed) {
  
  set.seed(seed)
  
  n_sp <- ncol(Y)
  
  result <- initialize_results(n_sp, n_lv)
  
  start_time <- Sys.time()
  
  fit <- gllvm::gllvm(y             = Y,
                      X             = data.frame(X1 = X1,
                                                 X2 = X2),
                      family        = binomial(link = "probit"),
                      method        = "LA",
                      num.lv        = n_lv,
                      formula       = ~ 1 + X1 + X2,
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
  result$alpha_standard_errors <- as.vector(fit$sd$beta0)
  Xcoef_estimates              <- fit$params$Xcoef
  Xcoef_standard_errors        <- fit$sd$Xcoef
  result$beta_estimates        <- as.vector(Xcoef_estimates[, 1])
  result$beta_standard_errors  <- as.vector(Xcoef_standard_errors[, 1])
  result$gamma_estimates       <- as.vector(Xcoef_estimates[, 2])
  result$gamma_standard_errors <- as.vector(Xcoef_standard_errors[, 2])
  
  result$latent_variables <- fit$lvs
  result$n_lv             <- ncol(result$latent_variables)
  result$loadings         <- fit$params$theta
  result$sigma_lv         <- if (!is.null(fit$params$sigma.lv)) 
                               as.numeric(fit$params$sigma.lv) else NA_real_
  
  result$convergence <- fit$convergence
  
  result
}

## JSDM.2 ######################################################################

fit_JSDM.2 <- function(Y,
                       X1,
                       X2,
                       n_lv,
                       seed) {
  
  set.seed(seed)
  
  n_obs <- nrow(Y)
  n_sp  <- ncol(Y)
  
  result <- initialize_results(n_sp, n_lv)
  
  start_time <- Sys.time()
  
  fit <- jSDM::jSDM_binomial_probit(burnin        = 100000,
                                    mcmc          = 100000,
                                    thin          = 5,
                                    presence_data = Y,
                                    site_formula  = ~ X1 + X2,
                                    site_data     = data.frame(X1 = X1,
                                                               X2 = X2),
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
  
  n_covars <- 3
  
  mean_beta   <- matrix(NA, nrow = n_sp, ncol = n_covars)
  se_beta     <- matrix(NA, nrow = n_sp, ncol = n_covars)
  mean_lambda <- matrix(NA, nrow = n_sp, ncol = n_lv)
  
  for (j in 1:n_sp) {
    sp_samples <- as.matrix(fit$mcmc.sp[[j]])
    
    for (cv in 1:n_covars) {
      mean_beta[j, cv] <- mean(sp_samples[, cv])
      se_beta[j, cv]   <- sd(sp_samples[, cv])
    }
    
    for (l in 1:n_lv) {
      mean_lambda[j, l] <- mean(sp_samples[, n_covars + l])
    }
  }
  
  result$alpha_estimates       <- mean_beta[, 1]
  result$alpha_standard_errors <- se_beta[, 1]
  result$beta_estimates        <- mean_beta[, 2]
  result$beta_standard_errors  <- se_beta[, 2]
  result$gamma_estimates       <- mean_beta[, 3]
  result$gamma_standard_errors <- se_beta[, 3]
  
  result$loadings <- mean_lambda
  
  if (!is.null(fit$mcmc.latent) && length(fit$mcmc.latent) > 0) {
    W_mean <- matrix(NA, nrow = n, ncol = n_lv)
    for (l in 1:n_lv) {
      lv_name <- paste0("lv_", l)
      if (lv_name %in% names(fit$mcmc.latent)) {
        lv_samples <- as.matrix(fit$mcmc.latent[[lv_name]])
        W_mean[, l] <- colMeans(lv_samples)
      }
    }
    result$latent_variables <- W_mean
    result$n_lv             <- n_lv
  }
  
  ess_alpha <- numeric(n_sp)
  ess_beta  <- numeric(n_sp)
  ess_gamma <- numeric(n_sp)
  
  for (j in 1:n_sp) {
    sp_mcmc <- fit$mcmc.sp[[j]]
    ess_all <- coda::effectiveSize(sp_mcmc)
    if (length(ess_all) >= 1) {
      ess_alpha[j] <- ess_all[1]
    }
    if (length(ess_all) >= 2) {
      ess_beta[j]  <- ess_all[2]
    }
    if (length(ess_all) >= 3) {
      ess_gamma[j] <- ess_all[3]
    }
  }
  
  result$ess_alpha <- mean(ess_alpha, na.rm = TRUE)
  result$ess_beta  <- mean(ess_beta, na.rm = TRUE)
  result$ess_gamma <- mean(ess_gamma, na.rm = TRUE)
  
  geweke_vals <- numeric(n_sp)
  for (j in 1:n_sp) {
    gw <- coda::geweke.diag(fit$mcmc.sp[[j]])$z
    if (is.numeric(gw) && length(gw) >= 3) {
      geweke_vals[j] <- max(abs(gw[1:3]), na.rm = TRUE)
    }
  }
  
  max_geweke         <- max(geweke_vals, na.rm = TRUE)
  result$geweke_max  <- max_geweke
  result$convergence <- max_geweke < 2.0
  
  result
}

## JSDM.3 ######################################################################

fit_JSDM.3 <- function(Y,
                       X1,
                       X2,
                       n_lv,
                       seed) {
  
  set.seed(seed)
  
  n_sp  <- ncol(Y)
  n_obs <- nrow(Y)
  
  result <- initialize_results(n_sp, n_lv)
  
  studyDesign <- data.frame(site = as.factor(1:n_obs))
  rL          <- Hmsc::HmscRandomLevel(units = studyDesign$site)
  rL$nfMin    <- n_lv
  rL$nfMax    <- n_lv
  
  model <- Hmsc::Hmsc(Y           = Y,
                      XData       = data.frame(X1 = X1, X2 = X2),
                      XFormula    = ~ 1 + X1 + X2,
                      XScale      = FALSE,
                      studyDesign = studyDesign,
                      ranLevels   = list(site = rL),
                      distr       = "probit")
  
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
  result$gamma_estimates <- postBeta$mean[3, ]
  
  mpost <- Hmsc::convertToCodaObject(fit)
  
  beta_samples  <- do.call(rbind, mpost$Beta)
  intercept_id <- seq(1, n_sp * 3, by = 3)
  x1_id        <- seq(2, n_sp * 3, by = 3)
  x2_id        <- seq(3, n_sp * 3, by = 3)
  
  result$alpha_standard_errors <- apply(beta_samples[, intercept_id, drop = FALSE],
                                        2, sd)
  result$beta_standard_errors  <- apply(beta_samples[, x1_id, drop = FALSE],
                                        2, sd)
  result$gamma_standard_errors <- apply(beta_samples[, x2_id, drop = FALSE],
                                        2, sd)
  
  postEta    <- Hmsc::getPostEstimate(fit, parName = "Eta")
  postLambda <- Hmsc::getPostEstimate(fit, parName = "Lambda")
  
  result$latent_variables <- postEta$mean
  result$n_lv             <- ncol(result$latent_variables)
  result$loadings         <- postLambda$mean
  
  ess_values       <- effectiveSize(mpost$Beta)
  result$ess_alpha <- mean(ess_values[intercept_id], na.rm = TRUE)
  result$ess_beta  <- mean(ess_values[x1_id], na.rm = TRUE)
  result$ess_gamma <- mean(ess_values[x2_id], na.rm = TRUE)
  
  if (length(mpost$Beta) > 1) {
    psrf_values        <- gelman.diag(mpost$Beta,
                                      multivariate = FALSE)$psrf[, "Point est."]
    result$psrf_max    <- max(psrf_values, na.rm = TRUE)
    result$convergence <- result$psrf_max < 1.05
  }
  
  result
}

## JSDM.4 ######################################################################

fit_JSDM.4 <- function(Y,
                       X1,
                       X2,
                       n_lv,
                       seed) {
  
  set.seed(seed)
  
  Y     <- as.matrix(Y)
  X_mat <- cbind(1, as.matrix(X1), as.matrix(X2))
  
  n_obs    <- nrow(Y)
  n_sp     <- ncol(Y)
  n_covars <- ncol(X_mat)
  
  result <- initialize_results(n_sp, n_lv)
  
  modelCode <- nimble::nimbleCode({
    
    for (i in 1:nvar) {
      meanenvsppar[i] ~ dnorm(0, 1)
    }
    V[1:nvar, 1:nvar] ~ dinvwish(S = V0[1:nvar, 1:nvar], df = f0)
    
    for (j in 1:nspecies) {
      envsppar[1:nvar, j] ~ dmnorm(meanenvsppar[1:nvar],
                                   cov = V[1:nvar, 1:nvar])
    }
    
    sigmalambdapar0 ~ dgamma(a[1], b[1])
    sigmalambdapar[1] <- sigmalambdapar0
    
    if (nlvbigger1) {
      for (i in 1:(nlv - 1)) {
        probslambdapar[i] ~ dgamma(a[2], b[2])
        sigmalambdapar[i + 1] <- sigmalambdapar[i] * probslambdapar[i]
      }
    }
    
    if (nlvbigger0) {
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
        lambdapar_raw_s_np[i] ~ dnorm(0,
                                      sigmalambdapar[indices_id_s[i, 2]] *
                                        phib[i])
      }
      
      for (i in 1:n_id_s) {
        lambdapar[indices_id_s[i, 1], indices_id_s[i, 2]] <-
          meanlambdapar[indices_id_s[i, 2]] + lambdapar_raw_s_np[i]
        lambdapar_abs_s_np[i] <-
          abs(lambdapar[indices_id_s[i, 1], indices_id_s[i, 2]])
      }
      
      if (runn_id_u) {
        for (i in 1:n_id_u) {
          phibb[i] ~ dgamma(nu / 2, nu / 2)
          lambdapar_raw_u_np[i] ~ dnorm(0,
                                        sigmalambdapar[indices_id_u[i, 2]] *
                                          phibb[i])
          lambdapar[indices_id_u[i, 1], indices_id_u[i, 2]] <-
            lambdapar_raw_u_np[i]
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
        CLel.nlv[i] <- sum(lambdapar[beta0lambda[i], 1:nlv] *
                             upar[alpha[i], 1:nlv])
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
  
  ModelConsts <- list(nobs         = n_obs,
                      nsites       = n_obs,
                      nspecies     = n_sp,
                      nvar         = n_covars,
                      nlv          = n_lv,
                      nlvbigger0   = (n_lv > 0),
                      nlvbigger1   = (n_lv > 1),
                      nvarbigger1  = (n_covars > 1),
                      beta0        = beta0,
                      beta0lambda  = beta0lambda,
                      alpha        = alpha,
                      a            = c(50, 50),
                      b            = c(1, 1),
                      nu           = 3,
                      V0           = diag(n_covars),
                      f0           = n_covars + 1,
                      indices_id_s = loading_info$indices_id_s,
                      n_id_s       = loading_info$n_id_s,
                      sum_id_s     = loading_info$sum_id_s,
                      indices_id_u = loading_info$indices_id_u,
                      n_id_u       = loading_info$n_id_u,
                      runn_id_u    = loading_info$runn_id_u)
  
  ModelData <- list(Y   = Y_long,
                    env = env_long)
  
  make_inits <- function() {
    inits <- list(meanenvsppar    = rnorm(n_covars, 0, 0.1),
                  V               = diag(n_covars) * 0.5,
                  envsppar        = matrix(rnorm(n_covars * n_sp, 0, 0.2),
                                           n_covars, n_sp),
                  sigmalambdapar0 = 1.0,
                  CL3             = rnorm(nobs, 0, 0.5))
    
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
                        loading_info$indices_id_s[i, 2]] <-
          inits$lambdapar_raw_s_np[i]
      }
    }
    
    inits$sigmalambdapar <- inits$sigmalambdapar0
    
    if (loading_info$runn_id_u) {
      inits$lambdapar_raw_u_np <- rnorm(loading_info$n_id_u, 0, 0.2)
      inits$phibb              <- rep(1.5, loading_info$n_id_u)
      for (i in 1:loading_info$n_id_u) {
        inits$lambdapar[loading_info$indices_id_u[i, 1],
                        loading_info$indices_id_u[i, 2]] <-
          inits$lambdapar_raw_u_np[i]
      }
    }
    inits
  }
  
  mcmc_params <- DEFAULT_PARAMS$mcmc
  
  inits_list <- replicate(mcmc_params$n_chains, make_inits(), simplify = FALSE)
  
  params <- c("meanenvsppar",
              "envsppar", 
              "V",
              "upar", 
              "lambdapar", 
              "sigmalambdapar")
  
  start_time <- Sys.time()
  
  out <- runMCMCbtadjust::runMCMC_btadjust(
    MCMC_language = "Nimble",
    code          = modelCode,
    constants     = ModelConsts,
    data          = ModelData,
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
    neff.min      = mcmc_params$neff_min,
    control       = list(time.max                   = 3600 * mcmc_params$time_max,
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
  
  envsp_cols <- grep("^envsppar\\[", colnames(combined))
  if (length(envsp_cols) > 0) {
    envsp_arr <- array(NA, dim = c(nrow(combined), n_covars, n_sp))
    for (id in envsp_cols) {
      ij <- as.numeric(strsplit(gsub("envsppar\\[|\\]", "",
                                     colnames(combined)[id]), ",")[[1]])
      envsp_arr[, ij[1], ij[2]] <- combined[, id]
    }
    
    for (sp in 1:n_sp) {
      result$alpha_estimates[sp]       <- mean(envsp_arr[, 1, sp], na.rm = TRUE)
      result$alpha_standard_errors[sp] <- sd(envsp_arr[, 1, sp], na.rm = TRUE)
      result$beta_estimates[sp]        <- mean(envsp_arr[, 2, sp], na.rm = TRUE)
      result$beta_standard_errors[sp]  <- sd(envsp_arr[, 2, sp], na.rm = TRUE)
      result$gamma_estimates[sp]       <- mean(envsp_arr[, 3, sp], na.rm = TRUE)
      result$gamma_standard_errors[sp] <- sd(envsp_arr[, 3, sp], na.rm = TRUE)
    }
    
    result$ess_alpha <- mean(sapply(1:n_sp, function(sp)
      effectiveSize(as.mcmc(envsp_arr[, 1, sp]))), na.rm = TRUE)
    result$ess_beta <- mean(sapply(1:n_sp, function(sp)
      effectiveSize(as.mcmc(envsp_arr[, 2, sp]))), na.rm = TRUE)
    result$ess_gamma <- mean(sapply(1:n_sp, function(sp)
      effectiveSize(as.mcmc(envsp_arr[, 3, sp]))), na.rm = TRUE)
  }
  
  u_cols <- grep("^upar\\[", colnames(combined))
  if (length(u_cols) > 0) {
    u_arr <- array(NA, dim = c(nrow(combined), n_obs, n_lv))
    for (id in u_cols) {
      jk <- as.numeric(strsplit(gsub("upar\\[|\\]", "", 
                                     colnames(combined)[id]), ",")[[1]])
      u_arr[, jk[1], jk[2]] <- combined[, id]
    }
    result$latent_variables <- apply(u_arr, c(2, 3), mean)
  }
  
  lambda_cols <- grep("^lambdapar\\[", colnames(combined))
  if (length(lambda_cols) > 0) {
    lambda_arr <- array(NA, dim = c(nrow(combined), n_sp, n_lv))
    for (id in lambda_cols) {
      jk <- as.numeric(strsplit(gsub("lambdapar\\[|\\]", "", 
                                     colnames(combined)[id]), ",")[[1]])
      lambda_arr[, jk[1], jk[2]] <- combined[, id]
    }
    result$loadings <- apply(lambda_arr, c(2, 3), mean)
  }
  
  result$convergence <- attrs$final.params$converged
  
  result
}

## JSDM.5 ######################################################################

fit_JSDM.5 <- function(Y,
                       X1,
                       X2,
                       n_lv,
                       seed) {
  
  set.seed(seed)
  
  Y     <- as.matrix(Y)
  X_mat <- cbind(as.matrix(X1), as.matrix(X2))
  
  n_obs    <- nrow(Y)
  n_sp     <- ncol(Y)
  n_covars <- ncol(X_mat)
  
  result <- initialize_results(n_sp, n_lv)
  
  model_string <- "
    model {
      for(i in 1:n) {
        for(j in 1:p) {
          eta[i, j] <- inprod(lambda[j, ], W[i, ]) + inprod(beta[j, ], X[i, ])
          probit(p_y[i, j]) <- beta0[j] + eta[i, j]
          y[i, j] ~ dbern(p_y[i, j])
        }
      }
      for(i in 1:n) {
        for(k in 1:n_lv) {
          W[i, k] ~ dnorm(0, 1)
        }
      }
      for(j in 1:p) {
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
      for(i in (n_lv + 1):p) {
        for(j in 1:n_lv) {
          lambda[i, j] ~ dnorm(0, 0.1)
        }
      }
      for(j in 1:p) {
        for(m in 1:n_covars) {
          beta[j, m] ~ dnorm(0, 0.1)
        }
      }
    }"
  
  data_list <- list(y        = Y,
                    X        = X_mat,
                    n        = n_obs,
                    p        = n_sp,
                    n_covars = n_covars,
                    n_lv     = n_lv)
  
  make_inits <- function() {
    inits <- list(beta0 = rnorm(n_sp, 0, 0.5),
                  beta  = matrix(rnorm(n_sp * n_covars, 0, 0.3), n_sp, n_covars),
                  W     = matrix(rnorm(n_obs * n_lv, 0, 0.5), n_obs, n_lv))
    inits
  }
  
  mcmc_params <- DEFAULT_PARAMS$mcmc
  
  inits_list <- replicate(mcmc_params$n_chains, make_inits(), simplify = FALSE)
  
  monitor_params <- c("beta0", "beta", "lambda", "W")
  
  random_id  <- format(Sys.time(), "%Y%m%d_%H%M%S")
  model_file <- file.path(tempdir(),
                          sprintf("jsdm5_%s_%04d.txt", random_id, sample.int(9999, 1)))
  writeLines(model_string, model_file)
  
  start_time <- Sys.time()
  
  out <- runMCMCbtadjust::runMCMC_btadjust(
    MCMC_language = "Jags",
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
    neff.min      = mcmc_params$neff_min,
    control       = list(time.max                   = 3600 * mcmc_params$time_max,
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
  beta0_mat  <- matrix(NA, nrow = nrow(combined), ncol = n_sp)
  for (id in beta0_cols) {
    j <- as.numeric(gsub("beta0\\[|\\]", "", colnames(combined)[id]))
    beta0_mat[, j] <- combined[, id]
  }
  
  beta_cols <- grep("^beta\\[", colnames(combined))
  beta_arr  <- array(NA, dim = c(nrow(combined), n_sp, n_covars))
  for (id in beta_cols) {
    ij <- as.numeric(strsplit(gsub("beta\\[|\\]", "", 
                                   colnames(combined)[id]), ",")[[1]])
    beta_arr[, ij[1], ij[2]] <- combined[, id]
  }
  
  lambda_cols <- grep("^lambda\\[", colnames(combined))
  lambda_arr  <- array(NA, dim = c(nrow(combined), n_sp, n_lv))
  for (id in lambda_cols) {
    ij <- as.numeric(strsplit(gsub("lambda\\[|\\]", "", 
                                   colnames(combined)[id]), ",")[[1]])
    lambda_arr[, ij[1], ij[2]] <- combined[, id]
  }
  result$loadings <- apply(lambda_arr, c(2, 3), mean, na.rm = TRUE)
  
  W_cols <- grep("^W\\[", colnames(combined))
  if (length(W_cols) > 0) {
    W_arr <- array(NA, dim = c(nrow(combined), n_obs, n_lv))
    for (id in W_cols) {
      ik <- as.numeric(strsplit(gsub("W\\[|\\]", "", 
                                     colnames(combined)[id]), ",")[[1]])
      W_arr[, ik[1], ik[2]] <- combined[, id]
    }
    result$latent_variables <- apply(W_arr, c(2, 3), mean)
  }
  
  for (sp in 1:n_sp) {
    result$alpha_estimates[sp]       <- mean(beta0_mat[, sp], na.rm = TRUE)
    result$alpha_standard_errors[sp] <- sd(beta0_mat[, sp], na.rm = TRUE)
    result$beta_estimates[sp]        <- mean(beta_arr[, sp, 1], na.rm = TRUE)
    result$beta_standard_errors[sp]  <- sd(beta_arr[, sp, 1], na.rm = TRUE)
    result$gamma_estimates[sp]       <- mean(beta_arr[, sp, 2], na.rm = TRUE)
    result$gamma_standard_errors[sp] <- sd(beta_arr[, sp, 2], na.rm = TRUE)
  }
  
  ess_alpha        <- sapply(1:n_sp, function(sp) effectiveSize(as.mcmc(beta0_mat[, sp])))
  ess_beta         <- sapply(1:n_sp, function(sp) effectiveSize(as.mcmc(beta_arr[, sp, 1])))
  ess_gamma        <- sapply(1:n_sp, function(sp) effectiveSize(as.mcmc(beta_arr[, sp, 2])))
  result$ess_alpha <- mean(ess_alpha, na.rm = TRUE)
  result$ess_beta  <- mean(ess_beta, na.rm = TRUE)
  result$ess_gamma <- mean(ess_gamma, na.rm = TRUE)
  
  result$convergence <- attrs$final.params$converged
  
  unlink(model_file)
  
  result
}

## JSDM.6 ######################################################################

fit_JSDM.6 <- function(Y,
                       X1,
                       X2,
                       n_lv,
                       seed) {
  
  set.seed(seed)
  
  Y     <- as.matrix(Y)
  X_mat <- cbind(as.matrix(X1), as.matrix(X2))
  
  n_obs    <- nrow(Y)
  n_sp     <- ncol(Y)
  n_covars <- ncol(X_mat)
  
  result <- initialize_results(n_sp, n_lv)
  
  model_string <- "
    model {
      for(i in 1:n) {
        for(j in 1:p) {
          eta[i, j] <- inprod(lambda[j, ], W[i, ]) + inprod(beta[j, ], X[i, ])
          probit(p_y[i, j]) <- beta0[j] + eta[i, j]
          y[i, j] ~ dbern(p_y[i, j])
        }
      }
      for(k in 1:n_lv) {
        tau_W[k] ~ dgamma(1, 1)
        sigma_W[k] <- pow(tau_W[k], -0.5)
      }
      for(i in 1:n) {
        for(k in 1:n_lv) {
          W[i, k] ~ dnorm(0, tau_W[k])
        }
      }
      for(j in 1:p) {
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
      for(i in (n_lv + 1):p) {
        for(j in 1:n_lv) {
          lambda[i, j] ~ dnorm(0, 0.1)
        }
      }
    }"
  
  data_list <- list(y        = Y,
                    X        = X_mat,
                    n        = n_obs,
                    p        = n_sp,
                    n_covars = n_covars,
                    n_lv     = n_lv)
  
  make_inits <- function() {
    inits <- list(beta0 = rnorm(n_sp, 0, 0.5),
                  beta  = matrix(rnorm(n_sp * n_covars, 0, 0.3),
                                 n_sp, n_covars),
                  W     = matrix(rnorm(n_obs * n_lv, 0, 0.5), n_obs, n_lv),
                  tau_W = rgamma(n_lv, shape = 1, rate = 1))
    inits
  }
  
  mcmc_params <- DEFAULT_PARAMS$mcmc
  
  inits_list <- replicate(mcmc_params$n_chains, make_inits(), simplify = FALSE)
  
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
  monitor_params <- c("beta0", "beta", "W", free_lambda, "sigma_W", "tau_W")
  
  random_id  <- format(Sys.time(), "%Y%m%d_%H%M%S")
  model_file <- file.path(tempdir(),
                          sprintf("jsdm6_%s_%04d.txt", random_id, sample.int(9999, 1)))
  writeLines(model_string, model_file)
  
  start_time <- Sys.time()
  
  out <- runMCMCbtadjust::runMCMC_btadjust(
    MCMC_language = "Jags",
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
    neff.min      = mcmc_params$neff_min,
    control       = list(time.max                   = 3600 * mcmc_params$time_max,
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
  beta0_mat  <- matrix(NA, nrow = nrow(combined), ncol = n_sp)
  for (id in beta0_cols) {
    j <- as.numeric(gsub("beta0\\[|\\]", "", colnames(combined)[id]))
    beta0_mat[, j] <- combined[, id]
  }
  
  beta_cols <- grep("^beta\\[", colnames(combined))
  beta_arr  <- array(NA, dim = c(nrow(combined), n_sp, n_covars))
  for (id in beta_cols) {
    ij <- as.numeric(strsplit(gsub("beta\\[|\\]", "", 
                                   colnames(combined)[id]), ",")[[1]])
    beta_arr[, ij[1], ij[2]] <- combined[, id]
  }
  
  W_arr  <- NULL
  W_cols <- grep("^W\\[", colnames(combined))
  if (length(W_cols) > 0) {
    W_arr <- array(NA, dim = c(nrow(combined), n_obs, n_lv))
    for (id in W_cols) {
      ik <- as.numeric(strsplit(gsub("W\\[|\\]", "", 
                                     colnames(combined)[id]), ",")[[1]])
      W_arr[, ik[1], ik[2]] <- combined[, id]
    }
    result$latent_variables <- apply(W_arr, c(2, 3), mean)
  }
  
  ## Three-stage fallback for sigma_W
  s_mat  <- NULL
  s_cols <- grep("^sigma_W\\[", colnames(combined))
  if (length(s_cols) > 0) {
    s_mat <- matrix(NA, nrow = nrow(combined), ncol = n_lv)
    for (id in s_cols) {
      k <- as.numeric(gsub("sigma_W\\[|\\]", "", colnames(combined)[id]))
      s_mat[, k] <- combined[, id]
    }
  } else {
    t_cols <- grep("^tau_W\\[", colnames(combined))
    if (length(t_cols) > 0) {
      t_mat <- matrix(NA, nrow = nrow(combined), ncol = n_lv)
      for (id in t_cols) {
        k <- as.numeric(gsub("tau_W\\[|\\]", "", colnames(combined)[id]))
        t_mat[, k] <- combined[, id]
      }
      s_mat <- 1 / sqrt(t_mat)
    } else if (!is.null(W_arr)) {
      s_mat <- matrix(NA, nrow = nrow(combined), ncol = n_lv)
      for (k in seq_len(n_lv)) {
        s_mat[, k] <- apply(W_arr[, , k], 1, sd, na.rm = TRUE)
      }
    }
  }
  if (!is.null(s_mat)) {
    result$sigma_W_estimate <- apply(s_mat, 2, mean, na.rm = TRUE)
    result$sigma_W_se       <- apply(s_mat, 2, sd,   na.rm = TRUE)
  }
  
  lambda_cols <- grep("^lambda\\[", colnames(combined))
  lambda_mean <- matrix(NA, nrow = n_sp, ncol = n_lv)
  if (length(lambda_cols) > 0) {
    for (id in lambda_cols) {
      jk <- as.numeric(strsplit(gsub("lambda\\[|\\]", "", 
                                     colnames(combined)[id]), ",")[[1]])
      lambda_mean[jk[1], jk[2]] <- mean(combined[, id], na.rm = TRUE)
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
  
  for (sp in 1:n_sp) {
    result$alpha_estimates[sp]       <- mean(beta0_mat[, sp], na.rm = TRUE)
    result$alpha_standard_errors[sp] <- sd(beta0_mat[, sp], na.rm = TRUE)
    result$beta_estimates[sp]        <- mean(beta_arr[, sp, 1], na.rm = TRUE)
    result$beta_standard_errors[sp]  <- sd(beta_arr[, sp, 1], na.rm = TRUE)
    result$gamma_estimates[sp]       <- mean(beta_arr[, sp, 2], na.rm = TRUE)
    result$gamma_standard_errors[sp] <- sd(beta_arr[, sp, 2], na.rm = TRUE)
  }
  
  ess_alpha        <- sapply(1:n_sp, function(sp) effectiveSize(as.mcmc(beta0_mat[, sp])))
  ess_beta         <- sapply(1:n_sp, function(sp) effectiveSize(as.mcmc(beta_arr[, sp, 1])))
  ess_gamma        <- sapply(1:n_sp, function(sp) effectiveSize(as.mcmc(beta_arr[, sp, 2])))
  result$ess_alpha <- mean(ess_alpha, na.rm = TRUE)
  result$ess_beta  <- mean(ess_beta, na.rm = TRUE)
  result$ess_gamma <- mean(ess_gamma, na.rm = TRUE)
  
  result$convergence <- attrs$final.params$converged
  
  unlink(model_file)
  
  result
}

## RUN #########################################################################

k_responding <- ((task_id - 1) %% 10) + 1
replicate_id <- ((task_id - 1) %/% 10 %%  100) + 1
model_id     <- ((task_id - 1) %/% 1000) + 1

model_names <- c("JSDM.1", "JSDM.2", "JSDM.3", "JSDM.4", "JSDM.5", "JSDM.6")
model_name  <- model_names[model_id]

data <- simulate_data(replicate_id = replicate_id,
                      k_responding = k_responding,
                      params       = DEFAULT_PARAMS)

fit_functions <- list(fit_JSDM.1,
                      fit_JSDM.2,
                      fit_JSDM.3,
                      fit_JSDM.4,
                      fit_JSDM.5,
                      fit_JSDM.6)

model_result <- fit_functions[[model_id]](Y    = data$Y,
                                          X1   = data$X1,
                                          X2   = data$X2)

results <- list(task_id      = task_id,
                model_id     = model_id,
                model_name   = model_name,
                replicate_id = replicate_id,
                k_responding = k_responding,
                seed         = data$seed,
                parameters = list(n          = DEFAULT_PARAMS$n,
                                  n_sp       = DEFAULT_PARAMS$n_sp,
                                  n_lv       = DEFAULT_PARAMS$n_lv,
                                  true_alpha = data$B["alpha", ],
                                  true_beta  = data$B["beta", ],
                                  true_gamma = data$B["gamma", ]),
                model_result = model_result)

file_name <- file.path(results_dir, sprintf("%s_k%02d_rep%03d.RData",
                                            gsub("\\.", "", model_name),
                                            k_responding,
                                            replicate_id))

save(results, file = file_name)

################################################################################
