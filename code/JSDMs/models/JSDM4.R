################################################################################
##                                                                            ##
##                         JSDM.4 FITTED WITH nimble                          ##
##                           (JSDMs/models/JSDM4.R)                           ##
##                                                                            ##
##       Addressing Missing Covariates in Species Distribution Models:        ##
##  Inferential Impacts and Mitigation via Joint Species Distribution Models  ##
##                                                                            ##
##                  Arthur F. Rossignol & Frédéric Gosselin                   ##
##                                                                            ##
##                                    2026                                    ##
##                                                                            ##
################################################################################

fit_JSDM.4 <- function(Y,
                       covariates,
                       n_lv,
                       mcmc_params,
                       seed) {

  set.seed(seed)

  Y     <- as.matrix(Y)
  X_mat <- cbind(1, do.call(cbind, lapply(covariates, as.matrix)))

  n_obs    <- nrow(Y)
  n_sp     <- ncol(Y)
  n_covars <- ncol(X_mat)

  fit_gamma <- (length(covariates) >= 2)

  result <- results_initialization(n_sp, n_lv)

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
  }
  else {
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

    inits <- list(meanenvsppar = rnorm(n_covars, 0, 0.1),
                  V            = diag(n_covars) * 0.5,
                  envsppar     = matrix(rnorm(n_covars * n_sp, 0, 0.2),
                                        n_covars, n_sp),
                  CL3          = rnorm(nobs, 0, 0.5))

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
                          loading_info$indices_id_s[i, 2]] <-
            inits$lambdapar_raw_s_np[i]
        }
      }

      if (n_lv > 1) {
        inits$probslambdapar    <- rep(0.9, n_lv - 1)
        inits$sigmalambdapar    <- rep(1, n_lv)
        inits$sigmalambdapar[1] <- inits$sigmalambdapar0
        for (i in 2:n_lv) {
          inits$sigmalambdapar[i] <- inits$sigmalambdapar[i - 1] *
                                     inits$probslambdapar[i - 1]
        }
      }
      else {
        inits$sigmalambdapar <- inits$sigmalambdapar0
      }

      if (loading_info$runn_id_u) {
        inits$lambdapar_raw_u_np <- rnorm(loading_info$n_id_u, 0, 0.2)
        inits$phibb              <- rep(1.5, loading_info$n_id_u)
        for (i in 1:loading_info$n_id_u) {
          inits$lambdapar[loading_info$indices_id_u[i, 1],
                          loading_info$indices_id_u[i, 2]] <-
            inits$lambdapar_raw_u_np[i]
        }
      }
    }

    return(inits)
  }

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

  time_max <- 3600 * mcmc_params$time_max_hours

  mcmc_control <- list(time.max                   = time_max,
                       round.thinmult             = TRUE,
                       print.diagnostics          = TRUE,
                       Ncycles.target             = 2,
                       check.convergence.firstrun = TRUE,
                       convtype                   = "Gelman",
                       seed                       = seed)

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
                                           control       = mcmc_control,
                                           control.MCMC  = list(parallelize = TRUE))

  result$computation_time <- as.numeric(difftime(Sys.time(),
                                                 start_time,
                                                 units = "secs"))

  combined <- do.call(rbind, out)
  attrs    <- attributes(out)

  neff_values <- attrs$final.diags$neff
  neff_names  <- names(neff_values)

  conv_values <- attrs$final.diags$conv
  conv_names  <- rownames(conv_values)
  coef_rows   <- grep("^envsppar\\[", conv_names)

  if (length(coef_rows) > 0) {
    result$psrf_max <- max(conv_values[coef_rows, "Point est."], na.rm = TRUE)
  }

  envsp_cols <- grep("^envsppar\\[", colnames(combined))

  if (length(envsp_cols) > 0) {
    envsp_arr <- array(NA, dim = c(nrow(combined), n_covars, n_sp))
    for (id in envsp_cols) {
      col_name      <- colnames(combined)[id]
      param_indices <- as.numeric(strsplit(gsub("envsppar\\[|\\]", "", col_name), ",")[[1]])
      envsp_arr[, param_indices[1], param_indices[2]] <- combined[, id]
    }

    for (species in 1:n_sp) {
      result$estimates$alpha[species]       <- mean(envsp_arr[, 1, species], na.rm = TRUE)
      result$standard_errors$alpha[species] <- sd(envsp_arr[, 1, species], na.rm = TRUE)

      if (n_covars >= 2) {
        result$estimates$beta[species]       <- mean(envsp_arr[, 2, species], na.rm = TRUE)
        result$standard_errors$beta[species] <- sd(envsp_arr[, 2, species], na.rm = TRUE)
      }

      if (fit_gamma && n_covars >= 3) {
        result$estimates$gamma[species]       <- mean(envsp_arr[, 3, species], na.rm = TRUE)
        result$standard_errors$gamma[species] <- sd(envsp_arr[, 3, species], na.rm = TRUE)
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
      result$scores <- apply(u_arr, c(2, 3), mean)
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

  result$converged <- attrs$final.params$converged

  return(result)
}

################################################################################
