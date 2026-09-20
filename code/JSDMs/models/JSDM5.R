################################################################################
##                                                                            ##
##            JSDM.5 FITTED WITH runjags AND UNIT LATENT VARIABLES            ##
##                           (JSDMs/models/JSDM5.R)                           ##
##                                                                            ##
##       Addressing Missing Covariates in Species Distribution Models:        ##
##  Inferential Impacts and Mitigation via Joint Species Distribution Models  ##
##                                                                            ##
##                  Arthur F. Rossignol & Frédéric Gosselin                   ##
##                                                                            ##
##                                    2026                                    ##
##                                                                            ##
################################################################################

fit_JSDM.5 <- function(Y,
                       covariates,
                       n_lv,
                       mcmc_params,
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

  result <- results_initialization(n_sp, n_lv)

  if (n_lv > 0) {
    # Build the lower-triangular lambda constraint block conditionally:
    # JAGS ranges like 1:(n_lv-1) are never empty (1:0 = c(1,0)), so each
    # loop must be omitted entirely when its index set is empty.
    lambda_upper_str <- ""
    if (n_lv > 1) {
      lambda_upper_str <- "
        for(i in 1:(n_lv - 1)) {
          for(j in (i + 1):n_lv) {
            lambda[i, j] <- 0
          }
        }"
    }

    lambda_lower_str <- ""
    if (n_lv > 1) {
      lambda_lower_str <- "
        for(i in 2:n_lv) {
          for(j in 1:(i - 1)) {
            lambda[i, j] ~ dnorm(0, 0.1)
          }
        }"
    }

    lambda_free_str <- ""
    if (n_sp > n_lv) {
      lambda_free_str <- "
        for(i in (n_lv + 1):n_sp) {
          for(j in 1:n_lv) {
            lambda[i, j] ~ dnorm(0, 0.1)
          }
        }"
    }

    model_str <- paste0("
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
        }", lambda_upper_str, "
        for(i in 1:n_lv) {
          lambda[i, i] ~ dnorm(0, 0.1) T(0,)
        }", lambda_lower_str, lambda_free_str, "
        for(j in 1:n_sp) {
          for(m in 1:n_covars) {
            beta[j, m] ~ dnorm(0, 0.1)
          }
        }
      }")
  }
  else {
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

  time_max <- 3600 * mcmc_params$time_max_hours

  mcmc_control <- list(time.max                   = time_max,
                       round.thinmult             = TRUE,
                       print.diagnostics          = FALSE,
                       Ncycles.target             = 2,
                       check.convergence.firstrun = FALSE,
                       convtype                   = "Gelman",
                       seed                       = seed)

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
  coef_rows   <- grep("^beta0\\[|^beta\\[", conv_names)

  if (length(coef_rows) > 0) {
    result$psrf_max <- max(conv_values[coef_rows, "Point est."], na.rm = TRUE)
  }

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
      result$scores <- apply(W_arr, c(2, 3), mean)
    }
  }

  for (species in 1:n_sp) {
    result$estimates$alpha[species]       <- mean(beta0_mat[, species], na.rm = TRUE)
    result$standard_errors$alpha[species] <- sd(beta0_mat[, species], na.rm = TRUE)
    result$estimates$beta[species]        <- mean(beta_arr[, species, 1], na.rm = TRUE)
    result$standard_errors$beta[species]  <- sd(beta_arr[, species, 1], na.rm = TRUE)

    if (fit_gamma) {
      result$estimates$gamma[species]       <- mean(beta_arr[, species, 2], na.rm = TRUE)
      result$standard_errors$gamma[species] <- sd(beta_arr[, species, 2], na.rm = TRUE)
    }
  }

  result$ess_alpha <- mean(neff_values[grep("^beta0\\[", neff_names)], na.rm = TRUE)
  result$ess_beta  <- mean(neff_values[grep("^beta\\[[0-9]+,1\\]", neff_names)], na.rm = TRUE)

  if (fit_gamma) {
    result$ess_gamma <- mean(neff_values[grep("^beta\\[[0-9]+,2\\]", neff_names)], na.rm = TRUE)
  }

  result$converged <- attrs$final.params$converged

  unlink(model_file)

  return(result)
}

################################################################################
