################################################################################
##                                                                            ##
##                          JSDM.2 FITTED WITH jSDM                           ##
##                           (JSDMs/models/JSDM2.R)                           ##
##                                                                            ##
##       Addressing Missing Covariates in Species Distribution Models:        ##
##  Inferential Impacts and Mitigation via Joint Species Distribution Models  ##
##                                                                            ##
##                  Arthur F. Rossignol & Frédéric Gosselin                   ##
##                                                                            ##
##                                    2026                                    ##
##                                                                            ##
################################################################################

fit_JSDM.2 <- function(Y,
                       covariates,
                       n_lv,
                       mcmc_params,
                       seed) {

  set.seed(seed)

  n_sp  <- ncol(Y)
  n_obs <- nrow(Y)

  fit_gamma <- (length(covariates) >= 2)

  result <- results_initialization(n_sp, n_lv)

  covariate_names <- paste(names(covariates), collapse = " + ")

  X_data       <- as.data.frame(covariates)
  site_formula <- as.formula(paste("~", covariate_names))

  n_chains <- mcmc_params$n_chains

  # overdispersed starting values and distinct seed for each chain

  beta_start_chain   <- seq(-0.5, 0.5, length.out = n_chains)
  lambda_start_chain <- seq(0.5, 1.5, length.out = n_chains)
  seed_chain         <- seed + (1:n_chains - 1) * 1000000

  fit_chain <- function(chain) {

    fit <- jSDM::jSDM_binomial_probit(burnin        = mcmc_params$burnin,
                                      mcmc          = mcmc_params$n_sample,
                                      thin          = mcmc_params$thin,
                                      presence_data = Y,
                                      site_formula  = site_formula,
                                      site_data     = X_data,
                                      n_latent      = n_lv,
                                      site_effect   = "none",
                                      beta_start    = beta_start_chain[chain],
                                      lambda_start  = lambda_start_chain[chain],
                                      W_start       = 0,
                                      alpha_start   = 0,
                                      V_alpha       = 1,
                                      mu_beta       = 0,
                                      V_beta        = 1,
                                      mu_lambda     = 0,
                                      V_lambda      = 1,
                                      shape_Valpha  = 0.5,
                                      rate_Valpha   = 0.0005,
                                      seed          = seed_chain[chain],
                                      verbose       = 0)

    return(fit)
  }

  start_time <- Sys.time()

  fits <- parallel::mclapply(1:n_chains,
                             fit_chain,
                             mc.cores       = n_chains,
                             mc.preschedule = FALSE)

  result$computation_time <- as.numeric(difftime(Sys.time(),
                                                 start_time,
                                                 units = "secs"))

  chain_failed <- sapply(fits, function(fit) inherits(fit, "try-error"))

  if (any(chain_failed)) {
    stop("Error: at least one jSDM chain failed.\n")
  }

  n_covars <- length(covariates) + 1

  mean_beta <- matrix(NA, nrow = n_sp, ncol = n_covars)
  se_beta   <- matrix(NA, nrow = n_sp, ncol = n_covars)

  if (n_lv > 0) {
    mean_lambda <- matrix(NA, nrow = n_sp, ncol = n_lv)
  }

  ess_alpha <- numeric(n_sp)
  ess_beta  <- numeric(n_sp)
  ess_gamma <- numeric(n_sp)
  psrf_sp   <- numeric(n_sp)

  for (j in 1:n_sp) {

    # chains of species j gathered in a single coda object

    sp_chains <- coda::mcmc.list(lapply(fits, function(fit) fit$mcmc.sp[[j]]))

    # posterior summaries computed on the pooled draws

    sp_samples <- as.matrix(sp_chains)

    for (covar_id in 1:n_covars) {
      mean_beta[j, covar_id] <- mean(sp_samples[, covar_id])
      se_beta[j, covar_id]   <- sd(sp_samples[, covar_id])
    }

    if (n_lv > 0) {
      for (l in 1:n_lv) {
        mean_lambda[j, l] <- mean(sp_samples[, n_covars + l])
      }
    }

    ess_all <- coda::effectiveSize(sp_chains)

    if (length(ess_all) >= 1) {
      ess_alpha[j] <- ess_all[1]
    }
    if (length(ess_all) >= 2) {
      ess_beta[j] <- ess_all[2]
    }
    if (fit_gamma && length(ess_all) >= 3) {
      ess_gamma[j] <- ess_all[3]
    }

    # Gelman-Rubin diagnostic on the species regression coefficients

    gelman_diagnostic <- coda::gelman.diag(sp_chains[, 1:n_covars],
                                           multivariate = FALSE,
                                           autoburnin   = FALSE)

    psrf_sp[j] <- max(gelman_diagnostic$psrf[, "Point est."], na.rm = TRUE)
  }

  result$estimates$alpha       <- mean_beta[, 1]
  result$standard_errors$alpha <- se_beta[, 1]
  result$estimates$beta        <- mean_beta[, 2]
  result$standard_errors$beta  <- se_beta[, 2]

  if (fit_gamma) {
    result$estimates$gamma       <- mean_beta[, 3]
    result$standard_errors$gamma <- se_beta[, 3]
  }

  if (n_lv > 0) {
    result$loadings <- mean_lambda

    latent_names <- names(fits[[1]]$mcmc.latent)

    if (!is.null(latent_names) && length(latent_names) > 0) {
      W_mean <- matrix(NA, nrow = n_obs, ncol = n_lv)

      for (l in 1:n_lv) {
        lv_name <- paste0("lv_", l)

        if (lv_name %in% latent_names) {
          lv_draws <- lapply(fits, function(fit) fit$mcmc.latent[[lv_name]])

          lv_chains   <- coda::mcmc.list(lv_draws)
          lv_samples  <- as.matrix(lv_chains)
          W_mean[, l] <- colMeans(lv_samples)
        }
      }

      result$scores <- W_mean
    }
  }

  result$ess_alpha <- mean(ess_alpha, na.rm = TRUE)
  result$ess_beta  <- mean(ess_beta, na.rm = TRUE)

  if (fit_gamma) {
    result$ess_gamma <- mean(ess_gamma, na.rm = TRUE)
  }

  result$psrf_max  <- max(psrf_sp, na.rm = TRUE)
  result$converged <- is.finite(result$psrf_max) &&
                      result$psrf_max < mcmc_params$conv_max

  return(result)
}

################################################################################
