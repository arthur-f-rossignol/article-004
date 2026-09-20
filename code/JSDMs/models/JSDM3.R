################################################################################
##                                                                            ##
##                          JSDM.3 FITTED WITH Hmsc                           ##
##                           (JSDMs/models/JSDM3.R)                           ##
##                                                                            ##
##       Addressing Missing Covariates in Species Distribution Models:        ##
##  Inferential Impacts and Mitigation via Joint Species Distribution Models  ##
##                                                                            ##
##                  Arthur F. Rossignol & Frédéric Gosselin                   ##
##                                                                            ##
##                                    2026                                    ##
##                                                                            ##
################################################################################

fit_JSDM.3 <- function(Y,
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

  X_data   <- as.data.frame(covariates)
  XFormula <- as.formula(paste("~ 1 +", covariate_names))

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
  }
  else {
    model <- Hmsc::Hmsc(Y        = Y,
                        XData    = X_data,
                        XFormula = XFormula,
                        distr    = "probit")
  }

  start_time <- Sys.time()

  fit <- Hmsc::sampleMcmc(model,
                          thin      = mcmc_params$thin,
                          samples   = mcmc_params$samples,
                          transient = mcmc_params$transient,
                          nChains   = mcmc_params$n_chains,
                          nParallel = mcmc_params$n_chains,
                          verbose   = 0,
                          initPar   = "fixed effects")

  result$computation_time <- as.numeric(difftime(Sys.time(),
                                                 start_time,
                                                 units = "secs"))

  postBeta <- Hmsc::getPostEstimate(fit, parName = "Beta")

  result$estimates$alpha <- postBeta$mean[1, ]
  result$estimates$beta  <- postBeta$mean[2, ]

  if (fit_gamma) {
    result$estimates$gamma <- postBeta$mean[3, ]
  }

  mpost <- Hmsc::convertToCodaObject(fit)

  beta_samples <- do.call(rbind, mpost$Beta)

  n_fixed <- length(covariates) + 1

  alpha_id <- seq(1, n_sp * n_fixed, by = n_fixed)
  beta_id  <- seq(2, n_sp * n_fixed, by = n_fixed)
  gamma_id <- seq(3, n_sp * n_fixed, by = n_fixed)

  alpha_samples <- beta_samples[, alpha_id, drop = FALSE]
  beta1_samples <- beta_samples[, beta_id, drop = FALSE]

  result$standard_errors$alpha <- apply(alpha_samples, 2, sd)
  result$standard_errors$beta  <- apply(beta1_samples, 2, sd)

  if (fit_gamma) {
    gamma_samples <- beta_samples[, gamma_id, drop = FALSE]
    result$standard_errors$gamma <- apply(gamma_samples, 2, sd)
  }

  if (n_lv > 0) {
    postEta    <- Hmsc::getPostEstimate(fit, parName = "Eta")
    postLambda <- Hmsc::getPostEstimate(fit, parName = "Lambda")

    if (!is.null(postEta)) {
      result$scores <- postEta$mean[[1]]
    }

    if (!is.null(postLambda)) {
      result$loadings <- t(postLambda$mean[[1]])
    }
  }

  ess_values <- effectiveSize(mpost$Beta)

  result$ess_alpha <- mean(ess_values[alpha_id], na.rm = TRUE)
  result$ess_beta  <- mean(ess_values[beta_id], na.rm = TRUE)

  if (fit_gamma) {
    result$ess_gamma <- mean(ess_values[gamma_id], na.rm = TRUE)
  }

  if (length(mpost$Beta) > 1) {
    gelman_diagnostic <- gelman.diag(mpost$Beta,
                                     multivariate = FALSE,
                                     autoburnin   = FALSE)
    psrf_values       <- gelman_diagnostic$psrf[, "Point est."]

    result$psrf_max  <- max(psrf_values, na.rm = TRUE)
    result$converged <- is.finite(result$psrf_max) &&
                        result$psrf_max < mcmc_params$conv_max
  }

  return(result)
}

################################################################################
