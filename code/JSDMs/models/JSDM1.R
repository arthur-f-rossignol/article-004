################################################################################
##                                                                            ##
##                          JSDM.1 FITTED WITH gllvm                          ##
##                           (JSDMs/models/JSDM1.R)                           ##
##                                                                            ##
##       Addressing Missing Covariates in Species Distribution Models:        ##
##  Inferential Impacts and Mitigation via Joint Species Distribution Models  ##
##                                                                            ##
##                  Arthur F. Rossignol & Frédéric Gosselin                   ##
##                                                                            ##
##                                    2026                                    ##
##                                                                            ##
################################################################################

fit_JSDM.1 <- function(Y,
                       covariates,
                       n_lv,
                       mcmc_params,
                       seed) {

  set.seed(seed)

  n_sp <- ncol(Y)

  fit_gamma <- (length(covariates) >= 2)

  result <- results_initialization(n_sp, n_lv)

  covariate_names <- paste(names(covariates), collapse = " + ")

  X_data  <- as.data.frame(covariates)
  formula <- as.formula(paste("~ 1 +", covariate_names))

  control       <- list(optimizer = "nlminb",
                        max.iter  = mcmc_params$max_iter)

  control_start <- list(starting.val = "res",
                        n.init       = mcmc_params$n_init,
                        n.init.max   = mcmc_params$n_init_max)

  start_time <- Sys.time()

  fit <- gllvm::gllvm(y             = Y,
                      X             = X_data,
                      family        = binomial(link = "probit"),
                      method        = "LA",
                      num.lv        = n_lv,
                      formula       = formula,
                      seed          = seed,
                      trace         = FALSE,
                      control       = control,
                      control.start = control_start)

  result$computation_time <- as.numeric(difftime(Sys.time(),
                                                 start_time,
                                                 units = "secs"))

  Xcoef_estimates       <- fit$params$Xcoef
  Xcoef_standard_errors <- fit$sd$Xcoef

  result$estimates$alpha       <- as.vector(fit$params$beta0)
  result$standard_errors$alpha <- as.vector(fit$sd$beta0)
  result$estimates$beta        <- as.vector(Xcoef_estimates[, 1])
  result$standard_errors$beta  <- as.vector(Xcoef_standard_errors[, 1])

  if (fit_gamma) {
    result$estimates$gamma       <- as.vector(Xcoef_estimates[, 2])
    result$standard_errors$gamma <- as.vector(Xcoef_standard_errors[, 2])
  }

  if (fit$num.lv > 0) {
    sigma_lv <- as.numeric(fit$params$sigma.lv)
    loadings <- as.matrix(fit$params$theta)

    for (k in 1:fit$num.lv) {
      loadings[, k] <- loadings[, k] * sigma_lv[k]
    }

    result$scores   <- as.matrix(fit$lvs)
    result$loadings <- loadings
    result$lv_scale <- sigma_lv
  }

  result$converged <- fit$convergence

  return(result)
}

################################################################################
