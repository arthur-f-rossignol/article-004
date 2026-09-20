################################################################################
##                                                                            ##
##                      SDMs WITH ONE MISSING COVARIATE                       ##
##                                (SDMs/SDM.R)                                ##
##                                                                            ##
##       Addressing Missing Covariates in Species Distribution Models:        ##
##  Inferential Impacts and Mitigation via Joint Species Distribution Models  ##
##                                                                            ##
##                  Arthur F. Rossignol & Frederic Gosselin                   ##
##                                                                            ##
##                                    2026                                    ##
##                                                                            ##
################################################################################

## PACKAGES ####################################################################

library(Matrix)
library(lme4)
library(robustbase)
library(glmmTMB)
library(nimble)
library(runjags)
library(R2jags)
library(coda)
library(runMCMCbtadjust)
library(GLMMadaptive)
library(mgcv)
library(spaMM)
library(cmdstanr)
library(brms)
library(INLA)

## ARGUMENTS FROM SLURM ########################################################

ARGS <- commandArgs(trailingOnly = TRUE)

seed <- as.integer(ARGS[1])
results_dir <- ARGS[2]

## PARAMETERIZATION ############################################################

n_obs       <- 1000       # number of sites
sigma_error <- 1          # standard deviation for error terms in Gaussian model

alpha_true <- 1           # intercept
beta_true  <- 1           # coefficient for covariate X1 (observed)
gamma_true <- 1           # coefficient for covariate X2 (unobserved / missing)

## LIST OF MODELS ##############################################################

data_models <- c("gaussian",               # Gaussian with identity link
                 "poisson",                # Poisson with log link
                 "bernoulli_probit",       # Bernoulli with probit link
                 "bernoulli_logit",        # Bernoulli with logit link
                 "bernoulli_cloglog")      # Bernoulli with cloglog link

fit_methods <- c("SDM.1_no_OLRE",
                 "SDM.1_OLRE",
                 "SDM.2_no_OLRE",
                 "SDM.2_OLRE",
                 "SDM.3_no_OLRE",
                 "SDM.3_OLRE",
                 "SDM.4",
                 "SDM.5_no_OLRE",
                 "SDM.5_OLRE",
                 "SDM.6_no_OLRE",
                 "SDM.6_OLRE",
                 "SDM.7_no_OLRE",
                 "SDM.7_OLRE",
                 "SDM.8_no_OLRE",
                 "SDM.8_OLRE",
                 "SDM.9_no_OLRE",
                 "SDM.10_no_OLRE",
                 "SDM.10_OLRE")

frequentist_OLRE_only_methods <- c("SDM.4", "SDM.9")

method_to_sdm <- c(SDM.1  = "SDM.1",
                   SDM.2  = "SDM.2",
                   SDM.3  = "SDM.3",
                   SDM.4  = "SDM.4",
                   SDM.5  = "SDM.5",
                   SDM.6  = "SDM.6",
                   SDM.7  = "SDM.7",
                   SDM.8  = "SDM.8",
                   SDM.9  = "SDM.9",
                   SDM.10 = "SDM.10")

method_compatibility <- list(SDM.1  = c("gaussian", 
                                        "poisson",
                                        "bernoulli_probit", 
                                        "bernoulli_logit",
                                        "bernoulli_cloglog"),
                             SDM.2  = c("gaussian", 
                                        "poisson",
                                        "bernoulli_probit",
                                        "bernoulli_logit",
                                        "bernoulli_cloglog"),
                             SDM.3  = c("gaussian", 
                                        "poisson",
                                        "bernoulli_probit", 
                                        "bernoulli_logit",
                                        "bernoulli_cloglog"),
                             SDM.4  = c("poisson",
                                        "bernoulli_probit", 
                                        "bernoulli_logit"),
                             SDM.5  = c("gaussian", 
                                        "poisson",
                                        "bernoulli_probit", 
                                        "bernoulli_logit"),
                             SDM.6  = c("gaussian", 
                                        "poisson",
                                        "bernoulli_probit", 
                                        "bernoulli_logit",
                                        "bernoulli_cloglog"),
                             SDM.7  = c("gaussian", 
                                        "poisson",
                                        "bernoulli_probit", 
                                        "bernoulli_logit",
                                        "bernoulli_cloglog"),
                             SDM.8  = c("gaussian", 
                                        "poisson",
                                        "bernoulli_probit", 
                                        "bernoulli_logit",
                                        "bernoulli_cloglog"),
                             SDM.9  = c("gaussian", 
                                        "poisson",
                                        "bernoulli_probit", 
                                        "bernoulli_logit",
                                        "bernoulli_cloglog"),
                             SDM.10 = c("gaussian",
                                        "poisson",
                                        "bernoulli_logit"))

is_olre_variant <- function(method_name) {
  (grepl("_OLRE$", method_name) && !grepl("_no_OLRE$", method_name)) ||
    (method_name %in% frequentist_OLRE_only_methods)
}

is_method_compatible <- function(method_name, data_model) {
  base_method <- sub("_(no_)?OLRE$", "", method_name)
  if (data_model == "gaussian" && is_olre_variant(method_name)) {
    return(FALSE)
  }
  if (base_method %in% names(method_compatibility)) {
    return(data_model %in% method_compatibility[[base_method]])
  }
  return(FALSE)
}

## DATA GENERATION #############################################################

covariates_generation <- function(n, seed) {
  
  set.seed(seed)
  
  X1 <- rnorm(n, 0, 1)
  X2 <- rnorm(n, 0, 1)
  
  return(list(X1 = X1, X2 = X2))
}

data_generation <- function(n, X1, X2, model_type, seed) {
  
  set.seed(seed)
  
  eta <- alpha_true + beta_true * X1 + gamma_true * X2
  
  Y <- switch(model_type,
              "gaussian"          = rnorm(n, eta, sigma_error),
              "poisson"           = rpois(n, exp(eta)),
              "bernoulli_probit"  = rbinom(n, 1, pnorm(eta)),
              "bernoulli_logit"   = rbinom(n, 1, plogis(eta)),
              "bernoulli_cloglog" = rbinom(n, 1, 1 - exp(-exp(eta))))
  
  return(Y)
}

## MODELS ######################################################################

script_args <- commandArgs(trailingOnly = FALSE)
file_arg    <- grep("^--file=", script_args, value = TRUE)
script_dir  <- getwd()

if (length(file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
}

models_dir <- script_dir

## RESULT ASSEMBLY #############################################################

assemble_results <- function(fits,
                             truth,
                             run,
                             settings) {

  scalar <- function(value) {

    if (length(value) == 0) {
      return(NA)
    }

    return(unname(value[1]))
  }

  coefficients <- data.frame()
  diagnostics  <- data.frame()

  for (fit in fits) {

    for (parameter in c("alpha", "beta")) {
      block <- data.frame(data_model = scalar(fit$data_model),
                          method     = scalar(fit$method),
                          SDM        = scalar(fit$SDM_label),
                          OLRE       = scalar(fit$OLRE),
                          parameter  = parameter,
                          truth      = scalar(truth[[parameter]]),
                          estimate   = scalar(fit[[paste0(parameter, "_est")]]),
                          std_error  = scalar(fit[[paste0(parameter, "_SE")]]),
                          converged  = scalar(fit$converged))

      coefficients <- rbind(coefficients, block)
    }

    diagnostics <- rbind(diagnostics,
                         data.frame(data_model = scalar(fit$data_model),
                                    method     = scalar(fit$method),
                                    SDM        = scalar(fit$SDM_label),
                                    OLRE       = scalar(fit$OLRE),
                                    converged  = scalar(fit$converged),
                                    sigma2_est = scalar(fit$sigma2_est),
                                    fit_time   = scalar(fit$fit_time)))
  }

  rownames(coefficients) <- NULL
  rownames(diagnostics)  <- NULL

  results <- list(run          = run,
                  settings     = settings,
                  truth        = truth,
                  coefficients = coefficients,
                  diagnostics  = diagnostics)

  return(results)
}

## MODEL FITTING FUNCTIONS #####################################################

source(file.path(models_dir, "models", "SDM1.R"))
source(file.path(models_dir, "models", "SDM2.R"))
source(file.path(models_dir, "models", "SDM3.R"))
source(file.path(models_dir, "models", "SDM4.R"))
source(file.path(models_dir, "models", "SDM5.R"))
source(file.path(models_dir, "models", "SDM6.R"))
source(file.path(models_dir, "models", "SDM7.R"))
source(file.path(models_dir, "models", "SDM8.R"))
source(file.path(models_dir, "models", "SDM9.R"))
source(file.path(models_dir, "models", "SDM10.R"))

## DISPATCHER ##################################################################

fit_model <- function(Y, X1, model_type, method, seed) {
  
  result <- list(method     = method,
                 SDM_label  = NA,
                 OLRE       = NA,
                 alpha_est  = NA,
                 alpha_SE   = NA,
                 beta_est   = NA,
                 beta_SE    = NA,
                 sigma2_est = NA,
                 converged  = NA)
  
  if (method %in% frequentist_OLRE_only_methods) {
    use_OLRE    <- TRUE
    base_method <- method
  } else if (grepl("_no_OLRE$", method)) {
    use_OLRE    <- FALSE
    base_method <- sub("_no_OLRE$", "", method)
  } else if (grepl("_OLRE$", method)) {
    use_OLRE    <- TRUE
    base_method <- sub("_OLRE$", "", method)
  } else {
    use_OLRE    <- FALSE
    base_method <- method
  }

  # SDM.4 always fits an OLRE; SDM.9 never does — override the generic label
  if (base_method == "SDM.4") {
    use_OLRE <- TRUE
  }
  if (base_method == "SDM.9") {
    use_OLRE <- FALSE
  }

  sub_result <- switch(base_method,
                       "SDM.1"  = fit_SDM.1(Y, X1, model_type, use_OLRE),
                       "SDM.2"  = fit_SDM.2(Y, X1, model_type, use_OLRE),
                       "SDM.3"  = fit_SDM.3(Y, X1, model_type, use_OLRE),
                       "SDM.4"  = fit_SDM.4(Y, X1, model_type),
                       "SDM.5"  = fit_SDM.5(Y, X1, model_type, use_OLRE, seed),
                       "SDM.6"  = fit_SDM.6(Y, X1, model_type, use_OLRE),
                       "SDM.7"  = fit_SDM.7(Y, X1, model_type, use_OLRE),
                       "SDM.8"  = fit_SDM.8(Y, X1, model_type, use_OLRE),
                       "SDM.9"  = fit_SDM.9(Y, X1, model_type),
                       "SDM.10" = fit_SDM.10(Y, X1, model_type, use_OLRE))
  
  result           <- modifyList(result, sub_result)
  result$SDM_label <- unname(method_to_sdm[base_method])
  result$OLRE      <- use_OLRE
  
  return(result)
}

## SINGLE REPLICATE ############################################################

single_replicate <- function(seed) {

  set.seed(seed)

  fits <- list()

  covariates <- covariates_generation(n_obs, seed)
  X1         <- covariates$X1
  X2         <- covariates$X2

  for (data_model in data_models) {

    Y <- data_generation(n_obs,
                         X1,
                         X2,
                         data_model,
                         seed + which(data_models == data_model))

    for (method in fit_methods) {

      if (!is_method_compatible(method, data_model)) {
        next
      }

      start_time <- Sys.time()

      fit_result <- fit_model(Y, X1, data_model, method, seed)

      end_time <- Sys.time()

      fit_result$data_model <- data_model
      fit_result$fit_time   <- as.numeric(difftime(end_time,
                                                   start_time,
                                                   units = "secs"))

      fits[[length(fits) + 1]] <- fit_result
    }
  }

  return(fits)
}

## RUN #########################################################################

start_time <- Sys.time()

fits <- single_replicate(seed)

end_time <- Sys.time()

truth <- list(alpha       = alpha_true,
              beta        = beta_true,
              gamma       = gamma_true,
              sigma_error = sigma_error)

run <- list(seed       = seed,
            timestamp  = Sys.time(),
            total_time = as.numeric(difftime(end_time,
                                             start_time,
                                             units = "secs")))

settings <- list(n_obs       = n_obs,
                 sigma_error = sigma_error)

results <- assemble_results(fits     = fits,
                            truth    = truth,
                            run      = run,
                            settings = settings)

output_file <- file.path(results_dir, sprintf("replicate_%06d.RData", seed))

save(results, file = output_file)

################################################################################