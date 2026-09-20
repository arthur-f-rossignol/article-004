################################################################################
##                                                                            ##
##          JSDMs WITH OPTIONAL MISSING COVARIATE AND LATENT VARIABLES        ##
##                               (JSDMs/JSDM.R)                               ##
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

n_lv              <- NULL
missing_covariate <- NULL

if (length(ARGS) >= 5) {
  n_lv <- as.integer(ARGS[5])
}

if (length(ARGS) >= 6) {
  missing_covariate <- as.logical(ARGS[6])
}

## MAIN PARAMETERS #############################################################

DEFAULT_MCMC <- list(n_init         = 50,
                     n_init_max     = 1e4,
                     max_iter       = 1e6,
                     burnin         = 100000,
                     n_sample       = 100000,
                     thin           = 5,
                     samples        = 50000,
                     transient      = 50000,
                     n_iter_min     = 25000,
                     n_burnin_min   = 25000,
                     thin_min       = 10,
                     n_chains       = 3,
                     conv_max       = 1.05,
                     n_eff_min      = 1000,
                     time_max_hours = 24)

DEFAULT_PARAMS <- list(n_obs             = 1000,     # number of sites
                       n_sp              = 10,       # number of species
                       n_covars          = 2,        # number of covariates (intercept excluded)
                       n_lv              = 1,        # number of latent variables
                       missing_covariate = TRUE,     # X2 is missing (unobserved)
                       mcmc              = DEFAULT_MCMC)

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
                          params) {

  set.seed(seed)

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

## MODELS ######################################################################

script_args <- commandArgs(trailingOnly = FALSE)
file_arg    <- grep("^--file=", script_args, value = TRUE)
script_dir  <- getwd()

if (length(file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
}

models_dir <- script_dir

## RESULT TEMPLATE INITIALIZATION ##############################################

results_initialization <- function(n_sp,
                                   n_lv) {

  species <- paste0("sp", 1:n_sp)

  estimates <- data.frame(species = species,
                          alpha   = rep(NA, n_sp),
                          beta    = rep(NA, n_sp),
                          gamma   = rep(NA, n_sp))

  standard_errors <- data.frame(species = species,
                                alpha   = rep(NA, n_sp),
                                beta    = rep(NA, n_sp),
                                gamma   = rep(NA, n_sp))

  if (n_lv > 0) {
    lv_scale    <- rep(NA, n_lv)
    lv_scale_se <- rep(NA, n_lv)
  }
  else {
    lv_scale    <- NULL
    lv_scale_se <- NULL
  }

  result <- list(estimates        = estimates,
                 standard_errors  = standard_errors,
                 scores           = NULL,
                 loadings         = NULL,
                 lv_scale         = lv_scale,
                 lv_scale_se      = lv_scale_se,
                 converged        = NA,
                 ess_alpha        = NA,
                 ess_beta         = NA,
                 ess_gamma        = NA,
                 psrf_max         = NA,
                 computation_time = NA)

  return(result)
}

## RESULT ASSEMBLY #############################################################

assemble_results <- function(model_results,
                             truth,
                             run,
                             settings) {

  model_names <- names(model_results)

  coefficients <- data.frame()
  diagnostics  <- data.frame()
  latent       <- list()

  for (model_name in model_names) {

    model <- model_results[[model_name]]

    for (parameter in c("alpha", "beta", "gamma")) {
      block <- data.frame(model     = model_name,
                          species   = model$estimates$species,
                          parameter = parameter,
                          truth     = truth[[parameter]],
                          estimate  = model$estimates[[parameter]],
                          std_error = model$standard_errors[[parameter]],
                          converged = model$converged)

      coefficients <- rbind(coefficients, block)
    }

    diagnostics <- rbind(diagnostics,
                         data.frame(model            = model_name,
                                    converged        = model$converged,
                                    ess_alpha        = model$ess_alpha,
                                    ess_beta         = model$ess_beta,
                                    ess_gamma        = model$ess_gamma,
                                    psrf_max         = model$psrf_max,
                                    computation_time = model$computation_time))

    latent[[model_name]] <- list(scores      = model$scores,
                                 loadings    = model$loadings,
                                 lv_scale    = model$lv_scale,
                                 lv_scale_se = model$lv_scale_se)
  }

  rownames(coefficients) <- NULL
  rownames(diagnostics)  <- NULL

  results <- list(run          = run,
                  settings     = settings,
                  truth        = truth,
                  coefficients = coefficients,
                  diagnostics  = diagnostics,
                  latent       = latent)

  return(results)
}

## MODEL FITTING FUNCTIONS #####################################################

source(file.path(models_dir, "models", "JSDM1.R"))
source(file.path(models_dir, "models", "JSDM2.R"))
source(file.path(models_dir, "models", "JSDM3.R"))
source(file.path(models_dir, "models", "JSDM4.R"))
source(file.path(models_dir, "models", "JSDM5.R"))
source(file.path(models_dir, "models", "JSDM6.R"))

## DATA ANALYSIS ###############################################################

run_JSDM <- function(seed,
                     results_dir,
                     scenario_id,
                     methods,
                     params) {

  scenarios <- build_scenarios(params$n_sp)
  scenario  <- scenarios[[scenario_id]]

  data <- simulate_data(seed     = seed,
                        scenario = scenario,
                        params   = params)

  covariates <- select_covariates(data, params$missing_covariate)

  truth <- data.frame(species = paste0("sp", 1:params$n_sp),
                      alpha   = as.vector(data$B["alpha", ]),
                      beta    = as.vector(data$B["beta", ]),
                      gamma   = as.vector(data$B["gamma", ]))

  fit_functions <- list(fit_JSDM.1,
                        fit_JSDM.2,
                        fit_JSDM.3,
                        fit_JSDM.4,
                        fit_JSDM.5,
                        fit_JSDM.6)

  model_results <- list()
  failures      <- list()

  total_start <- Sys.time()

  for (m in methods) {
    model_name <- paste0("JSDM.", m)

    fit_result <- try(fit_functions[[m]](Y           = data$Y,
                                         covariates  = covariates,
                                         n_lv        = params$n_lv,
                                         mcmc_params = params$mcmc,
                                         seed        = seed),
                      silent = TRUE)

    if (inherits(fit_result, "try-error")) {
      failures[[model_name]] <- as.character(fit_result)
      fit_result             <- results_initialization(params$n_sp, params$n_lv)
    }

    model_results[[model_name]] <- fit_result
  }

  total_time <- as.numeric(difftime(Sys.time(),
                                    total_start,
                                    units = "secs"))

  run <- list(seed                   = seed,
              timestamp              = Sys.time(),
              scenario_id            = scenario_id,
              scenario               = scenario$name,
              methods                = methods,
              failures               = failures,
              total_computation_time = total_time)

  settings <- list(n_obs             = params$n_obs,
                   n_sp              = params$n_sp,
                   n_covars          = params$n_covars,
                   n_lv              = params$n_lv,
                   missing_covariate = params$missing_covariate,
                   mcmc              = params$mcmc)

  results <- assemble_results(model_results = model_results,
                              truth         = truth,
                              run           = run,
                              settings      = settings)

  if (params$missing_covariate) {
    mc_tag <- "MC"
  }
  else {
    mc_tag <- "noMC"
  }

  lv_tag <- sprintf("%dLV", params$n_lv)

  if (length(methods) == 1) {
    output_file <- file.path(results_dir,
                             sprintf("JSDM%d_S%d_%s_%s_seed%06d.RData",
                                     methods, scenario_id, mc_tag, lv_tag, seed))
  }
  else {
    method_str  <- paste(methods, collapse = "_")
    output_file <- file.path(results_dir,
                             sprintf("JSDM_S%d_%s_%s_seed%06d_methods%s.RData",
                                     scenario_id, mc_tag, lv_tag, seed,
                                     method_str))
  }

  save(results, file = output_file)

  return(output_file)
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
         results_dir = results_dir,
         scenario_id = scenario_id,
         methods     = methods,
         params      = params)

################################################################################