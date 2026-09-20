################################################################################
##                                                                            ##
##   JSDMs WITH VARYING NUMBER OF SPECIES RESPONDING TO A MISSING COVARIATE   ##
##                (JSDMs/JSDM_number_species_responding_MC.R)                 ##
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

task_id <- as.integer(ARGS[1])
results_dir <- ARGS[2]

## PARAMETERS ##################################################################

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

DEFAULT_PARAMS <- list(n_obs      = 1000,
                       n_sp       = 10,
                       n_lv       = 1,
                       true_alpha = rep(1, 10),
                       true_beta  = rep(1, 10),
                       true_gamma = rep(1, 10),
                       mcmc       = DEFAULT_MCMC)

## DATA GENERATING PROCESS #####################################################

simulate_data <- function(replicate_id,
                          k_responding,
                          params) {
  
  seed <- replicate_id
  set.seed(seed)
  
  n_obs <- params$n_obs
  n_sp  <- params$n_sp
  
  X1 <- rnorm(n_obs, mean = 0, sd = 1)
  X2 <- rnorm(n_obs, mean = 0, sd = 1)
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
  Y <- matrix(NA, 
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

model_result <- fit_functions[[model_id]](Y           = data$Y,
                                          covariates  = list(X1 = data$X1),
                                          n_lv        = data$params$n_lv,
                                          mcmc_params = data$params$mcmc,
                                          seed        = data$seed)

model_results               <- list()
model_results[[model_name]] <- model_result

truth <- data.frame(species = paste0("sp", 1:DEFAULT_PARAMS$n_sp),
                    alpha   = as.vector(data$B["alpha", ]),
                    beta    = as.vector(data$B["beta", ]),
                    gamma   = as.vector(data$B["gamma", ]))

run <- list(task_id      = task_id,
            model_id     = model_id,
            model_name   = model_name,
            replicate_id = replicate_id,
            k_responding = k_responding,
            seed         = data$seed,
            timestamp    = Sys.time())

settings <- list(n_obs = DEFAULT_PARAMS$n_obs,
                 n_sp  = DEFAULT_PARAMS$n_sp,
                 n_lv  = DEFAULT_PARAMS$n_lv,
                 mcmc  = DEFAULT_PARAMS$mcmc)

results <- assemble_results(model_results = model_results,
                            truth         = truth,
                            run           = run,
                            settings      = settings)

file_name <- file.path(results_dir, sprintf("%s_k%02d_rep%03d.RData",
                                            gsub("\\.", "", model_name),
                                            k_responding,
                                            replicate_id))

save(results, file = file_name)

################################################################################
