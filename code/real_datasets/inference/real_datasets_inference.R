################################################################################
##                                                                            ##
##                    REAL ECOLOGICAL DATASETS • INFERENCE                    ##
##            (real_datasets/inference/real_datasets_inference.R)             ##
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

library(tidyverse)
library(jSDM)
library(vegan)
library(coda)
library(parallel)
library(readxl)

## ARGUMENTS FROM SLURM ########################################################

ARGS <- commandArgs(trailingOnly = TRUE)

dataset_id <- as.integer(ARGS[1])

base_dir    <- ARGS[2]
data_dir    <- file.path(base_dir, "datasets")
results_dir <- file.path(base_dir, "inference", "results")

## CONFIGURATION ###############################################################

DATASETS <- c("dataset.1",
              "dataset.2",
              "dataset.3",
              "dataset.4",
              "dataset.5",
              "dataset.6")

CONFIGS <- list(dataset.1 = list(n_pcs   = 5,
                                 min_occ = 0.05,
                                 max_occ = 0.95),
                dataset.2 = list(n_pcs   = 4,
                                 min_occ = 0.05,
                                 max_occ = 0.95),
                dataset.3 = list(n_pcs   = 5,
                                 min_occ = 0.05,
                                 max_occ = 0.95),
                dataset.4 = list(n_pcs   = 5,
                                 min_occ = 0.05,
                                 max_occ = 0.95),
                dataset.5 = list(n_pcs   = 3,
                                 min_occ = 0.05,
                                 max_occ = 0.95),
                dataset.6 = list(n_pcs   = 2,
                                 min_occ = 0,
                                 max_occ = 1))

N_LV <- 1

DEFAULT_PARAMS <- list(mcmc = list(n_iter   = 100000,
                                   n_burnin = 50000,
                                   thin     = 10,
                                   n_chains = 3),
                       prior = list(mu_beta   = 0,
                                    V_beta    = 10,
                                    mu_lambda = 0,
                                    V_lambda  = 10))

dataset_name <- DATASETS[dataset_id]
cfg          <- CONFIGS[[dataset_name]]

seed <- dataset_id

## HELPER FUNCTIONS ############################################################

species_filtering <- function(Y,
                              min_occ,
                              max_occ) {
  n     <- nrow(Y)
  occ   <- colSums(Y)
  min_c <- if (min_occ < 1) {
    ceiling(n * min_occ)
  } else {
    min_occ
  }
  max_c <- if (max_occ <= 1) {
    floor(n * max_occ)
  } else {
    max_occ
  }
  keep <- which(occ >= min_c & occ <= max_c)
  return(Y[, keep, drop = FALSE])
}

PCA <- function(X,
                n_pcs) {
  X    <- as.matrix(X)
  sds  <- apply(X, 2, sd)
  keep <- which(is.finite(sds) & sds > 0)
  if (length(keep) < ncol(X)) {
    X <- X[, keep, drop = FALSE]
  }
  pca   <- prcomp(scale(X), center = FALSE, scale. = FALSE)
  n_use <- min(n_pcs, ncol(X), nrow(X) - 1)
  pc    <- pca$x[, 1:n_use, drop = FALSE]
  colnames(pc) <- paste0("PC", 1:n_use)
  return(list(scores = pc, pc1 = pc[, 1], n_pcs = n_use, pca = pca))
}

## DATASET LOADING #############################################################

data_loading <- function(name) {
  
  if (name == "dataset.1") {
    
    spe <- read.csv(file.path(data_dir, "birds", "Birds_PA.csv"))
    env <- read.csv(file.path(data_dir, "birds", "Birds_Cov.csv"))
    
    Y <- ifelse(as.matrix(spe) > 0, 1, 0)
    X <- data.matrix(env)
    
  } else if (name == "dataset.2") {
    
    spe <- read.csv(file.path(data_dir, "butterflies", "Butterfly_PA.csv"))
    env <- read.csv(file.path(data_dir, "butterflies", "Butterfly_Cov.csv"))
    
    Y <- ifelse(as.matrix(spe) > 0, 1, 0)
    X <- data.matrix(env[, -1])
    
  } else if (name == "dataset.3") {
    
    data <- read_rds(file.path(data_dir, "mara", "mara_animal_compiled.rds")) %>%
      mutate(Site      = as.numeric(Site),
             sin_month = sin(2 * pi * Month / 12),
             cos_month = cos(2 * pi * Month / 12)) %>%
      filter(Yr_Mo != "2018-05",
             !is.na(Protein_lag1),
             !is.na(Height_lag1)) %>%
      dplyr::select(month_id, Name, x, y,
                    Cattle, Wildebeest, Zebra, Thompsons_Gazelle, Impala, Topi,
                    Eland, Buffalo, Grants_Gazelle, Waterbuck, Dikdik, Elephant,
                    Site, Pgrazed_lag1, Precip, Protein_lag1, Height_lag1,
                    sin_month, cos_month) %>%
      drop_na(.)
    
    Y <- data %>%
      dplyr::select(Cattle, Wildebeest, Zebra, Thompsons_Gazelle, Impala, Topi,
                    Eland, Buffalo, Grants_Gazelle, Waterbuck, Dikdik, Elephant) %>%
      as.matrix(.)
    Y <- ifelse(Y > 0, 1, 0)
    
    XData <- data %>%
      dplyr::select(Pgrazed_lag1, Precip, Protein_lag1, Height_lag1,
                    sin_month, cos_month)
    X <- data.matrix(XData)
    
  } else if (name == "dataset.4") {
    
    data("eucalypts", package = "jSDM")
    
    Y <- ifelse(as.matrix(eucalypts[, 1:12]) > 0, 1, 0)
    X <- data.matrix(eucalypts[, 13:19])
    
  } else if (name == "dataset.5") {
    
    data <- read.csv(file.path(data_dir, "kilpisjarvi",
                               "Kilpisjarvi_plant_data.csv"))
    
    sites <- unique(data$site)
    ny    <- length(sites)
    sp    <- unique(data$species)
    ns    <- length(sp)
    Y.pa  <- matrix(0, nrow = ny, ncol = ns)
    T3_GDD3 <- rep(NA, ny)
    T3_FDD  <- rep(NA, ny)
    moist_mean_summer <- rep(NA, ny)
    
    for (k in 1:nrow(data)) {
      i <- which(data$site[k] == sites)
      j <- which(data$species[k] == sp)
      Y.pa[i, j]           <- 1
      T3_GDD3[i]           <- data$T3_GDD3[k]
      moist_mean_summer[i] <- data$moist_mean_summer[k]
      T3_FDD[i]            <- data$T3_FDD[k]
    }
    
    sp.order    <- order(colSums(Y.pa), decreasing = TRUE)
    Y           <- Y.pa[, sp.order]
    colnames(Y) <- sp[sp.order]
    
    X <- data.frame(GDD = T3_GDD3, FDD = T3_FDD, SM = moist_mean_summer)
    X <- data.matrix(X)
    
  } else if (name == "dataset.6") {
    
    data("aravo", package = "jSDM")
    
    Y <- apply(aravo$spe > 0, 2, as.numeric)
    Y <- Y[, colSums(Y) >= 5]
    X <- data.matrix(cbind(poly(aravo$env$Snow, 2)))
    
  }
  
  return(list(Y = Y, X = X))
}

## jSDM MCMC ENGINE ############################################################

run_one_chain <- function(Y,
                          X_data,
                          site_formula,
                          n_lv,
                          n_burnin,
                          n_sample,
                          thin,
                          seed,
                          beta_start) {
  
  jSDM_args <- list(burnin        = n_burnin,
                    mcmc          = n_sample,
                    thin          = thin,
                    presence_data = Y,
                    site_formula  = site_formula,
                    site_data     = X_data,
                    n_latent      = n_lv,
                    site_effect   = "none",
                    beta_start    = beta_start,
                    mu_beta       = DEFAULT_PARAMS$prior$mu_beta,
                    V_beta        = DEFAULT_PARAMS$prior$V_beta,
                    seed          = seed,
                    verbose       = 0)
  
  if (n_lv > 0) {
    jSDM_args$lambda_start <- 0
    jSDM_args$W_start      <- 0
    jSDM_args$mu_lambda    <- DEFAULT_PARAMS$prior$mu_lambda
    jSDM_args$V_lambda     <- DEFAULT_PARAMS$prior$V_lambda
  }
  
  fit <- do.call(jSDM::jSDM_binomial_probit, jSDM_args)
  
  slim <- list(mcmc.sp       = fit$mcmc.sp,
               mcmc.latent   = fit$mcmc.latent,
               mcmc.Deviance = as.numeric(fit$mcmc.Deviance))
  rm(fit)
  gc(verbose = FALSE)
  
  return(slim)
}

## JSDM FITTING FUNCTION #######################################################

fit_JSDM <- function(Y,
                     X,
                     cfg,
                     seed) {
  
  Y <- as.matrix(Y)
  if (is.null(colnames(Y))) {
    colnames(Y) <- paste0("sp", seq_len(ncol(Y)))
  }
  colnames(Y) <- make.names(colnames(Y), unique = TRUE)
  Y <- as.data.frame(Y)
  
  X_mat <- as.matrix(X)
  if (is.null(dim(X_mat))) {
    X_mat <- matrix(X_mat, ncol = 1)
  }
  if (is.null(colnames(X_mat))) {
    colnames(X_mat) <- paste0("PC", seq_len(ncol(X_mat)))
  }
  colnames(X_mat) <- make.names(colnames(X_mat), unique = TRUE)
  X_data <- as.data.frame(X_mat)
  
  site_formula <- as.formula(paste("~", paste(colnames(X_data), collapse = " + ")))
  
  n_chains    <- DEFAULT_PARAMS$mcmc$n_chains
  beta_starts <- rep_len(c(0, -0.5, 0.5), n_chains)
  
  chain_job <- function(chain_id) {
    return(run_one_chain(Y            = Y,
                         X_data       = X_data,
                         site_formula = site_formula,
                         n_lv         = N_LV,
                         n_burnin     = DEFAULT_PARAMS$mcmc$n_burnin,
                         n_sample     = DEFAULT_PARAMS$mcmc$n_iter,
                         thin         = DEFAULT_PARAMS$mcmc$thin,
                         seed         = seed + chain_id,
                         beta_start   = beta_starts[chain_id]))
  }
  
  fits <- parallel::mclapply(seq_len(n_chains), chain_job, mc.cores = n_chains)
  
  nms         <- colnames(as.matrix(fits[[1]]$mcmc.sp[[1]]))
  beta_cols   <- grep("^beta", nms)
  lambda_cols <- grep("^lambda", nms)
  
  chains <- coda::mcmc.list(lapply(fits, function(fit) {
    out <- do.call(cbind, lapply(seq_along(fit$mcmc.sp), function(j) {
      m <- as.matrix(fit$mcmc.sp[[j]])[, beta_cols, drop = FALSE]
      colnames(m) <- c(sprintf("beta0[%d]", j),
                       if (ncol(m) > 1) {
                         sprintf("beta[%d,%d]", j, seq_len(ncol(m) - 1))
                       })
      return(m)
    }))
    return(coda::mcmc(out))
  }))
  rhat <- tryCatch(max(coda::gelman.diag(chains, multivariate = FALSE,
                                         autoburnin = FALSE)$psrf[, 1], na.rm = TRUE),
                   error = function(e) return(NA))
  neff <- tryCatch(min(coda::effectiveSize(chains), na.rm = TRUE),
                   error = function(e) return(NA))
  
  n_obs     <- nrow(Y)
  n_sp      <- length(fits[[1]]$mcmc.sp)
  n_covars  <- length(beta_cols) - 1
  per_chain <- nrow(as.matrix(fits[[1]]$mcmc.sp[[1]]))
  n_draws   <- per_chain * length(fits)
  
  beta0_draws  <- matrix(NA, n_draws, n_sp)
  beta_draws   <- array(NA, c(n_draws, n_sp, max(n_covars, 1)))
  lambda_draws <- array(NA, c(n_draws, n_sp, N_LV))
  
  for (c_id in seq_along(fits)) {
    rows <- ((c_id - 1) * per_chain + 1):(c_id * per_chain)
    for (j in seq_len(n_sp)) {
      m <- as.matrix(fits[[c_id]]$mcmc.sp[[j]])
      beta0_draws[rows, j] <- m[, beta_cols[1]]
      if (n_covars > 0) {
        beta_draws[rows, j, ] <- m[, beta_cols[-1], drop = FALSE]
      }
      lambda_draws[rows, j, ] <- m[, lambda_cols, drop = FALSE]
    }
  }
  lambda_draws[is.na(lambda_draws)] <- 0
  
  W_draws <- do.call(rbind, lapply(fits, function(fit) {
    return(as.matrix(fit$mcmc.latent[[1]]))
  }))
  W_mean <- matrix(colMeans(W_draws), n_obs, N_LV)
  
  return(list(beta0_draws   = beta0_draws,
              beta_draws    = beta_draws,
              lambda_draws  = lambda_draws,
              W_mean        = W_mean,
              n_draws       = n_draws,
              diagnostics   = list(rhat_max = rhat,
                                   neff_min = neff),
              mcmc_settings = list(n_burnin = DEFAULT_PARAMS$mcmc$n_burnin,
                                   n_sample = DEFAULT_PARAMS$mcmc$n_iter,
                                   thin     = DEFAULT_PARAMS$mcmc$thin,
                                   n_chains = n_chains)))
}

## SDM FITTING FUNCTION ########################################################

fit_SDM <- function(Y,
                    X,
                    cfg,
                    seed) {
  
  Y <- as.matrix(Y)
  if (is.null(colnames(Y))) {
    colnames(Y) <- paste0("sp", seq_len(ncol(Y)))
  }
  colnames(Y) <- make.names(colnames(Y), unique = TRUE)
  Y <- as.data.frame(Y)
  
  X_mat <- as.matrix(X)
  if (is.null(dim(X_mat))) {
    X_mat <- matrix(X_mat, ncol = 1)
  }
  if (is.null(colnames(X_mat))) {
    colnames(X_mat) <- paste0("PC", seq_len(ncol(X_mat)))
  }
  colnames(X_mat) <- make.names(colnames(X_mat), unique = TRUE)
  X_data <- as.data.frame(X_mat)
  
  site_formula <- as.formula(paste("~", paste(colnames(X_data), collapse = " + ")))
  
  n_chains    <- DEFAULT_PARAMS$mcmc$n_chains
  beta_starts <- rep_len(c(0, -0.5, 0.5), n_chains)
  
  chain_job <- function(chain_id) {
    return(run_one_chain(Y            = Y,
                         X_data       = X_data,
                         site_formula = site_formula,
                         n_lv         = 0,
                         n_burnin     = DEFAULT_PARAMS$mcmc$n_burnin,
                         n_sample     = DEFAULT_PARAMS$mcmc$n_iter,
                         thin         = DEFAULT_PARAMS$mcmc$thin,
                         seed         = seed + chain_id,
                         beta_start   = beta_starts[chain_id]))
  }
  
  fits <- parallel::mclapply(seq_len(n_chains), chain_job, mc.cores = n_chains)
  
  nms       <- colnames(as.matrix(fits[[1]]$mcmc.sp[[1]]))
  beta_cols <- grep("^beta", nms)
  
  chains <- coda::mcmc.list(lapply(fits, function(fit) {
    out <- do.call(cbind, lapply(seq_along(fit$mcmc.sp), function(j) {
      m <- as.matrix(fit$mcmc.sp[[j]])[, beta_cols, drop = FALSE]
      colnames(m) <- c(sprintf("beta0[%d]", j),
                       if (ncol(m) > 1) {
                         sprintf("beta[%d,%d]", j, seq_len(ncol(m) - 1))
                       })
      return(m)
    }))
    return(coda::mcmc(out))
  }))
  rhat <- tryCatch(max(coda::gelman.diag(chains, multivariate = FALSE,
                                         autoburnin = FALSE)$psrf[, 1], na.rm = TRUE),
                   error = function(e) return(NA))
  neff <- tryCatch(min(coda::effectiveSize(chains), na.rm = TRUE),
                   error = function(e) return(NA))
  
  n_sp      <- length(fits[[1]]$mcmc.sp)
  n_covars  <- length(beta_cols) - 1
  per_chain <- nrow(as.matrix(fits[[1]]$mcmc.sp[[1]]))
  n_draws   <- per_chain * length(fits)
  
  beta0_draws <- matrix(NA, n_draws, n_sp)
  beta_draws  <- array(NA, c(n_draws, n_sp, max(n_covars, 1)))
  
  for (c_id in seq_along(fits)) {
    rows <- ((c_id - 1) * per_chain + 1):(c_id * per_chain)
    for (j in seq_len(n_sp)) {
      m <- as.matrix(fits[[c_id]]$mcmc.sp[[j]])
      beta0_draws[rows, j] <- m[, beta_cols[1]]
      if (n_covars > 0) {
        beta_draws[rows, j, ] <- m[, beta_cols[-1], drop = FALSE]
      }
    }
  }
  
  return(list(beta0_draws   = beta0_draws,
              beta_draws    = beta_draws,
              n_draws       = n_draws,
              diagnostics   = list(rhat_max = rhat, neff_min = neff),
              mcmc_settings = list(n_burnin = DEFAULT_PARAMS$mcmc$n_burnin,
                                   n_sample = DEFAULT_PARAMS$mcmc$n_iter,
                                   thin     = DEFAULT_PARAMS$mcmc$thin,
                                   n_chains = n_chains)))
}

## RUN #########################################################################

t_start <- Sys.time()

data_raw <- data_loading(dataset_name)
Y        <- species_filtering(data_raw$Y, cfg$min_occ, cfg$max_occ)
pca      <- PCA(data_raw$X, cfg$n_pcs)
X_pca    <- pca$scores

n_sp          <- ncol(Y)
species_names <- colnames(Y)
if (is.null(species_names)) {
  species_names <- paste0("sp", 1:n_sp)
}

pca_var           <- pca$pca$sdev^2
pca_var_explained <- pca_var / sum(pca_var)

JSDM_results <- fit_JSDM(Y, X_pca, cfg, seed)
SDM_results  <- fit_SDM(Y, X_pca, cfg, seed)

model_info <- list(dataset_name      = dataset_name,
                   n_obs             = nrow(Y),
                   n_sp              = n_sp,
                   n_pcs             = pca$n_pcs,
                   n_lv              = N_LV,
                   species_names     = species_names,
                   pca_var_explained = pca_var_explained,
                   mcmc_params       = DEFAULT_PARAMS$mcmc,
                   jsdm_diagnostics  = JSDM_results$diagnostics,
                   sdm_diagnostics   = SDM_results$diagnostics,
                   seed              = seed,
                   total_time        = as.numeric(difftime(Sys.time(), t_start,
                                                           units = "secs")))

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

file_name <- file.path(results_dir, paste0(dataset_name, "_inference.RData"))

save(JSDM_results,
     SDM_results,
     model_info,
     Y,
     X_pca,
     file = file_name)

################################################################################