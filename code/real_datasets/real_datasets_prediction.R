################################################################################
##                                                                            ##
##                   REAL ECOLOGICAL DATASETS • PREDICTION                    ##
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

library(Hmsc)
library(tidyverse)
library(jSDM)
library(mvtnorm)
library(TMB)
library(runMCMCbtadjust)
library(coda)
library(parallel)
library(readxl)

## LOAD AND COMPILE TMB CODE ###################################################

TMB_code <- r"(
#include <TMB.hpp>
template<class Type>

Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(Y);
  DATA_MATRIX(X);
  DATA_MATRIX(beta);
  DATA_VECTOR(lambda);
  DATA_INTEGER(use_lv);

  PARAMETER_VECTOR(W);

  int n = Y.rows();
  int S = Y.cols();

  Type nll = Type(0.0);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < S; j++) {
      Type eta = Type(0.0);
      for (int k = 0; k < X.cols(); k++) {
        eta += X(i, k) * beta(k, j);
      }
      eta += lambda(j) * W(i);

      Type p = pnorm(eta, Type(0.0), Type(1.0));

      Type eps = Type(1.0e-12);
      p = CppAD::CondExpGt(p, Type(1.0) - eps, Type(1.0) - eps, p);
      p = CppAD::CondExpLt(p, eps, eps, p);

      nll -= Y(i, j) * log(p) + (Type(1.0) - Y(i, j)) * log(Type(1.0) - p);
    }
  }

  if (use_lv == 1) {
    for (int i = 0; i < n; i++) {
      nll -= dnorm(W(i), Type(0.0), Type(1.0), true);
    }
  }

  return nll;
})"

TMB_name <- "probit_marglik"
TMB_code <- paste0(TMB_name, ".cpp")
TMB::compile(TMB_code)
dyn.load(TMB::dynlib(TMB_name))

## CONFIGURATION ###############################################################

DATASETS <- c("butterflies", 
              "birds", 
              "eucalypts", 
              "kilpisjarvi", 
              "mara", 
              "aravo")

SPLITS <- c("interpolation", 
            "partial_extrapolation", 
            "full_extrapolation")

training_prop <- 0.5

CONFIG <- list(butterflies = list(n_pcs   = 4, 
                                  min_occ = 0.05, 
                                  max_occ = 0.95),
               birds       = list(n_pcs   = 5, 
                                  min_occ = 0.05, 
                                  max_occ = 0.95),
               eucalypts   = list(n_pcs   = 5, 
                                  min_occ = 0.05, 
                                  max_occ = 0.95),
               kilpisjarvi = list(n_pcs   = 3, 
                                  min_occ = 0.05, 
                                  max_occ = 0.95),
               mara        = list(n_pcs   = 5, 
                                  min_occ = 0.05, 
                                  max_occ = 0.95),
               aravo       = list(n_pcs   = 2, 
                                  min_occ = 0, 
                                  max_occ = 1))

DEFAULT_PARAMS <- list(mcmc = list(n_iter_min   = 10000,
                                   n_burnin_min = 10000,
                                   thin_min     = 10,
                                   conv_max     = 1.05,
                                   n_eff_min    = 5000,
                                   time_max     = 24))

## ARGUMENTS FROM SLURM ########################################################

ARGS <- commandArgs(trailingOnly = TRUE)

dataset_id <- as.integer(ARGS[1])
split_type <- as.integer(ARGS[2])
replicate  <- as.integer(ARGS[3])

dataset_name <- DATASETS[dataset_id]
split_name   <- SPLITS[split_type]
config       <- CONFIG[[dataset_name]]

## HELPER FUNCTIONS ############################################################

filter_species <- function(Y, min_occ, max_occ) {
  n <- nrow(Y)
  occ <- colSums(Y)
  min_c <- if (min_occ < 1) ceiling(n * min_occ) else min_occ
  max_c <- if (max_occ <= 1) floor(n * max_occ) else max_occ
  keep <- which(occ >= min_c & occ <= max_c)
  Y[, keep, drop = FALSE]
}

perform_pca <- function(X, n_pcs) {
  pca <- prcomp(scale(X), center = FALSE, scale. = FALSE)
  n_use <- min(n_pcs, ncol(X), nrow(X) - 1)
  pc <- pca$x[, 1:n_use, drop = FALSE]
  colnames(pc) <- paste0("PC", 1:n_use)
  list(scores = pc, 
       pc1    = pc[, 1], 
       n_pcs  = n_use)
}

interpolation_split <- function(Y, pc, prop, seed) {
  set.seed(seed)
  n <- nrow(Y)
  n_tr <- floor(n * prop)
  tr <- sample(1:n, n_tr)
  val <- setdiff(1:n, tr)
  list(Y_tr = Y[tr, ], 
       Y_val = Y[val, ], 
       X_tr = pc[tr, , drop = F],
       X_val = pc[val, , drop = F],
       tr_id = tr, 
       val_id = val)
}

partial_extrapolation_split <- function(Y, pc, pc1, prop, seed) {
  set.seed(seed)
  n <- nrow(Y)
  n_pairs <- floor(n / 2)
  shuf <- sample(1:n)
  tr <- val <- numeric(n_pairs)
  for (i in 1:n_pairs) {
    i1 <- shuf[2*i - 1]
    i2 <- shuf[2*i]
    if (pc1[i1] < pc1[i2]) { 
      tr[i] <- i1
      val[i] <- i2 
    } else { 
      tr[i] <- i2 
      val[i] <- i1 
    }
  }
  if (length(shuf) > 2 * n_pairs) {
    val <- c(val, shuf[2 * n_pairs + 1])
  }
  list(Y_tr = Y[tr, ], 
       Y_val = Y[val, ], 
       X_tr = pc[tr, , drop = F], 
       X_val = pc[val, , drop = F],
       tr_id = tr,
       val_id = val)
}

full_extrapolation_split <- function(Y, pc, pc1) {
  med <- median(pc1)
  tr <- which(pc1 < med)
  val <- which(pc1 >= med)
  list(Y_tr = Y[tr, ], 
       Y_val = Y[val, ], 
       X_tr = pc[tr, , drop = F], 
       X_val = pc[val, , drop = F],
       tr_id = tr, 
       val_id = val)
}

## DATASET LOADING #############################################################

load_data <- function(name) {

    if (name == "aravo") {
    
    data("aravo", package = "jSDM")
    
    Y <- apply(aravo$spe > 0, 2, as.numeric)
    Y <- Y[, colSums(Y) >= 5]
    X <- data.matrix(cbind(poly(aravo$env$Snow, 2)))
    
  } else if (name == "eucalypts") {
    
    data("eucalypts", package = "jSDM")
    
    Y <- ifelse(as.matrix(eucalypts[, 1:12]) > 0, 1, 0)
    X <- data.matrix(eucalypts[, 13:19])
    
  } else if (name == "birds") {
    
    spe <- read.csv('/path/to/datasets/datasets/birds/Birds_PA.csv')
    env <- read.csv('/path/to/datasets/datasets/birds/Birds_Cov.csv')
    
    Y <- ifelse(as.matrix(spe) > 0, 1, 0)
    X <- data.matrix(env)
    
  } else if (name == "butterflies") {
    
    spe <- read.csv('/path/to/datasets/butterflies/Butterfly_PA.csv')
    env <- read.csv('/path/to/datasets/butterflies/Butterfly_Cov.csv')
    
    Y <- ifelse(as.matrix(spe) > 0, 1, 0)
    X <- data.matrix(env[, -1])
    
  } else if (name == "kilpisjarvi") {
    
    data <- read.csv("/path/to/datasets/kilpisjarvi/Kilpisjarvi_plant_data.csv")
    
    sites <- unique(data$site)
    ny <- length(sites)
    sp <- unique(data$species)
    ns <- length(sp)
    Y.pa <- matrix(0, nrow = ny, ncol = ns)
    T3_GDD3 <- rep(NA, ny)
    T3_FDD <- rep(NA, ny)
    moist_mean_summer <- rep(NA, ny)
    
    for (k in 1:nrow(data)) {
      i <- which(data$site[k] == sites)
      j <- which(data$species[k] == sp)
      Y.pa[i, j] <- 1
      T3_GDD3[i] <- data$T3_GDD3[k]
      moist_mean_summer[i] <- data$moist_mean_summer[k]
      T3_FDD[i] <- data$T3_FDD[k]
    }
    
    sp.order    <- order(colSums(Y.pa), decreasing = TRUE)
    Y           <- Y.pa[, sp.order]
    colnames(Y) <- sp[sp.order]
    
    X <- data.frame(GDD = T3_GDD3, FDD = T3_FDD, SM = moist_mean_summer)
    X <- data.matrix(X)
    
  } else if (name == "mara") {
    
    data <- read_rds("/path/to/datasets/mara/mara_animal_compiled.rds") %>%
      mutate(Site = as.numeric(Site),
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
      dplyr::select(Pgrazed_lag1, Precip, Protein_lag1, Height_lag1, sin_month, cos_month)
    X <- data.matrix(XData)
    
    list(Y = Y, X = X)
  }
}

## JAGS MCMC ENGINE WITH runMCMCbtadjust() #####################################

run_mcmc <- function(model_string, 
                          data_list, 
                          inits_list, 
                          params,
                          seed, 
                          tag, 
                          mcmc_params = NULL, 
                           = NULL) {
  
  if (is.null(mcmc_params)) {
    mcmc_params <- DEFAULT_PARAMS$mcmc
  }
  
  model_file <- file.path(tempdir(), sprintf("%s.txt", tag))
  writeLines(model_string, model_file)
  
  out <- runMCMCbtadjust::runMCMC_btadjust(MCMC_language = "JAGS",
                                           code          = model_file,
                                           data          = data_list,
                                           inits         = inits_list,
                                           params        = params,
                                           params.conv   = params_conv,
                                           niter.min     = mcmc_params$n_iter_min,
                                           niter.max     = Inf,
                                           nburnin.min   = mcmc_params$n_burnin_min,
                                           nburnin.max   = Inf,
                                           thin.min      = mcmc_params$thin_min,
                                           thin.max      = Inf,
                                           Nchains       = length(inits_list),
                                           conv.max      = mcmc_params$conv_max,
                                           neff.min      = mcmc_params$n_eff_min,
                                           control       = list(time.max                   = 3600 * mcmc_params$time_max,
                                                                round.thinmult             = TRUE,
                                                                print.diagnostics          = FALSE,
                                                                Ncycles.target             = 2,
                                                                check.convergence.firstrun = FALSE,
                                                                convtype                   = "Gelman",
                                                                seed                       = seed),
                                           control.MCMC  = list(parallelize = TRUE))
  
  list(samples = out, attrs = attributes(out), file = model_file)
}

## JSDM FITTING FUNCTION #######################################################

fit_JSDM <- function(Y, X, config, seed) {
  
  Y     <- as.matrix(Y)
  X_mat <- as.matrix(X)
  if (is.vector(X_mat)) {
    X_mat <- matrix(X_mat, ncol = 1)
  }
  
  n_obs     <- nrow(Y)
  n_sp <- ncol(Y)
  n_covars  <- ncol(X_mat)
  n_lv    <- 1
  
  model_string <- "
    model {
      for(i in 1:n_obs) {
        for(j in 1:n_sp) {
          eta[i, j] <- inprod(lambda[j,], W[i,]) + inprod(beta[j,], X[i,])
          probit(p_y[i, j]) <- beta0[j] + eta[i, j]
          y[i, j] ~ dbern(p_y[i, j])
        }
      }

      for(i in 1:n_obs) {
        for(k in 1:n_lv) {
          W[i, k] ~ dnorm(0, 1)
        }
      }

      for(j in 1:n_sp) {
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

      for(j in 1:n_sp) {
        for(m in 1:n_covars) {
          beta[j, m] ~ dnorm(0, 0.1)
        }
      }
    }"
  
  data_list <- list(y        = Y,
                    X        = X_mat,
                    n_obs    = n_obs,
                    n_sp     = n_sp,
                    n_lv     = n_lv,
                    n_covars = n_covars)
  
  make_inits <- function() {
    list(beta0 = rnorm(n_sp, 0, 0.5),
         beta  = matrix(rnorm(n_sp * n_covars, 0, 0.3), n_sp, n_covars),
         W     = matrix(rnorm(n_obs * n_lv, 0, 0.5), n_obs, n_lv))
  }
  
  inits_list <- replicate(3, make_inits(), simplify = FALSE)
  
  run <- run_mcmc(model_string = model_string,
                  data_list    = data_list,
                  inits_list   = inits_list,
                  params       = c("beta0", "beta", "lambda", "W"),
                  params_conv  = c("beta0", "beta"),
                  seed         = seed,
                  tag          = "JSDM_inference")
  
  combined <- do.call(rbind, run$samples)
  attrs    <- run$attrs
  
  beta0_draws <- matrix(NA, nrow(combined), n_sp)
  beta0_cols <- grep("^beta0\\[", colnames(combined))
  for (id in beta0_cols) {
    j <- as.numeric(gsub("beta0\\[|\\]", "", colnames(combined)[id]))
    beta0_draws[, j] <- combined[, id]
  }
  
  beta_draws <- array(NA, dim = c(nrow(combined), n_sp, n_covars))
  beta_cols <- grep("^beta\\[", colnames(combined))
  for (id in beta_cols) {
    nm <- colnames(combined)[id]
    ij <- as.numeric(strsplit(gsub("beta\\[|\\]", "", nm), ",")[[1]])
    beta_draws[, ij[1], ij[2]] <- combined[, id]
  }
  
  lambda_draws <- array(NA, dim = c(nrow(combined), n_sp, n_lv))
  lambda_cols <- grep("^lambda\\[", colnames(combined))
  for (id in lambda_cols) {
    nm <- colnames(combined)[id]
    ij <- as.numeric(strsplit(gsub("lambda\\[|\\]", "", nm), ",")[[1]])
    lambda_draws[, ij[1], ij[2]] <- combined[, id]
  }
  lambda_draws[is.na(lambda_draws)] <- 0
  
  W_draws <- array(NA, dim = c(nrow(combined), n_obs, n_lv))
  W_cols <- grep("^W\\[", colnames(combined))
  for (id in W_cols) {
    nm <- colnames(combined)[id]
    ik <- as.numeric(strsplit(gsub("W\\[|\\]", "", nm), ",")[[1]])
    W_draws[, ik[1], ik[2]] <- combined[, id]
  }
  W_mean <- apply(W_draws, c(2, 3), mean, na.rm = TRUE)
  
  unlink(run$file)
  
  list(beta0_draws  = beta0_draws,
       beta_draws   = beta_draws,
       lambda_draws = lambda_draws,
       W_mean       = W_mean,
       n_draws      = nrow(combined),
       convergence  = attrs$final.params$converged)
}

## SDM FITTING FUNCTION ########################################################

fit_SDM <- function(Y, X, config, seed) {
  
  Y     <- as.matrix(Y)
  X_mat <- as.matrix(X)
  if (is.vector(X_mat)) {
    X_mat <- matrix(X_mat, ncol = 1)
  }
  
  n_obs     <- nrow(Y)
  n_sp <- ncol(Y)
  n_covars  <- ncol(X_mat)
  
  model_string <- "
    model {
      for(i in 1:n_obs) {
        for(j in 1:n_sp) {
          eta[i, j] <- inprod(beta[j,], X[i,])
          probit(p_y[i, j]) <- beta0[j] + eta[i, j]
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
  
  data_list <- list(y        = Y,
                    X        = X_mat,
                    n_obs    = n_obs,
                    n_sp     = n_sp,
                    n_covars = n_covars)
  
  make_inits <- function() {
    list(beta0 = rnorm(n_sp, 0, 0.5),
         beta  = matrix(rnorm(n_sp * n_covars, 0, 0.3), n_sp, n_covars))
  }
  
  inits_list <- replicate(3, make_inits(), simplify = FALSE)
  
  run <- run_mcmc(model_string = model_string,
                  data_list    = data_list,
                  inits_list   = inits_list,
                  params       = c("beta0", "beta"),
                  params_conv  = c("beta0", "beta"),
                  seed         = seed,
                  tag          = "SDM_inference")
  
  combined <- do.call(rbind, run$samples)
  attrs    <- run$attrs
  
  beta0_draws <- matrix(NA, nrow(combined), n_sp)
  beta0_cols <- grep("^beta0\\[", colnames(combined))
  for (id in beta0_cols) {
    j <- as.numeric(gsub("beta0\\[|\\]", "", colnames(combined)[id]))
    beta0_draws[, j] <- combined[, id]
  }
  
  beta_draws <- array(NA, dim = c(nrow(combined), n_sp, n_covars))
  beta_cols <- grep("^beta\\[", colnames(combined))
  for (id in beta_cols) {
    nm <- colnames(combined)[id]
    ij <- as.numeric(strsplit(gsub("beta\\[|\\]", "", nm), ",")[[1]])
    beta_draws[, ij[1], ij[2]] <- combined[, id]
  }
  
  unlink(run$file)
  
  list(beta0_draws = beta0_draws,
       beta_draws  = beta_draws,
       n_draws     = nrow(combined),
       convergence = attrs$final.params$converged)
}

## LOG-LIKELIHOOD COMPUTATION ##################################################

compute_loglik <- function(Y_test, 
                           X_test, 
                           posterior,
                           model_type, 
                           dll,
                           n_draws_target = 10000,
                           seed) {
  
  n_test      <- nrow(Y_test)
  S           <- ncol(Y_test)
  X_design    <- cbind(1, X_test)
  n_available <- posterior$n_draws
  
  if (n_available >= n_draws_target) {
    draw_id <- sample(1:n_available, n_draws_target, replace = FALSE)
  } else {
    draw_id <- 1:n_available
  }
  n_draws <- length(draw_id)
  
  if (model_type == "jsdm") {
    use_lv     <- 1L
    random_arg <- "W"
  } else {
    use_lv     <- 0L
    random_arg <- NULL
  }
  
  loglik_mat <- matrix(NA, n_draws, n_test)
  
  for (d_pos in seq_along(draw_id)) {
    d <- draw_id[d_pos]
    
    beta_coef <- posterior$beta_draws[d, , ]
    if (is.null(dim(beta_coef))) {
      beta_coef <- matrix(beta_coef, nrow = 1)
    } else {
      beta_coef <- t(beta_coef)
    }
    beta_d <- rbind(posterior$beta0_draws[d, ], beta_coef)
    
    if (model_type == "jsdm") {
      lambda_d <- posterior$lambda_draws[d, , 1]
    } else {
      lambda_d <- rep(0, S)
    }
    
    for (i in seq_len(n_test)) {
      data_i <- list(Y      = matrix(Y_test[i, ], nrow = 1),
                     X      = matrix(X_design[i, ], nrow = 1),
                     beta   = beta_d,
                     lambda = as.numeric(lambda_d),
                     use_lv = use_lv)
      pars_i <- list(W = 0)
      obj_i  <- MakeADFun(data_i, 
                          pars_i, 
                          random = random_arg,
                          DLL    = dll,
                          silent = TRUE)
      loglik_mat[d_pos, i] <- -as.numeric(obj_i$fn())
    }
  }
  
  loglik_site <- numeric(n_test)
  for (i in seq_len(n_test)) {
    ll_vec <- loglik_mat[, i]
    max_ll <- max(ll_vec)
    loglik_site[i] <- max_ll + log(mean(exp(ll_vec - max_ll)))
  }
  
  loglik_total <- sum(loglik_site)
  
  list(total    = loglik_total,
       per_site = loglik_site)
}

## RUN #########################################################################

data <- load_data(dataset_name)
Y    <- filter_species(data$Y, config$min_occ, config$max_occ)
pca  <- perform_pca(data$X, config$n_pcs)

n_sp          <- ncol(Y)
species_names <- colnames(Y)

if (is.null(species_names)) {
  species_names <- paste0("sp", 1:n_sp)
}

seed <- 1000 + (dataset_id - 1) * 10000 + (split_type - 1) * 1000 + replicate

if (split_type == 1) {
  split <- interpolation_split(Y, pca$scores, training_prop, seed)
} else if (split_type == 2) {
  split <- partial_extrapolation_split(Y, pca$scores, pca$pc1, training_prop, seed)
} else {
  split <- full_extrapolation_split(Y, pca$scores, pca$pc1)
}

JSDM_result <- fit_JSDM(split$Y_tr, split$X_tr, config, seed)
SDM_result  <- fit_SDM(split$Y_tr, split$X_tr, config, seed)

JSDM_loglik <- compute_loglik(split$Y_val,
                                    split$X_val,
                                    JSDM_result,
                                    "jsdm", 
                                    TMB_name,
                                    n_draws_target = 10000, 
                                    seed = seed)

SDM_loglik  <- compute_loglik(split$Y_val, 
                                    split$X_val,
                                    SDM_result,
                                    "sdm", 
                                    TMB_name,
                                    n_draws_target = 10000,
                                    seed = seed)

marginal_loglik <- list(JSDM = JSDM_loglik,
                        SDM  = SDM_loglik)

posterior_summary <- list(JSDM = list(beta0_mean  = colMeans(JSDM_result$beta0_draws),
                                      beta_mean   = apply(JSDM_result$beta_draws, c(2, 3), mean),
                                      lambda_mean = apply(JSDM_result$lambda_draws, c(2, 3), mean),
                                      n_draws     = JSDM_result$n_draws),
                          SDM = list(beta0_mean = colMeans(SDM_result$beta0_draws),
                                     beta_mean  = apply(SDM_result$beta_draws, c(2, 3), mean),
                                     n_draws    = SDM_result$n_draws))


dir.create("results", showWarnings = FALSE)
dir.create(file.path("results", dataset_name), showWarnings = FALSE)

file_name <- file.path("results", dataset_name,
                       sprintf("%s_rep%03d.RData", split_name, replicate))

save(dataset_name, 
     split_name, 
     replicate, 
     seed,
     marginal_loglik, 
     posterior_summary,
     n_sp, 
     species_names,
     file = file_name)

################################################################################
