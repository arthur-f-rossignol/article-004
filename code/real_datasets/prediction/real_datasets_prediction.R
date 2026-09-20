################################################################################
##                                                                            ##
##                   REAL ECOLOGICAL DATASETS • PREDICTION                    ##
##           (real_datasets/prediction/real_datasets_prediction.R)            ##
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
library(TMB)
library(coda)
library(parallel)
library(readxl)

## ARGUMENTS FROM SLURM ########################################################

ARGS <- commandArgs(trailingOnly = TRUE)

dataset_id <- as.integer(ARGS[1])
split_type <- as.integer(ARGS[2])
replicate  <- as.integer(ARGS[3])

base_dir    <- ARGS[4]
data_dir    <- file.path(base_dir, "datasets")
results_dir <- file.path(base_dir, "prediction", "results")

## CONFIGURATION ###############################################################

DATASETS <- c("dataset.1",
              "dataset.2",
              "dataset.3",
              "dataset.4",
              "dataset.5",
              "dataset.6")

SPLITS <- c("interpolation",
            "partial_extrapolation",
            "full_extrapolation")

TRAIN_PROP   <- 0.5
N_REPLICATES <- 100

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

DEFAULT_PARAMS <- list(mcmc = list(niter    = 100000,
                                   nburnin  = 50000,
                                   thin     = 10,
                                   n_chains = 3),
                       prior = list(mu_beta   = 0,
                                    V_beta    = 10,
                                    mu_lambda = 0,
                                    V_lambda  = 10))

N_DRAWS_TARGET <- 2000

dataset_name <- DATASETS[dataset_id]
split_name   <- SPLITS[split_type]
cfg          <- CONFIGS[[dataset_name]]

seed <- 1000 + (dataset_id - 1) * 10000 + (split_type - 1) * 1000 + replicate

## SDM TMB MODEL ###############################################################

TMB_SRC_SDM <- r"(
  #include <TMB.hpp>

  template<class Type>
  Type log_pnorm_safe(Type z) {
    Type cut  = Type(-20.0);
    Type z_hi = CppAD::CondExpLt(z, cut, cut, z);
    Type z_lo = CppAD::CondExpLt(z, cut, z, cut);

    Type direct = log(pnorm(z_hi));

    Type iz2  = Type(1.0) / (z_lo * z_lo);
    Type ser  = Type(1.0) - iz2 + Type(3.0) * iz2 * iz2
                - Type(15.0) * iz2 * iz2 * iz2;
    Type asym = dnorm(z_lo, Type(0.0), Type(1.0), true) - log(-z_lo) + log(ser);

    return CppAD::CondExpLt(z, cut, asym, direct);
  }

  template<class Type>
  Type objective_function<Type>::operator() ()
  {
    DATA_MATRIX(Y);
    DATA_MATRIX(X);

    PARAMETER_MATRIX(beta);

    int n = Y.rows();
    int S = Y.cols();
    int K = X.cols();

    vector<Type> nll_site(n);
    nll_site.setZero();

    Type nll = Type(0.0);

    for (int i = 0; i < n; i++) {

      Type nll_i = Type(0.0);

      for (int j = 0; j < S; j++) {

        Type eta = Type(0.0);
        for (int k = 0; k < K; k++) {
          eta += X(i, k) * beta(k, j);
        }

        Type s = Type(2.0) * Y(i, j) - Type(1.0);
        Type u = s * eta;

        nll_i -= log_pnorm_safe(u);
      }

      nll_site(i) = nll_i;
      nll        += nll_i;
    }

    REPORT(nll_site);

    return nll;
  })"

TMB_NAME_SDM <- "probit_marglik_SDM"

## JSDM TMB MODEL ##############################################################

TMB_SRC_JSDM <- r"(
  #include <TMB.hpp>

  template<class Type>
  Type log_pnorm_safe(Type z) {
    Type cut  = Type(-20.0);
    Type z_hi = CppAD::CondExpLt(z, cut, cut, z);
    Type z_lo = CppAD::CondExpLt(z, cut, z, cut);

    Type direct = log(pnorm(z_hi));

    Type iz2  = Type(1.0) / (z_lo * z_lo);
    Type ser  = Type(1.0) - iz2 + Type(3.0) * iz2 * iz2
                - Type(15.0) * iz2 * iz2 * iz2;
    Type asym = dnorm(z_lo, Type(0.0), Type(1.0), true) - log(-z_lo) + log(ser);

    return CppAD::CondExpLt(z, cut, asym, direct);
  }

  template<class Type>
  Type objective_function<Type>::operator() ()
  {
    DATA_MATRIX(Y);
    DATA_MATRIX(X);

    PARAMETER_MATRIX(beta);
    PARAMETER_VECTOR(lambda);
    PARAMETER_VECTOR(W);

    int n = Y.rows();
    int S = Y.cols();
    int K = X.cols();

    vector<Type> nll_site(n);
    vector<Type> hess_w(n);
    nll_site.setZero();
    hess_w.setZero();

    Type nll = Type(0.0);

    for (int i = 0; i < n; i++) {

      Type nll_i = Type(0.0);
      Type h_i   = Type(0.0);

      for (int j = 0; j < S; j++) {

        Type eta = Type(0.0);
        for (int k = 0; k < K; k++) {
          eta += X(i, k) * beta(k, j);
        }
        eta += lambda(j) * W(i);

        Type s = Type(2.0) * Y(i, j) - Type(1.0);
        Type u = s * eta;

        Type lp = log_pnorm_safe(u);
        nll_i -= lp;

        Type r = exp(dnorm(u, Type(0.0), Type(1.0), true) - lp);
        h_i += lambda(j) * lambda(j) * (u * r + r * r);
      }

      nll_i -= dnorm(W(i), Type(0.0), Type(1.0), true);
      h_i   += Type(1.0);

      nll_site(i) = nll_i;
      hess_w(i)   = h_i;
      nll        += nll_i;
    }

    REPORT(nll_site);
    REPORT(hess_w);

    return nll;
  })"

TMB_NAME_JSDM <- "probit_marglik_JSDM"

## TMB COMPILATION #############################################################

setup_tmb <- function(TMB_name,
                      TMB_src) {
  tmb_model_dir <- file.path(tempdir(), paste0("build_", Sys.getpid()))
  dir.create(tmb_model_dir, showWarnings = FALSE, recursive = TRUE)
  cpp <- file.path(tmb_model_dir, paste0(TMB_name, ".cpp"))
  writeLines(TMB_src, cpp)
  TMB::compile(cpp)
  dyn.load(TMB::dynlib(file.path(tmb_model_dir, TMB_name)))
  return(TMB_name)
}

setup_tmb(TMB_NAME_SDM, TMB_SRC_SDM)
setup_tmb(TMB_NAME_JSDM, TMB_SRC_JSDM)

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

log_sum_exp_rows <- function(M) {
  row_max <- apply(M, 1, max)
  return(row_max + log(rowSums(exp(M - row_max))))
}

interpolation_split <- function(Y,
                                pc,
                                prop,
                                seed) {
  set.seed(seed)
  n    <- nrow(Y)
  n_tr <- floor(n * prop)
  tr   <- sample(1:n, n_tr)
  val  <- setdiff(1:n, tr)
  return(list(Y_tr   = Y[tr, , drop = FALSE],
              Y_val  = Y[val, , drop = FALSE],
              X_tr   = pc[tr, , drop = FALSE],
              X_val  = pc[val, , drop = FALSE],
              tr_id  = tr,
              val_id = val))
}

partial_extrapolation_split <- function(Y,
                                        pc,
                                        pc1,
                                        prop,
                                        seed) {
  set.seed(seed)
  n       <- nrow(Y)
  n_pairs <- floor(n / 2)
  shuf    <- sample(1:n)
  tr <- val <- numeric(n_pairs)
  for (i in 1:n_pairs) {
    i1 <- shuf[2 * i - 1]
    i2 <- shuf[2 * i]
    if (pc1[i1] < pc1[i2]) {
      tr[i]  <- i1
      val[i] <- i2
    } else {
      tr[i]  <- i2
      val[i] <- i1
    }
  }
  if (length(shuf) > 2 * n_pairs) {
    val <- c(val, shuf[2 * n_pairs + 1])
  }
  return(list(Y_tr   = Y[tr, , drop = FALSE],
              Y_val  = Y[val, , drop = FALSE],
              X_tr   = pc[tr, , drop = FALSE],
              X_val  = pc[val, , drop = FALSE],
              tr_id  = tr,
              val_id = val))
}

full_extrapolation_split <- function(Y,
                                     pc,
                                     pc1) {
  med <- median(pc1)
  tr  <- which(pc1 < med)
  val <- which(pc1 >= med)
  return(list(Y_tr   = Y[tr, , drop = FALSE],
              Y_val  = Y[val, , drop = FALSE],
              X_tr   = pc[tr, , drop = FALSE],
              X_val  = pc[val, , drop = FALSE],
              tr_id  = tr,
              val_id = val))
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
  
  fit  <- do.call(jSDM::jSDM_binomial_probit, jSDM_args)
  
  slim <- list(mcmc.sp       = fit$mcmc.sp,
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
                         n_burnin     = DEFAULT_PARAMS$mcmc$nburnin,
                         n_sample     = DEFAULT_PARAMS$mcmc$niter,
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
  
  n_sp      <- length(fits[[1]]$mcmc.sp)
  n_covars  <- length(beta_cols) - 1
  per_chain <- nrow(as.matrix(fits[[1]]$mcmc.sp[[1]]))
  n_draws   <- per_chain * length(fits)
  
  beta0_draws  <- matrix(NA, n_draws, n_sp)
  beta_draws   <- array(NA, c(n_draws, n_sp, max(n_covars, 1)))
  lambda_draws <- array(NA, c(n_draws, n_sp, 1))
  
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
  
  return(list(beta0_draws   = beta0_draws,
              beta_draws    = beta_draws,
              lambda_draws  = lambda_draws,
              n_draws       = n_draws,
              diagnostics   = list(rhat_max = rhat,
                                   neff_min = neff),
              mcmc_settings = list(n_burnin = DEFAULT_PARAMS$mcmc$nburnin,
                                   n_sample = DEFAULT_PARAMS$mcmc$niter,
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
                         n_burnin     = DEFAULT_PARAMS$mcmc$nburnin,
                         n_sample     = DEFAULT_PARAMS$mcmc$niter,
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
              mcmc_settings = list(n_burnin = DEFAULT_PARAMS$mcmc$nburnin,
                                   n_sample = DEFAULT_PARAMS$mcmc$niter,
                                   thin     = DEFAULT_PARAMS$mcmc$thin,
                                   n_chains = n_chains)))
}

## JSDM LOG-LIKELIHOOD COMPUTATION #############################################

compute_loglik_JSDM <- function(Y_test,
                                X_test,
                                posterior,
                                n_draws_target,
                                seed = 1) {
  
  Y_test   <- as.matrix(Y_test)
  X_test   <- as.matrix(X_test)
  n_test   <- nrow(Y_test)
  n_sp     <- ncol(Y_test)
  n_covars <- ncol(X_test)
  
  set.seed(seed)
  n_available <- posterior$n_draws
  draw_id <- if (n_available > n_draws_target) {
    sort(sample.int(n_available, n_draws_target))
  } else {
    seq_len(n_available)
  }
  n_draws <- length(draw_id)
  
  obj <- TMB::MakeADFun(data       = list(Y = Y_test,
                                          X = cbind(1, X_test)),
                        parameters = list(beta   = matrix(0, n_covars + 1, n_sp),
                                          lambda = rep(0, n_sp),
                                          W      = rep(0, n_test)),
                        random     = "W",
                        DLL        = TMB_NAME_JSDM,
                        silent     = TRUE)
  
  par_init  <- obj$par
  is_beta   <- (names(par_init) == "beta")
  is_lambda <- (names(par_init) == "lambda")
  
  log_lik_draws <- matrix(NA, n_draws, n_test)
  
  for (d_pos in seq_len(n_draws)) {
    
    d           <- draw_id[d_pos]
    beta0_draw  <- posterior$beta0_draws[d, ]
    beta_draw   <- array(posterior$beta_draws[d, , ], c(n_sp, n_covars))
    lambda_draw <- as.numeric(posterior$lambda_draws[d, , 1])
    
    par            <- par_init
    par[is_beta]   <- as.vector(rbind(beta0_draw, t(beta_draw)))
    par[is_lambda] <- lambda_draw
    
    obj$fn(par)
    reported <- obj$report(obj$env$last.par)
    
    log_lik_draws[d_pos, ] <- -reported$nll_site + 0.5 * log(2 * pi) -
      0.5 * log(reported$hess_w)
  }
  
  loglik_site <- log_sum_exp_rows(t(log_lik_draws)) - log(n_draws)
  
  return(list(total        = sum(loglik_site),
              mean_by_site = mean(loglik_site),
              per_site     = loglik_site,
              n_draws_used = n_draws))
}

## SDM LOG-LIKELIHOOD COMPUTATION ##############################################

compute_loglik_SDM <- function(Y_test,
                               X_test,
                               posterior,
                               n_draws_target,
                               seed = 1) {
  
  Y_test   <- as.matrix(Y_test)
  X_test   <- as.matrix(X_test)
  n_test   <- nrow(Y_test)
  n_sp     <- ncol(Y_test)
  n_covars <- ncol(X_test)
  
  set.seed(seed)
  n_available <- posterior$n_draws
  draw_id <- if (n_available > n_draws_target) {
    sort(sample.int(n_available, n_draws_target))
  } else {
    seq_len(n_available)
  }
  n_draws <- length(draw_id)
  
  obj <- TMB::MakeADFun(data       = list(Y = Y_test,
                                          X = cbind(1, X_test)),
                        parameters = list(beta = matrix(0, n_covars + 1, n_sp)),
                        DLL        = TMB_NAME_SDM,
                        silent     = TRUE)
  
  log_lik_draws <- matrix(NA, n_draws, n_test)
  
  for (d_pos in seq_len(n_draws)) {
    
    d          <- draw_id[d_pos]
    beta0_draw <- posterior$beta0_draws[d, ]
    beta_draw  <- array(posterior$beta_draws[d, , ], c(n_sp, n_covars))
    
    obj$fn(as.vector(rbind(beta0_draw, t(beta_draw))))
    reported <- obj$report(obj$env$last.par)
    
    log_lik_draws[d_pos, ] <- -reported$nll_site
  }
  
  loglik_site <- log_sum_exp_rows(t(log_lik_draws)) - log(n_draws)
  
  return(list(total        = sum(loglik_site),
              mean_by_site = mean(loglik_site),
              per_site     = loglik_site,
              n_draws_used = n_draws))
}

## RUN #########################################################################

t_start <- Sys.time()

data <- data_loading(dataset_name)
Y    <- species_filtering(data$Y, cfg$min_occ, cfg$max_occ)
pca  <- PCA(data$X, cfg$n_pcs)

n_sp          <- ncol(Y)
species_names <- colnames(Y)
if (is.null(species_names)) {
  species_names <- paste0("sp", 1:n_sp)
}

split <- if (split_type == 1) {
  interpolation_split(Y, pca$scores, TRAIN_PROP, seed)
} else if (split_type == 2) {
  partial_extrapolation_split(Y, pca$scores, pca$pc1, TRAIN_PROP, seed)
} else {
  full_extrapolation_split(Y, pca$scores, pca$pc1)
}

JSDM_result <- fit_JSDM(split$Y_tr, split$X_tr, cfg, seed)
SDM_result  <- fit_SDM(split$Y_tr, split$X_tr, cfg, seed)

JSDM_loglik <- compute_loglik_JSDM(split$Y_val,
                                   split$X_val,
                                   JSDM_result,
                                   N_DRAWS_TARGET,
                                   seed)
SDM_loglik  <- compute_loglik_SDM(split$Y_val,
                                  split$X_val,
                                  SDM_result,
                                  N_DRAWS_TARGET,
                                  seed)

marginal_loglik <- list(JSDM = JSDM_loglik,
                        SDM = SDM_loglik,
                        diff_total = JSDM_loglik$total - SDM_loglik$total,
                        diff_per_site = JSDM_loglik$per_site - SDM_loglik$per_site)

posterior_summary <- list(
  JSDM = list(beta0_mean   = colMeans(JSDM_result$beta0_draws),
              beta0_sd     = apply(JSDM_result$beta0_draws, 2, sd),
              beta_mean    = apply(JSDM_result$beta_draws, c(2, 3), mean),
              beta_sd      = apply(JSDM_result$beta_draws, c(2, 3), sd),
              lambda_mean  = apply(JSDM_result$lambda_draws, c(2, 3), mean),
              lambda_sd    = apply(JSDM_result$lambda_draws, c(2, 3), sd),
              n_draws      = JSDM_result$n_draws),
  SDM  = list(beta0_mean = colMeans(SDM_result$beta0_draws),
              beta0_sd   = apply(SDM_result$beta0_draws, 2, sd),
              beta_mean  = apply(SDM_result$beta_draws, c(2, 3), mean),
              beta_sd    = apply(SDM_result$beta_draws, c(2, 3), sd),
              n_draws    = SDM_result$n_draws))

model_info <- list(engine           = "jSDM::jSDM_binomial_probit",
                   n_lv             = N_LV,
                   n_train          = nrow(split$Y_tr),
                   n_val            = nrow(split$Y_val),
                   jsdm_diagnostics = JSDM_result$diagnostics,
                   sdm_diagnostics  = SDM_result$diagnostics,
                   jsdm_settings    = JSDM_result$mcmc_settings,
                   sdm_settings     = SDM_result$mcmc_settings,
                   total_time       = as.numeric(difftime(Sys.time(), t_start,
                                                          units = "secs")))

dir.create(file.path(results_dir, dataset_name),
           showWarnings = FALSE, recursive = TRUE)

file_name <- file.path(results_dir, dataset_name,
                       sprintf("%s_rep%03d.RData", split_name, replicate))

save(dataset_name,
     split_name,
     replicate,
     seed,
     marginal_loglik,
     posterior_summary,
     model_info,
     n_sp,
     species_names,
     file = file_name)

################################################################################