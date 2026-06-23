################################################################################
##                                                                            ##
##                      SDMs WITH ONE MISSING COVARIATE                       ##
##                                                                            ##
##       Addressing Missing Covariates in Species Distribution Models:        ##
##  Inferential Impacts and Mitigation via Joint Species Distribution Models  ##
##                                                                            ##
##                  Arthur F. Rossignol & Frederic Gosselin                   ##
##                                                                            ##
##                                    2026                                    ##
##                                                                            ##
################################################################################

# SDM.1   lm | glm | lmer | glmer   stats | lme4    frequentist   no-ORLE + ORLE
# SDM.2   glmmTMB                   glmmTMB         frequentist   no-ORLE + ORLE
# SDM.3   fitme                     spaMM           frequentist   no-ORLE + ORLE
# SDM.4   mixed_model               GLMMadaptive    frequentist             ORLE
# SDM.5   brm                       brms            Bayesian      no-ORLE + ORLE
# SDM.6   inla                      INLA            Bayesian      no-ORLE + ORLE
# SDM.7                             nimble          Bayesian      no-ORLE + ORLE
# SDM.8                             runjags         Bayesian      no-ORLE + ORLE
# SDM.9   lmrob | glmrob            robustbase      frequentist   no-ORLE
# SDM.10  glm (quasi-likelihood)    stats           frequentist   no-ORLE 
# SDM.11  gam | gamm                mgcv            frequentist   no-ORLE + ORLE

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

args <- commandArgs(trailingOnly = TRUE)

seed_value <- as.integer(args[1])
output_dir <- args[2]

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
                 "SDM.9_OLRE",
                 "SDM.10-QP",
                 "SDM.10-QB",
                 "SDM.11_no_OLRE",
                 "SDM.11_OLRE")

method_to_sdm <- c(SDM.1     = "SDM.1",
                   SDM.2     = "SDM.2",
                   SDM.3     = "SDM.3",
                   SDM.4     = "SDM.4",
                   SDM.5     = "SDM.5",
                   SDM.6     = "SDM.6",
                   SDM.7     = "SDM.7",
                   SDM.8     = "SDM.8",
                   SDM.9     = "SDM.9",
                   SDM.10_QP = "SDM.10",
                   SDM.10_QB = "SDM.10",
                   SDM.11    = "SDM.11")

method_compatibility <- list(SDM.1     = c("gaussian", 
                                           "poisson",
                                           "bernoulli_probit", 
                                           "bernoulli_logit",
                                           "bernoulli_cloglog"),
                             SDM.2     = c("gaussian", 
                                           "poisson",
                                           "bernoulli_probit",
                                           "bernoulli_logit",
                                           "bernoulli_cloglog"),
                             SDM.3     = c("gaussian", 
                                           "poisson",
                                           "bernoulli_probit", 
                                           "bernoulli_logit",
                                           "bernoulli_cloglog"),
                             SDM.4     = c("poisson",
                                           "bernoulli_probit", 
                                           "bernoulli_logit"),
                             SDM.5     = c("gaussian", 
                                           "poisson",
                                           "bernoulli_probit", 
                                           "bernoulli_logit"),
                             SDM.6     = c("gaussian", 
                                           "poisson",
                                           "bernoulli_probit", 
                                           "bernoulli_logit",
                                           "bernoulli_cloglog"),
                             SDM.7     = c("gaussian", 
                                           "poisson",
                                           "bernoulli_probit", 
                                           "bernoulli_logit",
                                           "bernoulli_cloglog"),
                             SDM.8     = c("gaussian", 
                                           "poisson",
                                           "bernoulli_probit", 
                                           "bernoulli_logit",
                                           "bernoulli_cloglog"),
                             SDM.9     = c("gaussian", 
                                           "poisson",
                                           "bernoulli_probit", 
                                           "bernoulli_logit",
                                           "bernoulli_cloglog"),
                             SDM.10_QP = c("poisson"),
                             SDM.10_QB = c("bernoulli_probit", 
                                           "bernoulli_logit",
                                           "bernoulli_cloglog"),
                             SDM.11    = c("gaussian", 
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

generate_covariates <- function(n, seed) {
  
  set.seed(seed)
  
  X1 <- rnorm(n, 0, 1)
  X2 <- rnorm(n, 0, 1)
  
  return(list(X1 = X1, X2 = X2))
}

generate_data <- function(n, X1, X2, model_type, seed) {
  
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

## GET FAMILY FUNCTION #########################################################

get_family <- function(model_type) {
  switch(model_type,
         "gaussian"          = gaussian(),
         "poisson"           = poisson(link = "log"),
         "bernoulli_probit"  = binomial(link = "probit"),
         "bernoulli_logit"   = binomial(link = "logit"),
         "bernoulli_cloglog" = binomial(link = "cloglog"))
}

## SDM.1 #######################################################################

fit_SDM.1 <- function(Y, X1, model_type, use_OLRE) {
  
  obs_id <- factor(seq_along(Y))
  data_df <- data.frame(Y = Y, X1 = X1, obs_id = obs_id)
  
  family_spec <- get_family(model_type)
  
  if (use_OLRE) {
    if (model_type == "gaussian") {
      fit <- lmer(Y ~ X1 + (1 | obs_id), data = data_df)
    } else {
      fit <- glmer(Y ~ X1 + (1 | obs_id), family = family_spec, data = data_df)
    }
  } else {
    if (model_type == "gaussian") {
      fit <- lm(Y ~ X1, data = data_df)
    } else {
      fit <- glm(Y ~ X1, family = family_spec, data = data_df)
    }
  }
  
  summary_fit <- summary(fit)
  
  if (use_OLRE && (model_type != "gaussian" || inherits(fit, "lmerMod"))) {
    result <- list(alpha_est = summary_fit$coefficients[1, 1],
                   alpha_SE  = summary_fit$coefficients[1, 2],
                   beta_est  = summary_fit$coefficients[2, 1],
                   beta_SE   = summary_fit$coefficients[2, 2])
    if (inherits(fit, "lmerMod") || inherits(fit, "glmerMod")) {
      result$sigma2_est <- summary_fit$varcor$obs_id[1]
      ci <- confint(fit, parm = "theta_", method = "Wald", quiet = TRUE)
      theta_est <- sqrt(result$sigma2_est)
      theta_se  <- (ci[1, 2] - ci[1, 1]) / (2 * 1.96)
      result$sigma2_SE <- 2 * theta_est * theta_se
    }
    result$converged <- if (inherits(fit, "merMod")) fit@optinfo$conv$opt == 0 else TRUE
  } else {
    coef_summary <- summary(fit)$coefficients
    result <- list(alpha_est = coef_summary["(Intercept)", "Estimate"],
                   alpha_SE  = coef_summary["(Intercept)", "Std. Error"],
                   beta_est  = coef_summary["X1", "Estimate"],
                   beta_SE   = coef_summary["X1", "Std. Error"],
                   converged = fit$converged)
  }
  
  return(result)
}

## SDM.2 #######################################################################

fit_SDM.2 <- function(Y, X1, model_type, use_OLRE) {
  
  obs_id <- factor(seq_along(Y))
  data_df <- data.frame(Y = Y, X1 = X1, obs_id = obs_id)
  
  family_spec <- get_family(model_type)
  
  if (use_OLRE) {
    fit <- glmmTMB(Y ~ X1 + (1 | obs_id), family = family_spec, data = data_df)
  } else {
    fit <- glmmTMB(Y ~ X1, family = family_spec, data = data_df)
  }
  
  fit_summary <- summary(fit)
  result <- list(alpha_est = fit_summary$coefficients$cond[1, 1],
                 alpha_SE  = fit_summary$coefficients$cond[1, 2],
                 beta_est  = fit_summary$coefficients$cond[2, 1],
                 beta_SE   = fit_summary$coefficients$cond[2, 2])
  
  if (use_OLRE) {
    result$sigma2_est <- fit_summary$varcor$cond$obs_id[1]
    ci <- confint(fit, parm = "theta_", method = "Wald")
    if (nrow(ci) > 0) {
      log_sd_est <- log(sqrt(result$sigma2_est))
      log_sd_se  <- (ci[1, 2] - ci[1, 1]) / (2 * 1.96)
      result$sigma2_SE <- result$sigma2_est * 2 * log_sd_se
    } else {
      result$sigma2_SE <- NA
    }
  }
  
  result$converged <- fit$sdr$pdHess
  
  return(result)
}

## SDM.3 #######################################################################

fit_SDM.3 <- function(Y, X1, model_type, use_OLRE) {
  
  obs_id <- factor(seq_along(Y))
  data_df <- data.frame(Y = Y, X1 = X1, obs_id = obs_id)
  
  get_family_spaMM <- function(model_type) {
    switch(model_type,
           "gaussian"          = gaussian(),
           "poisson"           = poisson(),
           "bernoulli_probit"  = binomial(link = "probit"),
           "bernoulli_logit"   = binomial(link = "logit"),
           "bernoulli_cloglog" = binomial(link = "cloglog"))
  }
  
  family_spec <- get_family_spaMM(model_type)
  
  if (use_OLRE) {
    fit <- fitme(Y ~ X1 + (1 | obs_id), family = family_spec, data = data_df)
    
    fixed_coefs <- fixef(fit)
    fit_summary <- summary(fit, verbose = FALSE)
    se_values   <- fit_summary$beta_table[, "Cond. SE"]
    
    result <- list(alpha_est = fixed_coefs[1],
                   alpha_SE  = se_values[1],
                   beta_est  = fixed_coefs[2],
                   beta_SE   = se_values[2],
                   converged = NA_character_)
    
    lambda <- VarCorr(fit)
    if (!is.null(lambda) && length(lambda) > 0) {
      if ("obs_id" %in% names(lambda)) {
        result$sigma2_est <- as.numeric(lambda[["obs_id"]])
      } else if (length(lambda) > 0) {
        result$sigma2_est <- as.numeric(lambda[[1]])
      }
      ci <- confint(fit, parm = "lambda", verbose = FALSE)
      if (!is.null(ci) && is.matrix(ci) && nrow(ci) > 0) {
        result$sigma2_SE <- (ci[1, 2] - ci[1, 1]) / (2 * 1.96)
      } else {
        result$sigma2_SE <- NA
      }
    }
    
  } else {
    fit <- fitme(Y ~ X1, family = family_spec, data = data_df)
    
    fixed_coefs <- fixef(fit)
    fit_summary <- summary(fit, verbose = FALSE)
    se_values   <- fit_summary$beta_table[, "Cond. SE"]
    
    result <- list(alpha_est = fixed_coefs[1],
                   alpha_SE  = se_values[1],
                   beta_est  = fixed_coefs[2],
                   beta_SE   = se_values[2],
                   converged = NA_character_)
  }
  
  return(result)
}

## SDM.4 #######################################################################

fit_SDM.4 <- function(Y, X1, model_type) {
  
  obs_id <- factor(seq_along(Y))
  data_df <- data.frame(Y = Y, X1 = X1, obs_id = obs_id)
  
  get_family_GLMMadaptive <- function(model_type) {
    switch(model_type,
           "gaussian"         = gaussian(),
           "poisson"          = poisson(),
           "bernoulli_probit" = binomial(link = "probit"),
           "bernoulli_logit"  = binomial(link = "logit"))
  }
  
  family_spec <- get_family_GLMMadaptive(model_type)
  
  fit <- mixed_model(fixed  = Y ~ X1,
                     random = ~ 1 | obs_id,
                     data   = data_df,
                     family = family_spec)
  
  fixed_effects <- fixef(fit)
  se_fixed      <- sqrt(diag(vcov(fit)))
  
  result <- list(alpha_est = fixed_effects[1],
                 alpha_SE  = se_fixed[1],
                 beta_est  = fixed_effects[2],
                 beta_SE   = se_fixed[2],
                 converged = fit$converged)
  
  result$sigma2_est <- fit$D[1, 1]
  vcov_all <- vcov(fit, parm = "var-cov")
  result$sigma2_SE <- sqrt(vcov_all[1, 1])
  
  return(result)
}

## SDM.5 #######################################################################

fit_SDM.5 <- function(Y, X1, model_type, use_OLRE) {
  
  obs_id <- factor(seq_along(Y))
  data_df <- data.frame(Y = Y, X1 = X1, obs_id = obs_id)
  
  get_family_brms <- function(model_type) {
    switch(model_type,
           "gaussian"          = gaussian(),
           "poisson"           = poisson(link = "log"),
           "bernoulli_probit"  = bernoulli(link = "probit"),
           "bernoulli_logit"   = bernoulli(link = "logit"),
           "bernoulli_cloglog" = bernoulli(link = "cloglog"))
  }
  
  family_spec <- get_family_brms(model_type)
  
  if (use_OLRE) {
    formula_brms <- bf(Y ~ X1 + (1 | obs_id))
  } else {
    formula_brms <- bf(Y ~ X1)
  }
  
  prior_spec <- c(prior(normal(0, 10), class = Intercept),
                  prior(normal(0, 10), class = b))
  
  if (use_OLRE) {
    prior_spec <- c(prior_spec, prior(cauchy(0, 1), class = sd))
  }
  
  fit <- brm(formula = formula_brms,
             data    = data_df,
             family  = family_spec,
             prior   = prior_spec,
             chains  = 3,
             iter    = 40000,
             warmup  = 20000,
             thin    = 5,
             cores   = 3,
             control = list(adapt_delta = 0.95, max_treedepth = 12),
             seed    = seed)
  
  post_summary <- posterior_summary(fit, pars = c("b_Intercept", "b_X1"))
  
  result <- list(alpha_est = post_summary["b_Intercept", "Estimate"],
                 alpha_SE  = post_summary["b_Intercept", "Est.Error"],
                 beta_est  = post_summary["b_X1", "Estimate"],
                 beta_SE   = post_summary["b_X1", "Est.Error"],
                 converged = NA_character_)
  
  rhat_vals <- rhat(fit)
  if (any(rhat_vals > 1.1, na.rm = TRUE)) {
    result$converged <- FALSE
  }
  
  if (use_OLRE) {
    re_summary <- VarCorr(fit, summary = TRUE)
    sd_site <- re_summary$obs_id$sd[1, "Estimate"]
    result$sigma2_est <- sd_site^2
    sd_se <- re_summary$obs_id$sd[1, "Est.Error"]
    result$sigma2_SE <- 2 * sd_site * sd_se
  }
  
  return(result)
}

## SDM.6 #######################################################################

fit_SDM.6 <- function(Y, X1, model_type, use_OLRE = TRUE) {
  
  obs_id <- seq_along(Y)
  data_df <- data.frame(Y = Y, X1 = X1, obs_id = obs_id)
  
  get_family_inla <- function(model_type) {
    switch(model_type,
           "gaussian"          = "gaussian",
           "poisson"           = "poisson",
           "bernoulli_probit"  = "binomial",
           "bernoulli_logit"   = "binomial",
           "bernoulli_cloglog" = "binomial",
           NULL)
  }
  
  get_link_inla <- function(model_type) {
    switch(model_type,
           "gaussian"          = "identity",
           "poisson"           = "log",
           "bernoulli_probit"  = "probit",
           "bernoulli_logit"   = "logit",
           "bernoulli_cloglog" = "cloglog",
           "identity")
  }
  
  family_spec <- get_family_inla(model_type)
  link_spec   <- get_link_inla(model_type)
  
  if (use_OLRE) {
    formula_inla <- Y ~ X1 + f(obs_id, model = "iid")
  } else {
    formula_inla <- Y ~ X1
  }
  
  fit <- inla(formula           = formula_inla,
              family            = family_spec,
              data              = data_df,
              control.family    = list(link = link_spec),
              control.compute   = list(config = TRUE),
              control.predictor = list(compute = TRUE),
              verbose           = FALSE)
  
  fixed_effects <- fit$summary.fixed
  
  result <- list(alpha_est = fixed_effects["(Intercept)", "mean"],
                 alpha_SE  = fixed_effects["(Intercept)", "sd"],
                 beta_est  = fixed_effects["X1", "mean"],
                 beta_SE   = fixed_effects["X1", "sd"],
                 converged = TRUE)
  
  if (!is.null(fit$mode$mode.status) && fit$mode$mode.status != 0) {
    result$converged <- FALSE
  }
  
  if (use_OLRE && !is.null(fit$summary.hyperpar)) {
    hyper_summary     <- fit$summary.hyperpar
    precision_row     <- grep("Precision for obs_id", rownames(hyper_summary))
    precision_mean    <- hyper_summary[precision_row, "mean"]
    precision_sd      <- hyper_summary[precision_row, "sd"]
    result$sigma2_est <- 1 / precision_mean
    result$sigma2_SE  <- precision_sd / (precision_mean^2)
  }
  
  return(result)
}

## SDM.7 #######################################################################

fit_SDM.7 <- function(Y, X1, model_type, use_OLRE) {
  
  n_obs  <- length(Y)
  obs_id <- seq_len(n_obs)
  
  if (model_type == "gaussian") {
    nimble_code <- nimbleCode({
      alpha ~ dnorm(0, sd = 10)
      beta ~ dnorm(0, sd = 10)
      tau ~ dgamma(0.001, 0.001)
      sigma <- 1 / sqrt(tau)
      sigma2 <- sigma^2
      for(i in 1:n_obs) {
        mu[i] <- alpha + beta * X1[i]
        Y[i] ~ dnorm(mu[i], sd = sigma)
      }
    })
  } else if (model_type %in% c("bernoulli_probit")) {
    if (use_OLRE) {
      nimble_code <- nimbleCode({
        alpha ~ dnorm(0, sd = 10)
        beta ~ dnorm(0, sd = 10)
        tau ~ dgamma(0.001, 0.001)
        sigma <- 1 / sqrt(tau)
        sigma2 <- sigma^2
        for(j in 1:n_obs) {
          site_effect[j] ~ dnorm(0, sd = sigma)
        }
        for(i in 1:n_obs) {
          probit(p[i]) <- alpha + beta * X1[i] + site_effect[obs_id[i]]
          Y[i] ~ dbern(p[i])
        }
      })
    } else {
      nimble_code <- nimbleCode({
        alpha ~ dnorm(0, sd = 10)
        beta ~ dnorm(0, sd = 10)
        for(i in 1:n_obs) {
          probit(p[i]) <- alpha + beta * X1[i]
          Y[i] ~ dbern(p[i])
        }
      })
    }
  } else if (model_type %in% c("bernoulli_logit")) {
    if (use_OLRE) {
      nimble_code <- nimbleCode({
        alpha ~ dnorm(0, sd = 10)
        beta ~ dnorm(0, sd = 10)
        tau ~ dgamma(0.001, 0.001)
        sigma <- 1 / sqrt(tau)
        sigma2 <- sigma^2
        for(j in 1:n_obs) {
          site_effect[j] ~ dnorm(0, sd = sigma)
        }
        for(i in 1:n_obs) {
          logit(p[i]) <- alpha + beta * X1[i] + site_effect[obs_id[i]]
          Y[i] ~ dbern(p[i])
        }
      })
    } else {
      nimble_code <- nimbleCode({
        alpha ~ dnorm(0, sd = 10)
        beta ~ dnorm(0, sd = 10)
        for(i in 1:n_obs) {
          logit(p[i]) <- alpha + beta * X1[i]
          Y[i] ~ dbern(p[i])
        }
      })
    }
  } else if (model_type %in% c("bernoulli_cloglog")) {
    if (use_OLRE) {
      nimble_code <- nimbleCode({
        alpha ~ dnorm(0, sd = 10)
        beta ~ dnorm(0, sd = 10)
        tau ~ dgamma(0.001, 0.001)
        sigma <- 1 / sqrt(tau)
        sigma2 <- sigma^2
        for(j in 1:n_obs) {
          site_effect[j] ~ dnorm(0, sd = sigma)
        }
        for(i in 1:n_obs) {
          p[i] <- 1 - exp(-exp(alpha + beta * X1[i] + site_effect[obs_id[i]]))
          Y[i] ~ dbern(p[i])
        }
      })
    } else {
      nimble_code <- nimbleCode({
        alpha ~ dnorm(0, sd = 10)
        beta ~ dnorm(0, sd = 10)
        for(i in 1:n_obs) {
          p[i] <- 1 - exp(-exp(alpha + beta * X1[i]))
          Y[i] ~ dbern(p[i])
        }
      })
    }
  } else if (model_type %in% c("poisson")) {
    if (use_OLRE) {
      nimble_code <- nimbleCode({
        alpha ~ dnorm(0, sd = 10)
        beta ~ dnorm(0, sd = 10)
        tau ~ dgamma(0.001, 0.001)
        sigma <- 1 / sqrt(tau)
        sigma2 <- sigma^2
        for(j in 1:n_obs) {
          site_effect[j] ~ dnorm(0, sd = sigma)
        }
        for(i in 1:n_obs) {
          log(lambda[i]) <- alpha + beta * X1[i] + site_effect[obs_id[i]]
          Y[i] ~ dpois(lambda[i])
        }
      })
    } else {
      nimble_code <- nimbleCode({
        alpha ~ dnorm(0, sd = 10)
        beta ~ dnorm(0, sd = 10)
        for(i in 1:n_obs) {
          log(lambda[i]) <- alpha + beta * X1[i]
          Y[i] ~ dpois(lambda[i])
        }
      })
    }
  }
  
  nimble_data <- list(Y      = Y,
                      X1     = X1,
                      obs_id = obs_id)
  
  nimble_constants <- list(n_obs = n_obs)
  
  nimble_inits <- function(chain_id) {
    inits <- list(alpha = rnorm(1, 0, 1),
                  beta  = rnorm(1, 0, 1))
    if (model_type == "gaussian") {
      inits$sigma <- runif(1, 0.1, 2)
    }
    if (use_OLRE && model_type != "gaussian") {
      inits$sigma       <- runif(1, 0.1, 2)
      inits$site_effect <- rnorm(n_obs, 0, 0.1)
    }
    return(inits)
  }
  
  nimble_params <- c("alpha", "beta")
  
  if (use_OLRE || model_type == "gaussian") {
    nimble_params <- c(nimble_params, "sigma")
  }
  
  Nchains <- 3
  
  nimble_out <- runMCMC_btadjust(code         = nimble_code,
                                 constants    = nimble_constants,
                                 data         = nimble_data,
                                 inits        = lapply(1:Nchains, nimble_inits),
                                 params       = nimble_params,
                                 niter.min    = 10000,
                                 niter.max    = Inf,
                                 nburnin.min  = 10000,
                                 nburnin.max  = Inf,
                                 thin.min     = 1,
                                 thin.max     = Inf,
                                 Nchains      = Nchains,
                                 conv.max     = 1.05,
                                 neff.min     = 5000,
                                 control      = list(time.max                   = 36000,
                                                     round.thinmult             = TRUE,
                                                     print.diagnostics          = TRUE,
                                                     Ncycles.target             = 2,
                                                     check.convergence.firstrun = TRUE,
                                                     convtype                   = 'Gelman'),
                                 control.MCMC = list(parallelize = TRUE))
  
  samples <- do.call(rbind, nimble_out)
  attrs   <- attributes(nimble_out)
  
  result <- list(alpha_est = mean(samples[, "alpha"]),
                 alpha_SE  = sd(samples[, "alpha"]),
                 beta_est  = mean(samples[, "beta"]),
                 beta_SE   = sd(samples[, "beta"]),
                 converged = attrs$final.params$converged)
  
  if (use_OLRE) {
    result$sigma2_est <- mean(samples[, "sigma"]^2)
    result$sigma2_SE  <- sd(samples[, "sigma"]^2)
  }
  
  return(result)
}

## SDM.8 #######################################################################

fit_SDM.8 <- function(Y, X1, model_type, use_OLRE) {
  
  n_obs  <- length(Y)
  obs_id <- seq_len(n_obs)
  
  jags_models <- list()
  
  jags_models$gaussian_no_OLRE <- "
    model {
      alpha ~ dnorm(0, 0.01)
      beta ~ dnorm(0, 0.01)
      tau ~ dgamma(0.001, 0.001)
      sigma <- 1 / sqrt(tau)
      sigma2 <- sigma^2
      for (i in 1:n_obs) {
        mu[i] <- alpha + beta * X1[i]
        Y[i] ~ dnorm(mu[i], tau)
      }
    }"
  
  jags_models$bernoulli_probit_no_OLRE <- "
    model {
      alpha ~ dnorm(0, 0.01)
      beta ~ dnorm(0, 0.01)
      for (i in 1:n_obs) {
        z[i] <- alpha + beta * X1[i]
        p[i] <- phi(z[i])
        Y[i] ~ dbern(p[i])
      }
    }"
  
  jags_models$bernoulli_probit_OLRE <- "
    model {
      alpha ~ dnorm(0, 0.01)
      beta ~ dnorm(0, 0.01)
      tau ~ dgamma(0.001, 0.001)
      sigma <- 1 / sqrt(tau)
      sigma2 <- sigma^2
      for (j in 1:n_obs) {
        obs_effect[j] ~ dnorm(0, tau)
      }
      for (i in 1:n_obs) {
        z[i] <- alpha + beta * X1[i] + obs_effect[obs_id[i]]
        p[i] <- phi(z[i])
        Y[i] ~ dbern(p[i])
      }
    }"
  
  jags_models$bernoulli_logit_no_OLRE <- "
    model {
      alpha ~ dnorm(0, 0.01)
      beta ~ dnorm(0, 0.01)
      for (i in 1:n_obs) {
        logit(p[i]) <- alpha + beta * X1[i]
        Y[i] ~ dbern(p[i])
      }
    }"
  
  jags_models$bernoulli_logit_OLRE <- "
    model {
      alpha ~ dnorm(0, 0.01)
      beta ~ dnorm(0, 0.01)
      tau ~ dgamma(0.001, 0.001)
      sigma <- 1 / sqrt(tau)
      sigma2 <- sigma^2
      for (j in 1:n_obs) {
        obs_effect[j] ~ dnorm(0, tau)
      }
      for (i in 1:n_obs) {
        logit(p[i]) <- alpha + beta * X1[i] + obs_effect[obs_id[i]]
        Y[i] ~ dbern(p[i])
      }
    }"
  
  jags_models$bernoulli_cloglog_no_OLRE <- "
    model {
      alpha ~ dnorm(0, 0.01)
      beta ~ dnorm(0, 0.01)
      for (i in 1:n_obs) {
        p[i] <- 1 - exp(-exp(alpha + beta * X1[i]))
        Y[i] ~ dbern(p[i])
      }
    }"
  
  jags_models$bernoulli_cloglog_OLRE <- "
    model {
      alpha ~ dnorm(0, 0.01)
      beta ~ dnorm(0, 0.01)
      tau ~ dgamma(0.001, 0.001)
      sigma <- 1 / sqrt(tau)
      sigma2 <- sigma^2
      for (j in 1:n_obs) {
        obs_effect[j] ~ dnorm(0, tau)
      }
      for (i in 1:n_obs) {
        p[i] <- 1 - exp(-exp(alpha + beta * X1[i] + obs_effect[obs_id[i]]))
        Y[i] ~ dbern(p[i])
      }
    }"
  
  jags_models$poisson_no_OLRE <- "
    model {
      alpha ~ dnorm(0, 0.01)
      beta  ~ dnorm(0, 0.01)
      for (i in 1:n_obs) {
        log(lambda[i]) <- alpha + beta * X1[i]
        Y[i] ~ dpois(lambda[i])
      }
    }"
  
  jags_models$poisson_OLRE <- "
    model {
      alpha ~ dnorm(0, 0.01)
      beta  ~ dnorm(0, 0.01)
      tau ~ dgamma(0.001, 0.001)
      sigma <- 1 / sqrt(tau)
      sigma2 <- sigma^2
      for (j in 1:n_obs) {
        obs_effect[j] ~ dnorm(0, tau)
      }
      for (i in 1:n_obs) {
        log(lambda[i]) <- alpha + beta * X1[i] + obs_effect[obs_id[i]]
        Y[i] ~ dpois(lambda[i])
      }
    }"
  
  model_key         <- paste0(model_type, ifelse(use_OLRE, "_OLRE", "_no_OLRE"))
  jags_model_string <- jags_models[[model_key]]
  
  jags_model <- tempfile(fileext = ".txt")
  writeLines(jags_model_string, jags_model)
  
  jags_data <- list(Y      = as.numeric(Y),
                    X1     = as.numeric(X1),
                    obs_id = as.integer(obs_id),
                    n_obs  = as.integer(n_obs))
  
  jags_params <- c("alpha", "beta")
  
  if (use_OLRE || model_type == "gaussian") {
    jags_params <- c(jags_params, "sigma2")
  }
  
  make_inits <- function(chain_id) {
    inits <- list(alpha = rnorm(1, 0, 1),
                  beta  = rnorm(1, 0, 1))
    if (model_type == "gaussian") {
      sd0       <- runif(1, 0.5, 2)
      inits$tau <- 1 / (sd0^2)
    }
    if (use_OLRE) {
      sd_OLRE           <- runif(1, 0.1, 1.0)
      inits$tau         <- 1 / (sd_OLRE^2)
      inits$site_effect <- rnorm(n_obs, 0, sd = sd_OLRE * 0.5)
    }
    return(inits)
  }
  
  Nchains <- 3
  
  jags_out <- runMCMC_btadjust(MCMC_language = "SDM.8",
                               code          = jags_model,
                               data          = jags_data,
                               inits         = lapply(1:Nchains, make_inits),
                               params        = jags_params,
                               niter.min     = 10000,
                               niter.max     = Inf,
                               nburnin.min   = 10000,
                               nburnin.max   = Inf,
                               thin.min      = 1,
                               thin.max      = Inf,
                               Nchains       = Nchains,
                               conv.max      = 1.05,
                               neff.min      = 5000,
                               control       = list(time.max                   = 36000,
                                                    round.thinmult             = TRUE,
                                                    print.diagnostics          = TRUE,
                                                    Ncycles.target             = 2,
                                                    check.convergence.firstrun = TRUE,
                                                    convtype                   = 'Gelman'),
                               control.MCMC  = list(parallelize = TRUE))
  
  combined_samples <- do.call(rbind, jags_out)
  attrs            <- attributes(jags_out)
  
  result <- list(alpha_est = mean(combined_samples[, "alpha"]),
                 alpha_SE  = sd(combined_samples[, "alpha"]),
                 beta_est  = mean(combined_samples[, "beta"]),
                 beta_SE   = sd(combined_samples[, "beta"]),
                 converged = attrs$final.params$converged)
  
  if ("sigma2" %in% colnames(combined_samples)) {
    result$sigma2_est = mean(combined_samples[, "sigma2"])
    result$sigma2_SE  = sd(combined_samples[, "sigma2"])
  }
  
  unlink(jags_model)
  
  return(result)
}

## SDM.9 #######################################################################

fit_SDM.9 <- function(Y, X1, model_type) {
  
  data_df <- data.frame(Y = Y, X1 = X1)
  
  if (model_type == "gaussian") {
    fit <- lmrob(Y ~ X1, data = data_df)
  } else {
    family_spec <- get_family(model_type)
    fit <- glmrob(Y ~ X1, family = family_spec, data = data_df, method = "Mqle")
  }
  
  coef_summary <- summary(fit)$coefficients
  
  list(alpha_est = coef_summary["(Intercept)", "Estimate"],
       alpha_SE  = coef_summary["(Intercept)", "Std. Error"],
       beta_est  = coef_summary["X1", "Estimate"],
       beta_SE   = coef_summary["X1", "Std. Error"],
       converged = fit$converged)
}

## SDM.10 ######################################################################

fit_SDM.10 <- function(Y, X1, model_type) {
  
  data_df <- data.frame(Y = Y, X1 = X1)
  
  family_spec <- if (model_type == "poisson") quasipoisson() else quasibinomial()
  
  fit <- glm(Y ~ X1, family = family_spec, data = data_df)
  
  coef_summary <- summary(fit)$coefficients
  
  list(alpha_est = coef_summary["(Intercept)", "Estimate"],
       alpha_SE  = coef_summary["(Intercept)", "Std. Error"],
       beta_est  = coef_summary["X1", "Estimate"],
       beta_SE   = coef_summary["X1", "Std. Error"],
       converged = fit$converged)
}

## SDM.11 ######################################################################

fit_SDM.11 <- function(Y, X1, model_type, use_OLRE) {
  
  obs_id <- factor(seq_along(Y))
  data_df <- data.frame(Y = Y, X1 = X1, obs_id = obs_id)
  
  get_family_mgcv <- function(model_type) {
    switch(model_type,
           "gaussian"         = gaussian(),
           "poisson"          = poisson(),
           "bernoulli_probit" = binomial(link = "probit"),
           "bernoulli_logit"  = binomial(link = "logit"))
  }
  
  family_spec    <- get_family_mgcv(model_type)
  smooth_formula <- Y ~ s(X1, bs = "tp", k = 10)
  
  compute_ame <- function(gam_obj, newdata) {
    d_plus     <- newdata
    d_plus$X1  <- d_plus$X1 + 1e-5
    d_minus    <- newdata
    d_minus$X1 <- d_minus$X1 - 1e-5
    Xp_plus    <- predict(gam_obj, newdata = d_plus,  type = "lpmatrix")
    Xp_minus   <- predict(gam_obj, newdata = d_minus, type = "lpmatrix")
    Xd         <- (Xp_plus - Xp_minus) / (2 * 1e-5)
    grad       <- colMeans(Xd)
    ame        <- sum(grad * coef(gam_obj))
    ame_se     <- sqrt(as.numeric(t(grad) %*% vcov(gam_obj) %*% grad))
    list(est = ame, SE = ame_se)
  }
  
  if (use_OLRE) {
    fit <- gamm(formula = smooth_formula,
                random  = list(obs_id = ~ 1),
                family  = family_spec,
                data    = data_df)
    
    gam_summary <- summary(fit$gam)
    ame         <- compute_ame(fit$gam, data_df)
    
    result <- list(alpha_est = unname(coef(fit$gam)[1]),
                   alpha_SE  = gam_summary$se[1],
                   beta_est  = ame$est,
                   beta_SE   = ame$SE,
                   edf_X1    = unname(gam_summary$edf[1]),
                   converged = TRUE)
    
    vc                <- VarCorr(fit$lme)
    result$sigma2_est <- as.numeric(vc["obs_id", "Variance"])
    int               <- intervals(fit$lme, which = "var-cov")
    var_int           <- int$reStruct$obs_id
    result$sigma2_SE. <- (var_int[1, 3] - var_int[1, 1]) / (2 * 1.96)
    
  } else {
    fit <- gam(formula = smooth_formula,
               family  = family_spec,
               data    = data_df,
               method  = "REML")
    
    gam_summary <- summary(fit)
    ame         <- compute_ame(fit, data_df)
    
    result <- list(alpha_est = unname(coef(fit)[1]),
                   alpha_SE  = gam_summary$se[1],
                   beta_est  = ame$est,
                   beta_SE   = ame$SE,
                   edf_X1    = unname(gam_summary$edf[1]),
                   converged = fit$converged)
  }
  
  return(result)
}

## DISPATCHER ##################################################################

fit_model <- function(Y, X1, model_type, method) {
  
  result <- list(method     = method,
                 sdm_label  = NA,
                 olre       = NA,
                 converged  = NA,
                 alpha_est  = NA,
                 alpha_SE   = NA,
                 beta_est   = NA,
                 beta_SE    = NA,
                 sigma2_est = NA,
                 sigma2_SE  = NA,
                 error_msg  = NA)
  
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
  
  sub_result <- switch(base_method,
                       "SDM.1"     = fit_SDM.1(Y, X1, model_type, use_OLRE),
                       "SDM.2"     = fit_SDM.2(Y, X1, model_type, use_OLRE),
                       "SDM.3"     = fit_SDM.3(Y, X1, model_type, use_OLRE),
                       "SDM.4"     = fit_SDM.4(Y, X1, model_type),
                       "SDM.5"     = fit_SDM.5(Y, X1, model_type, use_OLRE),
                       "SDM.6"     = fit_SDM.6(Y, X1, model_type, use_OLRE),
                       "SDM.7"     = fit_SDM.7(Y, X1, model_type, use_OLRE),
                       "SDM.8"     = fit_SDM.8(Y, X1, model_type, use_OLRE),
                       "SDM.9"     = fit_SDM.9(Y, X1, model_type),
                       "SDM.10-QP" = fit_SDM.10(Y, X1, model_type),
                       "SDM.10-QB" = fit_SDM.10(Y, X1, model_type),
                       "SDM.11"    = fit_SDM.11(Y, X1, model_type, use_OLRE))
  
  result           <- modifyList(result, sub_result)
  result$sdm_label <- unname(method_to_sdm[base_method])
  result$olre      <- use_OLRE
  
  return(result)
}

## SINGLE REPLICATE FUNCTION ###################################################

run_single_replicate <- function(seed) {
  
  set.seed(seed)
  
  all_results <- list()
  
  covariates <- generate_covariates(n_obs, seed)
  X1         <- covariates$X1
  X2         <- covariates$X2
  
  for (data_model in data_models) {
    
    Y <- generate_data(n_obs,
                       X1,
                       X2,
                       data_model,
                       seed + which(data_models == data_model))
    
    model_results <- list(data_model  = data_model,
                          true_params = list(alpha       = alpha_true,
                                             beta        = beta_true,
                                             gamma       = gamma_true,
                                             sigma_error = sigma_error),
                          fits = list())
    
    for (method in fit_methods) {
      if (!is_method_compatible(method, data_model)) {
        next
      }
      
      start_time <- Sys.time()
      
      fit_result <- fit_model(Y, X1, data_model, method)
      
      end_time   <- Sys.time()
      
      fit_result$fit_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
      
      model_results$fits[[method]] <- fit_result
    }
    
    all_results[[data_model]] <- model_results
  }
  
  results <- list(seed       = seed,
                  timestamp  = Sys.time(),
                  parameters = list(n_obs     = n_obs,
                                    sigma_error = sigma_error),
                  results    = all_results)
  
  return(results)
}

## RUN #########################################################################

start_time <- Sys.time()

results <- run_single_replicate(seed_value)

end_time <- Sys.time()

results$total_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

output_file <- file.path(output_dir, sprintf("replicate_%06d.RData", seed_value))

save(results, file = output_file)

################################################################################
