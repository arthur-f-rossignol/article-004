################################################################################
##                                                                            ##
##                             DIAGNOSTIC METRICS                             ##
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

library(robustbase)
library(dplyr)

## COVERAGE RATE (CR) ##########################################################

compute_CR <- function(estimates,
                       SEs, 
                       true_value, 
                       converged = NULL, 
                       z = 1.96) {
  
  n <- length(estimates)
  if (is.null(converged)) converged <- rep(TRUE, n)

  CI_lower <- estimates - z * SEs
  CI_upper <- estimates + z * SEs

  covered <- ifelse(converged & !is.na(SEs),
                    true_value >= CI_lower & true_value <= CI_upper,
                    NA)
  
  mean(covered, na.rm = TRUE)
}

## ROOT MEAN SQUARE ERROR (RMSE) ###############################################

compute_RMSE <- function(estimates, true_value) {
  sqrt(mean((estimates - true_value)^2, na.rm = TRUE))
}

## ROOT MEAN SQUARE RANDOM ERROR (RMSRE) #######################################

compute_RMSRE <- function(estimates, ses, true_value) {
  sqrt(mean((estimates - true_value)^2 + ses^2, na.rm = TRUE))
}

## BIAS SIGNIFICANCE ###########################################################

compute_bias_significance <- function(estimates, 
                                      true_value) {
  
  if (length(estimates) < 2 || all(is.na(estimates))) {
    return("")
  }

  test <- summary(lmrob(estimates - true_value ~ 1,
                          data = as.data.frame(estimates)))$coefficients[4]

  if (test < 0.001) {
    return("***")
  } else if (test < 0.01) {
    return("**")
  } else if (test < 0.05) {
    return("*")
  } else {
    return("")
  }
}

## BIAS MAGNITUDE ###########################################################

get_magnitude_levels <- function(data_model) {
  if (grepl("gaussian", data_model, ignore.case = TRUE)) {
    c(0.2, 0.1, 0.05)
  } else if (grepl("poisson", data_model, ignore.case = TRUE)) {
    c(0.2, 0.1, 0.05)
  } else if (grepl("bernoulli", data_model, ignore.case = TRUE)) {
    c(0.2, 0.1, 0.05)
  } else {
    c(0.2, 0.1, 0.05)
  }
}

compute_bias_magnitude <- function(estimates, true_value, data_model) {
  
  if (length(estimates) < 2 || all(is.na(estimates))) {
    return("")
  }
  
  pop <- estimates - true_value
  pop <- pop[is.finite(pop)]
  
  if (length(pop) < 1) {
    return("")
  }
  
  levels <- get_magnitude_levels(data_model)
  magnitude <- ""

  for (level in levels) {
    if (mean(pop > -level & pop < level) >= 0.95) {
      magnitude <- paste0(magnitude, "0")
    }
  }
  
  for (level in levels) {
    if (mean(pop > level) >= 0.95) {      
      magnitude <- paste0(magnitude, "+")
    }
  }
  
  for (level in levels) {
    if (mean(pop < -level) >= 0.95) {
      magnitude <- paste0(magnitude, "-")
    }
  }
  
  magnitude
}

## COVERAGE RATE SIGNIFICANCE ##################################################

compute_CR_significance <- function(true_value, 
                                    CI_lower, 
                                    CI_upper,
                                    target_CR = 0.05) {
  
  valid <- is.finite(CI_lower) & is.finite(CI_upper)
  
  lo <- CI_lower[valid]
  hi <- CI_upper[valid]
  
  N  <- length(lo)
  if (N < 1) {
    return("")
  }

  m.neg <- sum(true_value < lo)
  m.pos <- sum(true_value > hi)

  set.seed(1)
  u <- runif(1)

  crude <- pbinom(m.neg + m.pos - 1, size = N, prob = target_CR) +
    u * dbinom(m.neg + m.pos, size = N, prob = target_CR)
  test <- min(crude, 1 - crude) * 2

  if (test < 0.001) {
    "***"
  } else if (test < 0.01) {
    "**"
  } else if (test < 0.05) {
    "*"
  } else {
    ""
  }
}

## USAGE SKETCH FOR A PER-REPLICATE DATA FRAME #################################

d %>%
  group_by(data_model, method) %>%
  summarise(CR = compute_CR(coef_est,
                            coef_se,
                            beta_true,
                            converged),
            RMSE = compute_RMSE(coef_est, 
                                beta_true),
            RMSRE = compute_RMSRE(coef_est, 
                                  coef_se,
                                  beta_true),
            bias_sign = compute_bias_significance(coef_est, 
                                                  beta_true),
            bias_magn = compute_bias_magnitude(coef_est, 
                                               beta_true, 
                                               first(data_model)),
            CR_sign = compute_CR_significance(beta_true,
                                              coef_est - 1.96 * coef_se,
                                              coef_est + 1.96 * coef_se),
            .groups = "drop")

################################################################################