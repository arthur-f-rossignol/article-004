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
                       converged = NULL) {
  
  n <- length(estimates)
  if (is.null(converged)) {
    converged <- rep(TRUE, n)
  }

  CI_lower <- estimates - 1.96 * SEs
  CI_upper <- estimates + 1.96 * SEs

  covered <- ifelse(converged & !is.na(SEs),
                    true_value >= CI_lower & true_value <= CI_upper,
                    NA)
  
  mean(covered, na.rm = TRUE)
}

## ROOT MEAN SQUARE ERROR (RMSE) ###############################################

compute_RMSE <- function(estimates, 
                         true_value) {
  sqrt(mean((estimates - true_value)^2, na.rm = TRUE))
}

## ROOT MEAN SQUARE RANDOM ERROR (RMSRE) #######################################

compute_RMSRE <- function(estimates, 
                          SEs,
                          true_value) {
  sqrt(mean((estimates - true_value)^2 + SEs^2, na.rm = TRUE))
}

## BIAS SIGNIFICANCE ###########################################################

compute_bias_significance <- function(estimates, 
                                      true_value) {
  
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
    c(0.5, 0.1, 0.05)
  } else if (grepl("poisson", data_model, ignore.case = TRUE)) {
    c(0.5, 0.1, 0.05)
  } else if (grepl("bernoulli", data_model, ignore.case = TRUE)) {
    c(0.5, 0.1, 0.05)
  } else {
    c(0.5, 0.1, 0.05)
  }
}

compute_bias_magnitude <- function(estimates, true_value, data_model) {
  
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

compute_CR_significance <- function(estimates,
                                    SEs, 
                                    true_value) {

  CI_lower <- estimates - 1.96 * SEs
  CI_upper <- estimates + 1.96 * SEs

  m.neg <- sum(true_value < CI_lower)
  m.pos <- sum(true_value > CI_upper)

  set.seed(1)
  u <- runif(1)

  crude <- pbinom(m.neg + m.pos - 1, size = N, prob = target_CR0.05) +
    u * dbinom(m.neg + m.pos, size = N, prob = 0.05)
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
            bias_significance = compute_bias_significance(coef_est, 
                                                          beta_true),
            bias_magnitude = compute_bias_magnitude(coef_est, 
                                                    beta_true, 
                                                    data_model),
            CR_significance = compute_CR_significance(beta_true,
                                                      coef_est,
                                                      coef_se),
            .groups = "drop")

################################################################################
