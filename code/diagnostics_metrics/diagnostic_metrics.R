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

coverage_rate <- function(estimates,
                          SEs, 
                          true_value, 
                          converged = NULL) {
  
  n <- length(estimates)

  CI_lower <- estimates - 1.96 * SEs
  CI_upper <- estimates + 1.96 * SEs

  covered <- true_value >= CI_lower & true_value <= CI_upper)
  
  mean(covered)
}

## ROOT MEAN SQUARE ERROR (RMSE) ###############################################

RMSE <- function(estimates,
                 true_value) {
  sqrt(mean((estimates - true_value)^2, na.rm = TRUE))
}

## ROOT MEAN SQUARE RANDOM ERROR (RMSRE) #######################################

RMSRE <- function(estimates, 
                  SEs,
                  true_value) {
  sqrt(mean((estimates - true_value)^2 + SEs^2, na.rm = TRUE))
}

## BIAS SIGNIFICANCE ###########################################################

bias_significance <- function(estimates,
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

bias_magnitude <- function(estimates, 
                           true_value, 
                           data_model) {
  
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

coverage_rate_significance <- function(estimates,
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

df %>%
  group_by(data_model, method) %>%
  summarise(CR = coverage_rate(coef_est,
                               coef_se,
                               beta_true),
            RMSE = RMSE(coef_est, 
                        beta_true),
            RMSRE = RMSRE(coef_est,
                          coef_se,
                          beta_true),
            bias_significance = bias_significance(coef_est, 
                                                   beta_true),
            bias_magnitude = bias_magnitude(coef_est,
                                            beta_true, 
                                            data_model),
            CR_significance = coverage_rate_significance(beta_true,
                                                         coef_est,
                                                         coef_se),
            .groups = "drop")

################################################################################
