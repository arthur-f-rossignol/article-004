################################################################################
##                                                                            ##
##                                    BIAS                                    ##
##                        (diagnostics_metrics/bias.R)                        ##
##                                                                            ##
##       Addressing Missing Covariates in Species Distribution Models:        ##
##  Inferential Impacts and Mitigation via Joint Species Distribution Models  ##
##                                                                            ##
##                  Arthur F. Rossignol & Frédéric Gosselin                   ##
##                                                                            ##
##                                    2026                                    ##
##                                                                            ##
################################################################################

compute_bias <- function(estimates, true_value, data_model) {

  if (length(estimates) < 2 || all(is.na(estimates))) {
    return(list(value        = NA,
                significance = "",
                magnitude    = ""))
  }

  pop <- estimates - true_value
  pop <- pop[is.finite(pop)]

  value_mean_bias <- mean(pop)

  value_median_bias <- median(pop)

  test <- summary(lmrob(pop ~ 1,
                        data = as.data.frame(pop)))$coefficients[4]

  if (test <= 0.001) {
    significance <- "***"
  } else if (test <= 0.01) {
    significance <- "**"
  } else if (test <= 0.05) {
    significance <- "*"
  } else {
    significance <- ""
  }
  
  magnitude_levels  <- c(0.05, 0.1, 0.5, 1)
  n0      <- sum(vapply(magnitude_levels,
                        function(L) mean(pop >= -L & pop <= L) >= 0.95,
                        logical(1)))
  nplus   <- sum(vapply(magnitude_levels,
                        function(L) mean(pop >= L) >= 0.95,
                        logical(1)))
  nminus  <- sum(vapply(magnitude_levels,
                        function(L) mean(pop <= -L) >= 0.95,
                        logical(1)))

  magnitude <- paste0(strrep("0", n0),
                      strrep("+", nplus),
                      paste(rep("-", nminus), collapse = " "))

  return(list(mean_bias    = value_mean_bias,
              median_bias  = value_median_bias,
              significance = significance,
              magnitude    = magnitude))
}

################################################################################
