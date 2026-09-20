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

## Single entry point returning the bias value, its significance stars and
## its magnitude code for a vector of replicate estimates.

compute_bias <- function(estimates, true_value, data_model) {

  if (length(estimates) < 2 || all(is.na(estimates))) {
    return(list(value        = NA,
                significance = "",
                magnitude    = ""))
  }

  pop <- estimates - true_value
  pop <- pop[is.finite(pop)]

  if (length(pop) < 1) {
    return(list(value        = NA,
                significance = "",
                magnitude    = ""))
  }

  ## computation: mean bias over the replicates

  value <- mean(pop)

  ## significance: robust one-sample test of a zero location

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

  ## magnitude: "0" if 95% of biases fall within the closed band [-L, L],
  ## "+" / "-" if 95% exceed it upward / downward; one symbol per level,
  ## minus signs separated by spaces so they do not merge visually

  levels  <- get_magnitude_levels(data_model)
  n0      <- sum(vapply(levels,
                        function(L) mean(pop >= -L & pop <= L) >= 0.95,
                        logical(1)))
  nplus   <- sum(vapply(levels,
                        function(L) mean(pop >= L) >= 0.95,
                        logical(1)))
  nminus  <- sum(vapply(levels,
                        function(L) mean(pop <= -L) >= 0.95,
                        logical(1)))

  magnitude <- paste0(strrep("0", n0),
                      strrep("+", nplus),
                      paste(rep("-", nminus), collapse = " "))

  list(value        = value,
       significance = significance,
       magnitude    = magnitude)
}

## magnitude thresholds per data model

magnitude_levels <- c(0.05, 0.1, 0.5, 1)

get_magnitude_levels <- function(data_model) {
  ## kept for API compatibility: the levels no longer depend on the family
  magnitude_levels
}

################################################################################
