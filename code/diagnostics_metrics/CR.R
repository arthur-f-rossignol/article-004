################################################################################
##                                                                            ##
##                             COVERAGE RATE (CR)                             ##
##                         (diagnostics_metrics/CR.R)                         ##
##                                                                            ##
##       Addressing Missing Covariates in Species Distribution Models:        ##
##  Inferential Impacts and Mitigation via Joint Species Distribution Models  ##
##                                                                            ##
##                  Arthur F. Rossignol & Frédéric Gosselin                   ##
##                                                                            ##
##                                    2026                                    ##
##                                                                            ##
################################################################################

## Single entry point returning the coverage rate and its significance stars
## for a vector of replicate estimates and standard errors.

compute_CR <- function(estimates,
                       SEs,
                       true_value,
                       converged = NULL,
                       z = 1.96,
                       target_CR = 0.05) {

  n <- length(estimates)
  if (is.null(converged)) converged <- rep(TRUE, n)

  CI_lower <- estimates - z * SEs
  CI_upper <- estimates + z * SEs

  ## computation: coverage rate over the converged replicates

  covered <- ifelse(converged & !is.na(SEs),
                    true_value >= CI_lower & true_value <= CI_upper,
                    NA)

  coverage <- mean(covered, na.rm = TRUE)

  ## significance: randomized exact binomial test of the non-coverage rate
  ## against the target error rate

  valid <- is.finite(CI_lower) & is.finite(CI_upper)

  lo <- CI_lower[valid]
  hi <- CI_upper[valid]

  N <- length(lo)
  if (N < 1) {
    return(list(coverage     = coverage,
                significance = ""))
  }

  m.neg <- sum(true_value < lo)
  m.pos <- sum(true_value > hi)

  set.seed(1)
  u <- runif(1)

  crude <- pbinom(m.neg + m.pos - 1, size = N, prob = target_CR) +
    u * dbinom(m.neg + m.pos, size = N, prob = target_CR)
  test <- min(crude, 1 - crude) * 2

  if (test < 0.001) {
    significance <- "***"
  } else if (test < 0.01) {
    significance <- "**"
  } else if (test < 0.05) {
    significance <- "*"
  } else {
    significance <- ""
  }

  list(coverage     = coverage,
       significance = significance)
}

################################################################################
