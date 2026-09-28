################################################################################
##                                                                            ##
##                   ROOT MEAN SQUARE RANDOM ERROR (RMSRE)                    ##
##                       (diagnostics_metrics/RMSRE.R)                        ##
##                                                                            ##
##       Addressing Missing Covariates in Species Distribution Models:        ##
##  Inferential Impacts and Mitigation via Joint Species Distribution Models  ##
##                                                                            ##
##                  Arthur F. Rossignol & Frédéric Gosselin                   ##
##                                                                            ##
##                                    2026                                    ##
##                                                                            ##
################################################################################

compute_RMSRE <- function(estimates, ses, true_value) {
  return(sqrt(mean((estimates - true_value)^2 + ses^2, na.rm = TRUE)))
}

################################################################################
