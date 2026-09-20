################################################################################
##                                                                            ##
##                       ROOT MEAN SQUARE ERROR (RMSE)                        ##
##                        (diagnostics_metrics/RMSE.R)                        ##
##                                                                            ##
##       Addressing Missing Covariates in Species Distribution Models:        ##
##  Inferential Impacts and Mitigation via Joint Species Distribution Models  ##
##                                                                            ##
##                  Arthur F. Rossignol & Frédéric Gosselin                   ##
##                                                                            ##
##                                    2026                                    ##
##                                                                            ##
################################################################################

compute_RMSE <- function(estimates, true_value) {
  sqrt(mean((estimates - true_value)^2, na.rm = TRUE))
}

################################################################################
