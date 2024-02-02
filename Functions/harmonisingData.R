# ============================================================================ #
#                               harmonisingData                                #
#                            Marta Alcalde-Herraiz                             #
# ---------------------------------------------------------------------------- #
# Summary:                                                                     #
# This scripts harmonises data between exposure and outcome to compare         #
# associations between exposure-positive control                               #
# ============================================================================ #

harmonisingData <- function(exposure,
                            outcome,
                            gwasID = NULL){
  if(nrow(outcome) != nrow(exposure)){
    source(here::here('Functions','findProxies.R'))
    outcome <- findProxies(exposure,outcome, gwasID = gwasID)
  }

  dat_harmonised <- TwoSampleMR::harmonise_data(exposure,outcome)
  return(dat_harmonised)
}
