#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GET POSTERIOR DERIVED VARIABLES ARRANGED BY VARIABLE OR BY SAMPLE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get_posterior_derived_variables

# Called from run-model-v1.R

# Author: Robyn Forrest
# Date created:  August 2, 2024
# Last Modified: August 2, 2024

# This function handles getting out the posterior samples of derived model variables
#  i.e., REPORT() objects from model

#  See readme at https://github.com/kaskr/tmbstan:
#    "What if you want a posterior for derived quantities in the report? Just
#    loop through each posterior sample (row) and call the report function"

# Version 1: posteriors_by_variable
# For reporting, graphs etc
# This is a list where each element is a variable of interest
# Make a list for putting posterior estimates of biomass,
#    recruits and q
# Note that parameters, rec devs and logf are
#    already reported in the mc object
#   (but have added Ft to REPORT to simplify projection model)

# Version 2: posteriors_by_sample
#  for projections
# This is a list where each element is a posterior sample
#  containing all the variables needed for the proj model
#  Gets passed to projection model using purrr::map2_df()

# arguments:
# obj = rtmb object from model
# mc = tmbstan output from fitmcmc in run_model_v1.R
# type = sets whether to return posteriors by variable (type=byvariable)
#  or posteriors by sample (type=bysample)
# proj_years = number of projection years

# Returns:
# if type=byvariable:
#   A list with length 6 (biomass, numbers, recruits, surv, Ft, q)
#   Each list item contains a dataframe with nrow=number of posterior samples (minus warmup)
# if type=bysample:
#   A list with length = number of posterior samples (minus warmup)
#   Each list item contains a list with variables needed for projections
#   Also contains posteriors for the objective function components
# Note the burn in (warmup) samples have already been removed

get_posterior_derived_variables <- function(obj,mc,type="bysample", proj_years){

  # Version 1: for reporting, graphs etc
  # posteriors_by_variable
  if(type=="byvariable"){
    # prepare the list
    posteriors_by_variable <- list()
    posteriors_by_variable$biomass  <- matrix(NA, ncol=nyrs+1, nrow=nrow(mcdf))
    posteriors_by_variable$numbers <- matrix(NA, ncol=nyrs+1, nrow=nrow(mcdf))
    posteriors_by_variable$recruits <- matrix(NA, ncol=nyrs-dat$sage, nrow=nrow(mcdf))
    posteriors_by_variable$surv <- matrix(NA, ncol=nyrs, nrow=nrow(mcdf))
    posteriors_by_variable$Ft <- matrix(NA, ncol=nyrs, nrow=nrow(mcdf))
    posteriors_by_variable$q <- matrix(NA, ncol=dat$nit, nrow=nrow(mcdf))

    # loop over posterior samples
    for(i in 1:nrow(mcdf)){
      r <- obj$report(mcdf[i,])

      posteriors_by_variable$biomass[i,]  <- r$biomass
      posteriors_by_variable$numbers[i,]  <- r$numbers
      posteriors_by_variable$recruits[i,] <- r$rt
      posteriors_by_variable$surv[i,]  <- r$surv
      posteriors_by_variable$Ft[i,] <- r$Ft
      posteriors_by_variable$q[i,]  <- r$q
    }# end for i
    output <- posteriors_by_variable
  } # end if

  # Version 2: for projections
  # posteriors_by_sample
  if(type=="bysample"){
    # prepare the list
    posteriors_by_sample <- list()

    # loop over posterior samples
    for(i in 1:nrow(mcdf)){
      r <- obj$report(mcdf[i,])

      posteriors_by_sample[[i]] <- r
      # Append 3 leading parameters and proj_years
      posteriors_by_sample[[i]]$log_ro <- mcdf$log_ro[i]
      posteriors_by_sample[[i]]$h      <- mcdf$h[i]
      posteriors_by_sample[[i]]$log_m  <- mcdf$log_m[i]
      posteriors_by_sample[[i]]$proj_years  <- proj_years # for now add projection years here
    } # end for i
    output <- posteriors_by_sample
  } # end if

 return(output)
}
