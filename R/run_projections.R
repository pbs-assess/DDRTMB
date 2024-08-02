#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN PROJECTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run_projections
# Called from run-model-v1.R

# Author: Robyn Forrest
# Date created:  August 2, 2024
# Last Modified: August 2, 2024

# This function handles calling the projection model under the range
#  of future tacs and calls the reference points functions

# arguments:
# posteriors = posterior outputs from the MCMC, includes derived variables

# Returns:
# A list with length ntac (number of projected tacs)
# Each list item contains a dataframe with all the projection variables
#   for each posterior sample
# Note the burn in (warmup) samples have already been removed

# NOTE: as for all the function scripts, this function assumes that dat, ctl and pfc
#  are already in the global space

run_projections <- function(posteriors){

  # List object for projection outputs
  # These will be the inputs for decision tables
  proj_out <- list()

  # Need to loop over future TACs but do not need to loop
  #  over posterior samples. Let purrr do that.
  ntac <- pfc$num.tac

  for(i in 1:ntac){
    tac <- pfc$tac.vec[i]
    message(paste("Getting projections for TAC",tac))

    # Run the projection model for tac[i]
    proj_out[[i]] <- purrr::map2_df(posteriors_by_sample, tac, project_model)
  }
  names(proj_out) <- pfc$tac.vec

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Get reference points
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 1. Historical reference points first
  # Note that the function depends on dat and pfc being in the global space
  message("Calculating reference points\n")
  histrefpts <- purrr::map_df(posteriors_by_sample,calc_hist_refpts)

  # 2. Do MSY-based reference points
  # Note that the function depends on dat and ctl being in the global space
  msyrefpts <- purrr::map_df(posteriors_by_sample,ddiff_msy)

  # testrefpts should be in global space
  if(testrefpts == TRUE){
    # For testing ddiff_msy - run out a delay difference model for 100 years
    #   to make sure eqm calcs in ddiff_msy are actually returning eqm values
    #   just pick one posterior sample (too slow to do all of them!)
    message("Checking reference points against long model\n")
    samp <- 1
    msyrefpts_long <- ddiff_msy_long(posteriors_by_sample[[samp]])

    if(msyrefpts_long[1]==msyrefpts$msy[samp]&&
       msyrefpts_long[2]==msyrefpts$fmsy[samp]&&
       msyrefpts_long[3]==msyrefpts$bmsy[samp]&&
       msyrefpts_long[4]==msyrefpts$bo[samp]){
      message("Equilibrium reference points match long model!\n")
    }else{
      message("Equilibrium reference points DO NOT match long model!\n")
      print("Long model")
      print(msyrefpts_long)
      print("Equilibrium model")
      print(c(msyrefpts$msy[samp],msyrefpts$fmsy[samp],msyrefpts$bmsy[samp],msyrefpts$bo[samp]))
    } # end ifelse
  } # end if

  # Add reference points to proj_out (they are the same for each tac)
  for(i in 1:ntac){
    proj_out[[i]]<- cbind(proj_out[[i]],histrefpts,msyrefpts)
  }

  # Now we have a lovely list of all the metrics we need for decision tables
  #  one list object per tac :-)
  return(proj_out)
}
