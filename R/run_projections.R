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

  npyr <- posteriors$proj_years # number of projection years (set by user)
  pyr  <- indat$nyr+1+npyr # actual final projection year
  pyrs <- (indat$nyr+1+1):(indat$nyr+1+npyr) # actual years of projection period

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

    if(round(msyrefpts_long[1],2)==round(msyrefpts$msy[samp],2)&&
       round(msyrefpts_long[2],2)==round(msyrefpts$fmsy[samp],2)&&
       round(msyrefpts_long[3],2)==round(msyrefpts$bmsy[samp],2)&&
       round(msyrefpts_long[4],2)==round(msyrefpts$bo[samp],2)){
      message("Equilibrium reference points match long model!\n")
    }else{
      message("Equilibrium reference points DO NOT match long model!\n")
    } # end ifelse
    print(paste("Long model for posterior sample",samp))
    print(msyrefpts_long)
    print(paste("Equilibrium model for posterior sample",samp))
    print(c(round(msyrefpts$msy[samp],2),round(msyrefpts$fmsy[samp],2),
            round(msyrefpts$bmsy[samp],2),round(msyrefpts$bo[samp],2)))
  } # end if testrefpts

  # Add reference points to proj_out (they are the same for each tac)
  for(i in 1:ntac){
    proj_out[[i]]<- cbind(proj_out[[i]],histrefpts,msyrefpts)
  }

  # Now we need to get stock status estimates
  # (dividing F or Biomass by refpoints or benchmarks)
  # As in project_model(), need to be robust to proj_year >1
  if(npyr==1){
    tmp1 <- data.frame((biomass[nyrs+2]/biomass[nyrs+1]),
                      (biomass[nyrs+2]/bavg),
                      (biomass[nyrs+2]/bmin),
                      (biomass[nyrs+2]/(0.8*bmsy)),
                      (biomass[nyrs+2]/(0.4*bmsy)),
                      (biomass[nyrs+2]/bo),
                      (ft[nyrs+1]/ft[nyrs]),
                      (ft[nyrs+1]/favg),
                      (ft[nyrs+1]/fmsy))
    colnames(output)<-c(paste0("B",nyr+2,"B",nyr+1),
                        paste0("B",nyr+2,"Bavg"),
                        paste0("B",nyr+2,"Bmin"),
                        paste0("B",nyr+2,"08Bmsy"),
                        paste0("B",nyr+2,"04Bmsy"),
                        paste0("B",nyr+2,"B0"),
                        paste0("F",nyr+1,"F",nyr),
                        paste0("F",nyr+1,"Favg"),
                        paste0("F",nyr+1,"Fmsy"))
   proj_out[[i]]<- cbind(proj_out[[i]],tmp1)
  }else{
    tmp1 <- data.frame((biomass[nyrs+2]/biomass[nyrs+1]),
                       (biomass[nyrs+2]/bavg),
                       (biomass[nyrs+2]/bmin),
                       (biomass[nyrs+2]/(0.8*bmsy)),
                       (biomass[nyrs+2]/(0.4*bmsy)),
                       (biomass[nyrs+2]/bo),
                       (ft[nyrs+1]/ft[nyrs]),
                       (ft[nyrs+1]/favg),
                       (ft[nyrs+1]/fmsy))
    cols<-c(paste0("B",nyr+2,"B",nyr+1),
            paste0("B",nyr+2,"Bavg"),
            paste0("B",nyr+2,"Bmin"),
            paste0("B",nyr+2,"08Bmsy"),
            paste0("B",nyr+2,"04Bmsy"),
            paste0("B",nyr+2,"B0"),
            paste0("F",nyr+1,"F",nyr),
            paste0("F",nyr+1,"Favg"),
            paste0("F",nyr+1,"Fmsy"))

    for(ii in 2:npyr){
      tmp2 <- data.frame((biomass[nyrs+1+ii]/biomass[nyrs+1]),
                         (biomass[nyrs+1+ii]/bavg),
                         (biomass[nyrs+1+ii]/bmin),
                         (biomass[nyrs+1+ii]/(0.8*bmsy)),
                         (biomass[nyrs+1+ii]/(0.4*bmsy)),
                         (biomass[nyrs+1+ii]/bo),
                         (ft[nyrs+ii]/ft[nyrs]),
                         (ft[nyrs+ii]/favg),
                         (ft[nyrs+ii]/fmsy))
      tmp1 <- cbind(tmp1,tmp2)
      cols<-c(paste0("B",nyrs+1+ii,"B",nyr+1),
              paste0("B",nyrs+1+ii,"Bavg"),
              paste0("B",nyrs+1+ii,"Bmin"),
              paste0("B",nyrs+1+ii,"08Bmsy"),
              paste0("B",nyrs+1+ii,"04Bmsy"),
              paste0("B",nyrs+1+ii,"B0"),
              paste0("F",nyr+ii,"F",nyr),
              paste0("F",nyr+ii,"Favg"),
              paste0("F",nyr+ii,"Fmsy"))

    } # end for ii
    colnames(tmp3) <- cols
    proj_out[[i]]<- cbind(proj_out[[i]],tmp1)
  }# end ifelse

  # Now we have a lovely list of all the metrics we need for decision tables
  #  one list object per tac :-)
  return(proj_out)
}
