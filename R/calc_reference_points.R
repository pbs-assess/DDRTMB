#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE REFERENCE POINTS FOR THE DELAY DIFFERENCE MODEL
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Author: Robyn Forrest
# Date created:  August 1, 2024
# Last Modified: August 1, 2024

# Functions in this script:
# 1. ddiff_msy based on ddiff_msy in devs/iscam.tpl
# Use equilibrium calculations for delay difference model

# calcReferencePoints (based on iscam calcReferencePoints
# (see devs\iscam.tpl and devs\msy.cpp)

# 2. FOR TESTING ONLY:
#  ddiff_msy_slow
# The slow way, run out the model for a long time under a fixed F
# arguments and check the final values match the ye and be values from
# the first way

# 3. NOT IMPLEMENTED YET:
#  COULD ALSO TRY IMPLEMENTING THE AGE STRUCTURED VERSION WITH KNIFE-EDGED
#  SELECTIVITY AND MATURITY AT KAGE. SHOULD RETURN SAME REF POINTS AS
#  DDIFF_MSY AND *MIGHT BE QUICKER*

# 1. ddiff_msy - delay-difference equilibrium calculations
# arguments:
# posteriors = a list of posterior outputs from one mcmc posterior sample
# alpha_g, rho_g and wk = growth parameters from the dat file

# Returns:
# a list with msy, fmsy and bmsy for that posterior sample

ddiff_msy <- function(posteriors,
                           alpha_g,
                           rho_g,
                           wk){
   ftest <- seq(0,3,by=.001)
   rec_a <- posteriors$alpha.sr
   rec_b <- posteriors$beta.sr
   M <- exp(posteriors$log_m)

   ye <- numeric(length(ftest))
   be <- numeric(length(ftest))

  # Calculate equilibrium survivorship as function of FMSY
  for(k in 1:length(ftest)){
    # se = eqm survival
    # we = eqm average weight
    se <- exp(-M - ftest[k])
    we <- (se*alpha_g+wk *(1.-se))/(1.-rho_g*se)
    be[k] = -1 * ((-we + se * alpha_g + se * rho_g * we +
              wk * rec_a * we)/(rec_b * (-we + se * alpha_g + se * rho_g * we)))
    ye[k] = be[k] * (1.0 - exp(-ftest[k] - M)) * (ftest[k] /(ftest[k] + M))
    if(ye[k] < 0){
      ye[k] = 0
    }
    if(be[k] < 0){
      be[k] = 0
    }
  }# end for k
  msy = max(ye);
  fmsy = ftest[which(ye==msy)]
  bmsy = be[which(ye==msy)]

  out <- list("msy"=msy, "fmsy"=fmsy, "bmsy"=bmsy)

  return(out)
} # end function


# 2. ddiff_msy_slow
# arguments:
# posteriors = a list of posterior outputs from one mcmc posterior sample
# alpha_g, rho_g and wk = growth parameters from the dat file

# Returns:
# a list with msy, fmsy and bmsy for that posterior sample

ddiff_msy_slow() <- function(posteriors,
                      alpha_g,
                      rho_g,
                      wk){

  yrs <- 1:1000 # a long time
  ftest <- seq(0,3,by=.001)
  rec_a <- posteriors$alpha.sr
  rec_b <- posteriors$beta.sr
  M <- exp(posteriors$log_m)

  ye <- numeric(length(ftest))
  be <- numeric(length(ftest))

  # Calculate equilibrium biomass and yield as function of F
  for(k in 1:length(ftest)){

    surv <- numeric(yrs)
    numbers <- numeric(yrs)
    biomass <- numeric(yrs)
    recruits <- numeric(yrs)
    annual_mean_wt <- numeric(yrs)

    # Initial Conditions
    snat <- exp(-M) # natural survival rate
    # Eqm mean weight
    wbar <- (snat*alpha_g + wk*(1-snat))/(1-snat*rho_g) # RF: calculation checked against rep file
    # Unfished numbers and biomass
    no <- ro/(1-snat)
    bo <- no*wbar # RF: calculation checked against rep file
    # initialize at posterior initial conditions
    surv[i] <- exp(-M-ftest[k]) # fished survival rate
    numbers[1] <- posteriors$numbers[1]
    biomass[1] <- posteriors$biomass[1]
    annual_mean_wt[1] <- biomass[1]/numbers[1]
    recruits[1] <- posteriors$recruits[1]

    for(j in 2:yrs){


    } # end j loop
   }# end for k

  msy = max(ye);
  fmsy = ftest[which(ye==msy)]
  bmsy = be[which(ye==msy)]


  out <- list("msy"=msy, "fmsy"=fmsy, "bmsy"=bmsy)

  return(out)
} # end function

