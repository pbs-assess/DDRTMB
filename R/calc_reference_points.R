#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE REFERENCE POINTS FOR THE DELAY DIFFERENCE MODEL
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Author: Robyn Forrest
# Date created:  August 1, 2024
# Last Modified: August 1, 2024

# Functions in this script:
# 1. ddiff_msy based on ddiff_msy in devs/iscam.tpl
# Use equilibrium calculations for delay difference model

# 2. ddiff_msy_long (FOR TESTING ONLY)
# The slow way, run out the model for a long time under a sequence of fixed F
#  and check the final year values match the ye and be values from ddiff_msy

# 3. calc_hist_refpts
# Gets the historical reference points for biomass and F
#  based on the posterior estimates of biomass and F,
#   and settings in the ctl file that set the end years for averaging and
#   the year for bmin
# TODO: Make the settings more flexible

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. ddiff_msy - delay-difference equilibrium calculations
# arguments:
# posteriors = a list of posterior outputs from one mcmc posterior sample

# Returns:
# a df with msy, fmsy, bmsy and bo for that posterior sample
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ddiff_msy <- function(posteriors){

   ftest <- seq(0,3,by=.001)
   nft <- length(ftest)
   rec_a <- posteriors$alpha.sr
   rec_b <- posteriors$beta.sr
   M <- exp(posteriors$log_m)
   ro <- exp(posteriors$log_ro)
   # alpha_g, rho_g and wk = growth parameters from the dat file - assumes dat is in global space
   alpha_g <- dat$alpha.g
   rho_g <- dat$rho.g
   wk <- dat$wk

   ye <- numeric(nft)
   be <- numeric(nft)

  # Calculate equilibrium survivorship as function of FMSY
  for(k in 1:nft){
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

  msy <- max(ye)
  fmsy <- ftest[which(ye==msy)]
  bmsy <- be[which(ye==msy)]

  out <- data.frame("msy"=msy, "fmsy"=fmsy, "bmsy"=bmsy, "bo"=be[1])

  return(out)
} # end function

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. ddiff_msy_slow - just for testing eqm calcs in ddiff_msy
# arguments:
# posteriors = a list of posterior outputs from one mcmc posterior sample

# Returns:
# a df with msy, fmsy, bmsy and bo for that posterior sample
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ddiff_msy_long <- function(posteriors){

  yrs <- 1:100 # a long time
  nyrs <- length(yrs)
  ftest <- seq(0,3,by=.001)
  nft <- length(ftest)
  rec_a <- posteriors$alpha.sr
  rec_b <- posteriors$beta.sr
  M <- exp(posteriors$log_m)
  ro <- exp(posteriors$log_ro)
  # alpha_g, rho_g and wk = growth parameters from the dat file
  # kage = age at recrtuitment from the dat file
  # srr = stock recruit relationship 1=BH, 2=Ricker, from the ctl file
  alpha_g <- dat$alpha.g
  rho_g <- dat$rho.g
  wk <- dat$wk
  kage <- dat$kage
  srr <- ctl$misc[2]

  ye <- numeric(nft)
  be <- numeric(nft)

  # Calculate equilibrium biomass and yield as function of F  1:nft
  for(k in 1:nft){

    numbers <- numeric(nyrs)
    biomass <- numeric(nyrs)
    recruits <- numeric(nyrs)

    # Initial Conditions
    snat <- exp(-M) # natural survival rate
    surv <- exp(-M-ftest[k]) # fished survival rate (constant bc M and ft are constant)
    # Eqm mean weight
    wbar <- (snat*alpha_g + wk*(1-snat))/(1-snat*rho_g) # RF: calculation checked against rep file
    # Unfished numbers and biomass
    no <- ro/(1-snat)
    bo <- no*wbar # RF: calculation checked against rep file

    # initialize at bo
    numbers[1] <- no
    biomass[1] <- bo
    recruits[1] <- ro

    for(i in 2:nyrs){
       if(i <= dat$kage){
         recruits[i] <- ro
       }else{
         sbt <- biomass[i-kage]
         if(srr==1)recruits[i] <- (posteriors$alpha.sr*sbt/(1.+posteriors$beta.sr*sbt))
         if(srr==2)recruits[i] <- (posteriors$alpha.sr*sbt*exp(-posteriors$beta.sr*sbt))
       }

      biomass[i] <- surv*(rho_g*biomass[i-1]+alpha_g*numbers[i-1]) +wk*recruits[i] # eq. 9.2.5 in HW
      numbers[i] <- surv*numbers[i-1]+recruits[i]

    } # end i loop (yrs)

    # get final equilibrium values
    # just calculate yield in final year
    be[k] <- biomass[nyrs] # take last biomass value as eqm value
    ye[k] <- (ftest[k]/(ftest[k]+M)) * (1-exp(-(ftest[k]+M)))*be[k]

  }# end k loop (ftest)

  msy <- max(ye)
  fmsy <- ftest[which(ye==msy)]
  bmsy <- be[which(ye==msy)]

  out <- data.frame("msy"=msy, "fmsy"=fmsy, "bmsy"=bmsy, "bo"=bo)

  return(out)

} # end function

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. calc_hist_refpts
# Gets the historical reference points for biomass and F
# arguments:
# posteriors = a list of posterior outputs from one mcmc posterior sample

# Returns:
# a df with bavg, favg and bmin for that posterior sample
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

calc_hist_refpts <- function(posteriors){

  yrs  <- dat$syr:dat$nyr # actual historical years
  nyrs <- length(yrs) # number of historical years

  biomass <- posteriors$biomass
  ft <- posteriors$Ft

  # Create a lookup for years because ADMB works with actual years as index
  year_lookup <- as.data.frame(cbind(yrs, 1:nyrs))
  colnames(year_lookup) <- c("year", "year_index")

  endyr_refpt <- pfc$ctl.options[7]
  endyr_refpt_ind <- year_lookup[which(year_lookup[,1]==endyr_refpt),2]
  bminyr <- pfc$ctl.options[9]
  bminyr_ind <- year_lookup[which(year_lookup[,1]==bminyr),2]

  bavg <- mean(biomass[1:endyr_refpt_ind])
  favg <- mean(ft[1:endyr_refpt_ind])
  bmin  <- biomass[bminyr_ind]

  out <- data.frame("bavg"=bavg, "favg"=favg, "bmin"=bmin)
  return(out)
}

