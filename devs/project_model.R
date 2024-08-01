#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PROJECTION MODEL
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# projection_model (based on iscam project_model_dd (see devs\iscam.tpl))

# Author: Robyn Forrest
# Date created:  July 26, 2024
# Last Modified: August 1, 2024

# This model takes one sample of posterior parameters and historical estimates
#  and projects biomass and F forward the number of projection years (proj_years)
#  It also calculates ref points based on estimated parameters so that
#   stock status relative to ref points under each tac can be reported

# This model is called [in run-model-v1]* using purrr::map2df() which loops
#  over posterior samples and future TACs, and builds the data-frame
#  needed to make decision tables
#   * This model will be called by a yet-to-be-written function

# arguments:
# posteriors = a list of posterior outputs from one mcmc posterior sample
# tac = the future Total Allowable Catch for this projection

# Returns:
# A one-row data frame for the tac
#  all the variables of interest for the decision tables
#  If pyr>1, projections for subsequent years will be appended columnwise


# TODO:
# Add reference point calcs and calculate projected stock status

library(here)
library(RTMB)
source(here("R","get_ftdd.R"))

# Load documentation and inputs
devtools::document()
devtools::load_all()

# Get global inputs
dat <- pcod2020dat # Data inputs. Use ?pcod2020dat to see definitions
ctl <- pcod2020ctl # Control inputs. Use ?pcod2020ctl to see definitions. Not all used in d-d model
pfc <- pcod2020pfc # Control inputs for projections. Use ?pcod2020pfc to see definitions

########################################################################
# For testing only. DELETE AFTER TESTING!
# These are the function arguments
# Take first posterior sample
posteriors <- readRDS(here("outputs","MCMC_outputs_bysample.rda"))[[1]]
tac <- 100
########################################################################

project_model <- function(posteriors,
                          tac){

  getAll(dat,ctl,pfc) # RTMB function. Puts arguments into global space

  # Get the years for the projection under the tac
  # The model projects biomass into nyr+1 based on catches to nyr.
  #  pyr determines the number of projection years *after* this for
  #    decision tables
  npyr <- posteriors$proj_years # number of projection years (set by user)
  pyr  <- dat$nyr+1+posteriors$proj_years # actual final projection year
  pyrs <- (dat$nyr+1+1):(dat$nyr+1+posteriors$proj_years) # actual years of projection period
  nyrs <- length(dat$syr:dat$nyr) # number of historical years

  # 1. Get posterior estimates of leading parameters
  ro        <- exp(posteriors$log_ro)
  steepness <- posteriors$h
  m         <- exp(posteriors$log_m)
  tau       <- posteriors$tau

  # 2. Set up the time series vectors for model variables of interest
  # Fill the historical period and add projections years
  # For biomas and numbers, the historical period includes the one year
  #   projection arising from catch in the final historical year
  biomass <- c(posteriors$biomass, rep(NA,npyr)) # this is p_bt in iscam
  numbers <- c(posteriors$numbers, rep(NA,npyr)) # this is p_N in iscam
  recruits <- c(posteriors$recruits, rep(NA,npyr)) # this is p_rt in iscam
  ft <- c(posteriors$Ft, rep(NA,npyr)) # this is p_ft in iscam
  surv <- c(posteriors$surv, rep(NA,npyr)) # survival rate, this is p_S in iscam

  # 3. Simulate population into the future under constant tac
  for(i in (nyrs+1):(nyrs+1+npyr)){

    # recruits
    set.seed(99+i)
    xx  <- rnorm(1)*tau
    sbt <- biomass[i-dat$kage]

    # BH (1) or Ricker (2)
    if(ctl$misc[2]==1)recruits[i] <- (posteriors$alpha.sr*sbt/(1.+posteriors$beta.sr*sbt))*exp(xx-0.5*tau*tau)
    if(ctl$misc[2]==2)recruits[i] <- (posteriors$alpha.sr*sbt*exp(-posteriors$beta.sr*sbt))*exp(xx-0.5*tau*tau)

    #numbers and biomass
    # for dat$nyr+1 the projected values using final historical year catch
    #   are already included in the biomass and numbers vectors
    if(i>(nyrs+1)){
       biomass[i] <- (surv[i-1]*(dat$rho.g*biomass[i-1]+alpha.g*numbers[i-1])+dat$wk*recruits[i])
       numbers[i] <- surv[i-1]*numbers[i-1]+recruits[i]
     }

      # get_ftdd is defined in the get_ftdd.R
      # Newton-Raphson algorithm solves baranov eqn given a catch
      ft[i] = get_ftdd(tac,m,biomass[i])

      # Test get_ftdd with Baranov equation
      # put in the Newton solution for F and make sure you get the tac back
      # testf = ft[i]
      # testc = (testf/(testf+m))*(1-exp(-m-testf))*biomass[i]
      # tac
      # testc --> should be identical to tac


      #Calculate survival for next projection year
      surv[i] = exp(-(m+ft[i]))
    }


  # REPORT_SECTION
  output <- as.data.frame(matrix(nrow=1,ncol=3))
  output[1,1] <- tac
  output[1,2] <- pyr
  output[1,4] <- posteriors$biomass[1]


  colnames(output) <- c("TAC","pyr","B")

  return(output)

  } # end model
