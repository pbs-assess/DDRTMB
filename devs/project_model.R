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

# Returns
# A one-row data frame for the tac
#  all the variables of interest for the decision tables
#  If pyr>1, projections for subsequent years will be appended columnwise


# TODO:
# Add reference point calcs and calculate projected stock status

library(here)

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
tac <- 0
########################################################################

project_model <- function(posteriors,
                          tac){

  getAll(dat,ctl,pfc) # RTMB function. Puts arguments into global space

  # Get the years for the projection under the tac
  # The model projects biomass into nyr+1 based on catches to nyr.
  #  pyr determines the number of projection years after this for
  #    decision tables
  pyrs <- (dat$nyr+1):dat$nyr+1+posteriors$proj_years
  hyrs <- dat$syr:dat$nyr # historical years

  # 1. initParameters
  ro        <- exp(theta[1])
  steepness <- theta[2]
  m         <- exp(theta[3])
  tau       <- sqrt(1.0-rho)*varphi # 0.8 for P. cod

  # A decision was made in 2018 to fix these two parameters to be the same
  #   as log_ro (P. Starr). May revisit this later
  log_avgrec <- log(ro)
  log_recinit <- log(ro)

  # Get the Goodyear Compensation Ratio

  # This is the ratio of juvenile survival at low stock size (alpha) to unfished juvenile survival
  # Calculation differs between BH and Ricker S-R functions (see Forrest et al. 2010 and refs therein https://doi.org/10.1139/F10-077)
  # **IMPORTANT: Steve called this kappa in iscam**
  # **To avoid confusion with the leading variance parameter kappa**
  # **RF has renamed this parameter CR**
  if(ctl$misc[2]==1){
    # Beverton-Holt formulation
    CR <- (4*steepness)/(1-steepness)
  }
  if(ctl$misc[2]==2){
    # Ricker formulation
    CR <- (5*steepness)^1.25
  }
  if(!ctl$misc[2] %in% 1:2){
    stop("S-R relationship not specified! ctl$misc[2] must be set to 1 (BH) or 2 (Ricker).")
  }

  # End initParameters
  #|---------------------------------------------------------------------|

  #|---------------------------------------------------------------------|
  # 2. calcTotalMortality_deldiff();
  # Purpose: This function calculates fishing mortality, total mortality and annual
  #          survival rates S=exp(-Z) for each year.
  #          Z also is updated with time-varying natural mortality rates if
  #          specified by user.

  # Set up F objects
  # In the one fleet case for the delay-diff model, these are the same.
  ft <- matrix(0, nrow=ngear, ncol=nyrs) # a matrix of fishing mortality for each gear and year
  Ft <- numeric(nyrs) # a vector of total F for each year (in the ASM this would account for selectivity but for delay-diff, selectivity is knife-edged)

  # Work down catch obs matrix to associate Fs with catch for each gear
  # The catch matrix in the data file has 4 columns:
  #   year, gear, type (1=catch in weight (tonnes); 2=catch in numbers), catch
  # nctobs is from the dat list and indicates the total number of catch obs
  for(ii in 1:nctobs){
    # Set up counters
    iyear <- catch[ii,1] # actual year
    i <- as.integer(year_lookup[which(year_lookup[,1]==iyear),2]) # year index
    k <- catch[ii,2] # gear

    ftmp  <- exp(log_ft_pars[ii]) # log_ft_pars has length nctobs
    # ft is a matrix with ngear rows and nyr columns
    ft[k,i] <- ftmp # fishing mortality for gear k in year i (in ASM this is modified by selectivity)
    Ft[i]   <- Ft[i] + ftmp # Total fishing mortality in year i
  } # end for ii

  Mt <- rep(m,nyrs) # natural mortality
  Zt <- Ft+Mt       # total mortality
  surv <- exp(-Zt) # fished survival rate

  # End calcTotalMortality_deldiff
  # #|---------------------------------------------------------------------|

  # #|---------------------------------------------------------------------|
  # 3. calcNumbersBiomass_deldiff();
  # Purpose: This function calculates  total biomass and total numbers
  # according to the delay difference model equations from Hilborn and Walters.
  # Quantities are calculated for each year.

  # DD initialization
  # Equilibrium mean weight - obtained from weq = Beq/Neq and solving for weq
  # i.e., weq = [surv(alpha.Neq + rho.Beq + wk.R] / [surv.Neq + R]
  #  with substitutions Neq = Beq/weq and R = Neq(1 - surv)
  # From SJDM, also used by Sinclair in 2005 p cod assessment

  # Initialize numbers and biomass
  numbers <- numeric(length=nyrs+1)
  biomass <- numeric(length=nyrs+1)
  sbt <- numeric(length=nyrs+1) # can eliminate this eventually. It's the same as biomass in the dd model
  # FIXME: I think this should just be from syr:(nyr-sage)
  log_rt <- numeric(length=nyrs) # estimated log recruits
  annual_mean_wt <- numeric(length=nyrs)
  # Add recruitment for projection year ... assume it is average
  rnplus=exp(log_avgrec)

  # Initial Conditions
  snat <- exp(-m) # natural survival rate
  # Eqm mean weight
  wbar <- (snat*alpha.g + wk*(1-snat))/(1-snat*rho.g) # RF: calculation checked against rep file

  # Unfished numbers and biomass
  no <- ro/(1-snat)
  bo <- no*wbar # RF: calculation checked against rep file
  # added sbo for delay diff model - delete later
  sbo <- bo

  # Parameters of Stock-Recruit relationship
  # Maximum juvenile survival rate (same for BH and Ricker)
  alpha.sr <- CR*(ro/bo) # NOTE: *Steve called this so*. RF thinks this is confusing notation because the o implies unfished, whereas this is max juv survival rate when the biomass is close to zero

  # Beta parameter (density-dependence in juvenile survival)
  if(ctl$misc[2]==1){
    # Beverton-Holt
    beta.sr <- (CR -1)/bo
  }
  if(ctl$misc[2]==2){
    # Ricker
    beta.sr <- log(CR)/bo
  }

  # Time dynamics
  for(i in 2:nyrs){

    # Update recruits
    log_rt[i] <- log_avgrec +log_rec_devs[i]

    # Update biomass and numbers
    biomass[i] <- surv[i-1]*(rho.g*biomass[i-1]+alpha.g*numbers[i-1]) +wk*exp(log_rt[i]) # eq. 9.2.5 in HW
    numbers[i] <- surv[i-1]*numbers[i-1]+exp(log_rt[i])
    annual_mean_wt[i] <- biomass[i]/numbers[i]		# calculate predicted weight in dynamics - possible option to fit to it
    sbt[i] <- biomass[i]
  }

  # RF confirmed numbers, biomass and mean wt calculations with rep file
  biomass[nyrs+1]  <- (surv[nyrs]*(rho.g*biomass[nyrs]+alpha.g*numbers[nyrs]) + wk*rnplus)
  numbers[nyrs+1]  <- surv[nyrs]*numbers[nyrs]+rnplus
  sbt[nyrs+1]      <- biomass[nyrs+1] # set spawning biomass to biomass

  # End calcNumbersBiomass_deldiff
  #|---------------------------------------------------------------------|

  #|---------------------------------------------------------------------|
  # 4. calcFisheryObservations_deldiff()
  # Purpose: This function calculates commercial catches for each year and gear
  # The catch matrix in the data file has 4 columns:
  #   year, gear, type (1=catch in weight (tonnes); 2=catch in numbers), catch
  # nctobs is from the dat list and indicates the total number of catch obs

  # Set up predicted catch by gear and year
  # like log_ft_pars, this is one long vector of all catches from all gears
  ct <- numeric(length=nctobs) # predicted catch
  eta <- numeric(length=nctobs) # catch residuals

  for(ii in 1:nctobs){
    # Set up counters
    iyear <- catch[ii,1] # actual year
    i <- as.integer(year_lookup[which(year_lookup[,1]==iyear),2]) # year index
    k <- catch[ii,2] # gear
    mm <- catch[ii,3] # type: 1=catch in weight; 2=catch in numbers
    d_ct <- catch[ii,4] # observed catch

    # Baranov catch equation
    if(mm==1) {
      # catch in weight
      ct[ii] = (ft[k,i]/(Mt[i] + ft[k,i]))*(1-exp(-Mt[i]-ft[k,i]))*biomass[i]

    }
    if(mm==2){
      # catch in numbers
      ct[ii] = (ft[k,i]/(Mt[i] + ft[k,i]))*(1-exp(-Mt[i]-ft[k,i]))*numbers[i]
    }
    if(!mm %in% 1:2){
      stop("Catch type must be 1 (weight) or 2 (numbers). Set in column 3 of dat$catch")
    }

    } # end for ii

  # End calcFisheryObservations_deldiff
  #|---------------------------------------------------------------------|

  #|---------------------------------------------------------------------|
  # 6. calcStockRecruitment_deldiff()
  # Purpose:
  # This function is used to derive the underlying stock-recruitment
  # relationship that is ultimately used in determining MSY-based reference
  # points.  The objective of this function is to determine the appropriate
  # Ro, Bo and steepness values of either the Beverton-Holt or Ricker  Stock-
  #   Recruitment Model:
  #  Beverton-Holt Model
  #  Rt=k*Ro*St/(Bo+(k-1)*St)*exp(delta-0.5*tau*tau) \f$
  #
  #   Ricker Model
  #  Rt=so*St*exp(-beta*St)*exp(delta-0.5*tau*tau) \f$

  # Set up vectors
  # Don't need to declare first two
  #rt     <- numeric(length=nyrs-sage)  # estimated recruits from calcNumbersBiomass_deldiff()
  #delta  <- numeric(length=nyrs-sage)  # residuals between estimated R and R from S-R curve (process err)
  tmp_rt <- numeric(length=nyrs) # recruits derived from stock-recruit model
  tmp_rt[1:nyrs] <- 0

  # get the process error term from the errors in variables parameters
  # [For the P. cod assessment, rho and varphi are set to give tau=0.8 and obs error=0.2]
  tau <- sqrt(1-rho)*varphi

  # counter
  iicount <- 0

  for(i in 1:nyrs){
    iicount <- iicount+1

    if(ctl$misc[2]==1){
      # Beverton-Holt recruitment
      if(iicount <= kage){
        tmp_rt[i] =  alpha.sr*sbt[1]/(1.+beta.sr*sbt[1])
      }else{
        tmp_rt[i] = alpha.sr*sbt[i-kage]/(1.+beta.sr*sbt[i-kage])
      }
    }
    if(ctl$misc[2]==2){
      # Ricker recruitment
      if(iicount <= kage){
        tmp_rt[i] =  alpha.sr*sbt[1]*exp(-beta.sr*sbt[1])
      }else{
        tmp_rt[i] = alpha.sr*sbt[i-kage]*exp(-beta.sr*sbt[i-kage])
      }
    }
    if(!ctl$misc[2] %in% 1:2){
      stop("Recruitment model must be 1 (Beverton-Holt) or 2 (Ricker). Set in ctl$misc[2]")
    }

  } # end year loop i

  # estimated recruits from calcNumbersBiomass_deldiff()
  # RF Checked against rep file
  rt <- exp(log_rt[(sage+1):nyrs])

  # Calculate delta: process errors (deviations from S-R function)
  # RF: ADD A NOTE WHY BIAS CORRECTION ADDED HERE
  delta <- log(rt)-log(tmp_rt[(sage+1):nyrs])+0.5*tau*tau # RF Checked against rep file

  # End calcStockRecruitment_deldiff
  #|---------------------------------------------------------------------|

  # REPORT_SECTION
  output <- as.data.frame(matrix(nrow=1,ncol=3))
  output[1,1] <- tac
  output[1,2] <- pyr
  output[1,4] <- posteriors$biomass[1]


  colnames(output) <- c("TAC","pyr","B")

  return(output)

  } # end model
