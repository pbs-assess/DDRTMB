# This is a development script for a delay-difference model in rtmb
# The first draft attempts to reproduce the iscam delay difference model from
#   the 2018 Pacific Cod stock assessment for Haida Gwaii/Queen Charlotte
#   Sound stock in Pacific Canada (Area 5ABCD) [published in 2020]

# **For this first version, this is a ONE group, ONE area, ONE sex model
# Therefore, we can strip out a lot of the indexing in the iscam inputs**
# M is fixed, not time-varying - adapt later for tvm

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# IMPORTANT: BEFORE RUNNING THIS CODE, PLEASE SOURCE devs/filter-inputs.R

# This will load the raw 2020 5ABCD Pacific Cod inputs into the package
# Need to do this while package is in development
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Package name: DDRTMB

# Authors: Robyn Forrest (RF) and Catarina Wor (CW) (Pacific Biological Station, Nanaimo, Canada)

# Date created:  May 15, 2024
# Last Modified: June 28, 2024

# Notes:
# - The iscam input files are already loaded into the package:
#   1. pcod2020dat (the input data)
#   2. pcod2020ctl (controls for running the model,
#      including prior descriptions and misc settings)
#   3. pcod2020pfc (controls for projections and decision tables)
#
# - The input files include some inputs for the iscam age structured model
#       that are not used in the delay-difference implementation. The input
#       files will eventually be customised
# - Iscam uses an errors-in-variables approach to partition observation error,
#       with multiplicative weighting of each index observation using annual CVs
# - Want to change this to additive weightings (as per SS3), then also explore
#       state space options. But first try to reproduce the iscam results!

# TODO (move these to Issues on gitHub repo):
# 1. Translation of iscam to RTMB:
#  - CHECK all estimated parameters are included in pars
#  - CHECK length of log_rec_devs and log_init_rec_devs
#  - CHECK all test calcs. Catches are a little bit off. Could be rounding in iscam rep and par files
#  - Can probably relax requirement of setting nmeanwt to 1 when no mean weight data
#  - Probably don't need all the counters
#  - Check that pars list is complete
#  - Move some of model sections into separate functions
#  - Implement MCMC
#  - Check predicted rt against rep file
#  - Check source of fished equilibrium equation and check code (in calcNumbersBiomass_deldiff)

# 2. Potential model changes:
#  - Tidy up the three recruitment parameters - currently set to all be the same as per Paul Starr's request
#  - Transform normal space parameters to log?
#  - Need for Jacobian transformations?
#  - Look at bias correction (see Thorson and Kristensen paper - but not needed for bayesian)
#  - Look at alternate settings for Errors in Variables (e.g., weights additive instead of multiplicative)
#  - Look at state-space implementation
#  - Add time-varying M
#  - Think about how to make this a multi-species, multi-area model? (MICE)
#  - Some of the architecture is already in iscam
#  - consider a stan version
#  - consider a delay differential model in continuous time (see CJW correspondence)

# 3. Graphic outputs and diagnostics
#  - Coordinate with Sean, Nick, Catarina, others, ... for standardized set of visualizations of outputs and diagnostics

# Document and build package (these are also buttons in Rstudio)
#    this will incorporate new functions saved in the R and data folders

# Load documentation and inputs
devtools::document()
devtools::load_all()

library(here)
library(tidyverse)
library(RTMB)

# FOR TESTING
# source("devs/load-models.R")
# pcod2020rep<-read.report.file("data-raw/iscam.rep")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set up data and parameter controls

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Data inputs and controls
# DATA_SECTION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Input lists (*run filter-inputs.R first*)
dat <- pcod2020dat # Data inputs. Use ?pcod2020dat to see definitions
ctl <- pcod2020ctl # Control inputs. Use ?pcod2020ctl to see definitions
pfc <- pcod2020pfc # Control inputs for projections. Use ?pcod2020pfc to see definitions
nyrs <- dat$nyr-dat$syr+1
yrs <-  dat$syr:dat$nyr
ages <-  dat$sage:dat$nage

# get the number of and index for commercial (fishery) fleets
nfleet <- 0 # number of fishing fleets (not surveys)
for(i in 1:dat$ngear){
  nfleet <- ifelse(dat$alloc[i]>0, nfleet+1, nfleet)
}

# Index to identfy which gears are fishing fleets
fleetindex <- which(dat$alloc>0)

# get number of estimated log_ft parameters
ft_count = dat$nctobs # one estimated ft per catch obs

# Create a lookup for years because ADMB works with actual years as index
year_lookup <- as.data.frame(cbind(yrs, 1:nyrs))
colnames(year_lookup) <- c("year", "year_index")

# Get length and weight at age for first year (for initializing population at non-equilibrium, ctl$misc[5]==0)
la <- dat$linf*(1. - exp(-dat$k*(ages-dat$to)))
wa <- dat$lwscal*la^dat$lwpow
d3_wt_avg <- wa # just to be consistent with rep file

# Get settings for priors
# Leading parameters
num_params <- ctl$num.params
theta_control <- ctl$params
prior_settings_q <- ctl$surv.q

# Tiny number to stop logs breaking in some places
TINY <- 1.e-08

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Parameters
# PARAMETER_SECTION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par <- list()
# Begin by using initial values specified in pcod2020ctl
# Leading parameters
par$log_ro <- theta_control[1,1] # log unfished recruitment (syr)
par$h <- theta_control[2,1] # steepness
par$log_m <- theta_control[3,1] # log natural mortality
# In the current version of the P cod model, these are set the same as ro so do not estimate
#   par$log_avgrec <- theta_control[4,1] # log average recruitment (syr+1 to nyr)
#   par$log_recinit <- ptheta_control[5,1] #l og average of initial recruitments to fill first year if population if population is not unfished at syr
par$rho <- theta_control[6,1] # Errors in Variables: fraction of the total variance associated with observation error
par$kappa <- theta_control[7,1] # Errors in Variables: total precision (inverse of variance) of the total error.
# TODO: Check these are dimensioned and initialized correctly
par$log_ft_pars <- rep(0,nyrs) # estimated log fishing mortalities (total across fleets)
par$init_log_rec_devs <- rep(0,(dat$nage-dat$sage + 1)) # vector(length=dat$nage - dat$sage + 1) # I think this is length nage-sage+1 (i.e., length 2:9)
par$log_rec_devs <- rep(0,nyrs) # vector(length=nyrs)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. Model
# Eventually define some components as functions and move to separate R scripts
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
model <- function(par){

  getAll(par,dat, pfc) # RTMB function. Puts arguments into global space
  jnll <- 0 # initialize joint neg log likelihood
  pll  <- 0 # initialize priors component of neg log likelihood

  # Pseudocode from iscam
  # 1. initParameters()
  # 2. calcTotalMortality_deldiff()
  # 3. calcNumbersBiomass_deldiff()
  # 4. calcFisheryObservations_deldiff()
  # 5. calcSurveyObservations_deldiff()
  # 6. calcStockRecruitment_deldiff()
  # 7. calcAnnualMeanWeight_deldiff() //RF added this for P cod - only gets added to objective function if cntrl(15)==1
  # 8. calcObjectiveFunction()

#|---------------------------------------------------------------------|
  # 1. initParameters
  theta <- c(exp(par$log_ro),
             par$h,
             exp(par$log_m),
             par$rho,
             par$kappa)

  ro        <- theta[1]
  steepness <- theta[2]
  m         <- theta[3]
  rho       <- theta[4]
  varphi    <- sqrt(1.0/theta[5])
  sig       <- sqrt(rho)*varphi
  tau       <- sqrt(1.0-rho)*varphi

  # Fixed parameters
  # A decision was made in 2018 to fix these two parameters to be the same
  #   as log_ro (P. Starr). May revisit this later
  log_avgrec <- par$log_ro
  log_recinit <- par$log_ro # RF: I don't think it makes sense to use this one
  sig_c <- ctl$misc[4] # sd in catch likelihood
  sig_w <- ctl$weight.sig # sd in mean weight likelihood (called weight_sig in iscam)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~ TESTING VALUES FROM data-raw/iscam.rep and data-raw/iscam.par (2020 PCOD RESULTS) ~
  # DELETE THIS ONCE MODEL EQUATIONS ARE TESTED
  ro        <- 3376.55
  steepness <- 0.809901
  m         <- 0.306306
  rho       <- 0.058824
  varphi    <- 0.824621
  log_avgrec <- log(ro)
  log_recinit <- log(ro)
  # from iscam.par file
  log_ft_pars <- c(-2.29361, -1.81516, -1.52265, -1.75990, -1.99738, -2.38600, -1.96696, -2.56048, -2.04428, -1.65530, -1.44354, -1.56225, -1.68573, -1.89923, -2.53602, -2.00452, -1.99661, -2.16381, -1.87796, -1.72526, -1.78542, -1.97597, -2.04832, -1.38198, -1.55781, -1.92887, -1.98627, -2.11772, -2.25886, -2.71592, -1.82257, -1.02055, -1.41985, -1.95312, -1.88346, -1.08038, -1.18403, -1.29623, -2.10493, -2.29368, -2.19147, -2.11395, -2.28967, -2.53718, -2.71805, -3.39986, -3.27054, -2.89189, -2.72193, -2.52443, -2.59201, -3.34633, -3.53683, -2.94841, -2.34753, -2.56502, -2.88457, -2.87729, -2.74026, -2.64641, -3.14998, -3.54423, -3.96561, -3.51804, -3.65543)
  init_log_rec_devs <- c(-0.297834, -0.195054, -0.126748, -0.0955628, -0.0934589, -0.108359, -0.130301, 1.04734)
  log_rec_devs <- c(1.05722, 1.10583, -0.139089, -0.165389, -0.298059, -0.336892, -0.173463, 2.84111, 0.284625, 0.163418, -0.0760200, -0.352092, -0.626335, -0.538303, -0.320139, -0.0816409, 2.69634, 0.0765257, 0.524992, 0.510128, 0.356662, 0.953328, 0.574398, 0.840802, 0.173325, 0.402038, 0.278233, -0.103700, 0.166054, 0.213154, 1.49743, 2.13800, -0.221516, -0.0713425, 0.874159, 1.27436, -0.245994, -0.775609, -0.898877, -0.701367, -0.142345, -0.829222, -0.954500, -1.11217, -1.11537, 0.209017, 0.409310, -0.409217, -0.845547, -1.24699, -1.39305, -1.25216, -0.294358, 0.668812, 0.131646, -0.489765, -0.691204, -0.667682, -0.629868, -0.792061, -0.796493, -0.646523, 0.347852, -0.110935, -0.232896)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
  Ft <- rep(0,nyrs) # a vector of total F for each year (in the ASM this would account for selectivity but for delay-diff, selectivity is knife-edged)

  # Work down catch obs matrix to associate Fs with catch for each gear
  # The catch matrix in the data file has 4 columns:
  #   year, gear, type (1=catch in weight (tonnes); 2=catch in numbers), catch
  # nctobs is from the dat list and indicates the total number of catch obs
   for(ii in 1:nctobs){
      # Set up counters
      iyear <- catch[ii,1] # actual year
      i <- as.integer(year_lookup[which(year_lookup[,1]==iyear),2]) # year index
      k <- catch[ii,2] # gear

      ftmp    <- exp(log_ft_pars[ii]) # log_ft_pars has length nctobs
      # ft is a matrix with ngear rows and nyr columns
      ft[k,i] <- ftmp # fishing mortality for gear k in year i (in ASM this is modified by selectivity)
      Ft[i]   <- Ft[i] + ftmp # Total fishing mortality in year i
   } # end for ii

  Mt <- rep(m,nyrs) # natural mortality
  Zt <- Ft+Mt       # total mortality
  surv <- exp(-Zt) # fished survival rate

  # End calcTotalMortality_deldiff
#|---------------------------------------------------------------------|

#|---------------------------------------------------------------------|
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
  numbers <- vector(length=nyrs+1)
  biomass <- vector(length=nyrs+1)
  sbt <- vector(length=nyrs+1) # can eliminate this eventually. It's the same as biomass in the dd model
  # FIXME: I think this should just be from syr:(nyr-sage)
  log_rt <- vector(length=nyrs) # estimated log recruits
  annual_mean_wt <- vector(length=nyrs)
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

  # Three options for how first year is initialized (set in ctrl$misc[5]):
  # 0. Unfished and not at equilibrium
  # 1. Unfished equilibrium
  # 2. Fished and at equilibrium with Ft in first year
  if(ctl$misc[5]==0){
    # Unfished and not at equilibrium - initialize with age structure as in ASM
    # Set up a vector of n-at-age for the first year (length sage:nage) with devs
    # Decay exp(log_recinit + init_log_rec_devs[j-1]) by exp(-M) to fill n-at-age in first year
    tmp_nAge <- vector(length=length(sage:nage)) # Numbers at age in first year

    # First age is exp(log_avgrec + log_rec_devs[1]) because this is recruits in year 1
    tmp_nAge[1] <- exp(log_avgrec + log_rec_devs[1])

    # Now fill in the other numbers at age for year 1
    for(j in 2:length(tmp_nAge)){
      tmp_nAge[j] <- exp(log_recinit + init_log_rec_devs[j-1]) * exp(-Mt[1]*(j-1))
    }
    tmp_nAge[j] <- tmp_nAge[j]/(1 - exp(-Mt[1])) # plus group

    ## TEST (see iscam calcNumbersAtAge function) - make sure tr is the same as log(tmp_nAge)
    ## DELETE THIS IN LATER VERSIONS
    ## A little tricky  because ADMB indexes actual ages or years rather than indexes
    ## **YES, calculations are identical**
    # tr <- vector(length=length(sage:nage)) # iscam ASM's version of log(tmp_n_Age)
    # lx <- vector(length=length(sage:nage)) # unfished survivorship at age
    # lx[1] <- 1
    # for(j in 2:length(sage:nage)){
    #   lx[j] = lx[j-1] * exp(-Mt[1])
    # }
    # lx[j] <- lx[j]/(1-exp(-Mt[1]))
    #
    # tr[1]  <- log_avgrec+log_rec_devs[1]
    # tr[2:length(tr)] <- log_recinit +init_log_rec_devs
    # tr[2:length(tr)] <- tr[2:length(tr)]+log(lx[2:length(tr)])
    # exp(tr)

    # Add up the numbers at age to get numbers in year 1
    # RF confirmed numbers, biomass and mean wt calculations with rep file
    numbers[1] <- sum(tmp_nAge)
    biomass[1] <- sum(tmp_nAge*d3_wt_avg) # d3_wt_avg calculated in data section from vonB parameters
    annual_mean_wt[1] <- biomass[1]/numbers[1]

    #  Initialise log recruits
    log_rt[1] <- log_avgrec+log_rec_devs[1] # this is just the same as log(tmp_nAge[1])
  }
  if(ctl$misc[5]==1){
    # Unfished equilibrium - **not tested for P. cod**
    numbers[1] <- no
    biomass[1] <- bo
    annual_mean_wt[1] <- biomass[1]/numbers[1]
    log_rt[1] <- log(ro) # Shouldn't this be plus log_rec_devs[1]???
  }
  if(ctl$misc[5]==2){
    # Fished at equilibrium with Ft in first year **not tested for P. cod**
    # NEED TO CHECK THESE CALCS. SOURCE IS 2004 P.COD ASSESSMENT OR H&W 1992 - need to check
    # RF: I think I have some code from SJDM somewhere with these calcs
    sfished = surv[1] # equilibrium survivorship at initial fishing mortality (gear 1 commercial fishery)
    annual_mean_wt[1] = (sfished*alpha.g + wk*(1.-sfished))/(1-rho.g*sfished)

    biomass[1] = -(annual_mean_wt[1]*(wk*alpha.sr-1)
                   +sfished*(alpha.g + rho.g * annual_mean_wt[1]))/
                  (beta.sr*(sfished*alpha.g + sfished*rho.g*annual_mean_wt[1]-
                  annual_mean_wt[1]))
    numbers[1] = biomass[1]/annual_mean_wt[1]

   # log rt originally missing from this option
   # chose log_avgrec as placeholder-- dangerous if fishing in first year and before was very high.
   log_rt[1] <- log_avgrec[1] #+ log_rec_devs[1]  # Shouldn't this be plus log_rec_devs[1]???
  }
  if(!ctl$misc[5] %in% 0:2){
    stop("Starting conditions not specified! ctl$misc[5] must be set to 0, 1 or 2 (see ?pcod2020ctl for definitions).")
  }

  # Set sbt[1] <- biomass[1]
  # Eventually can delete this duplication but for now will be easier
  # bc sbt appears elsewhere in code as hangover from ASM
  sbt[1] <- biomass[1]

  # Time dynamics
  for(i in 2:nyrs){

    # Update recruits
    log_rt[i] <- log_avgrec +log_rec_devs[i]

    # Update biomass and numbers
    biomass[i] <- surv[i-1]*(rho.g*biomass[i-1]+alpha.g*numbers[i-1])+
								wk*exp(log_rt[i]) # eq. 9.2.5 in HW
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
  ct <- vector(length=nctobs) # predicted catch
  eta <- vector(length=nctobs) # catch residuals

  for(ii in 1:nctobs){
    # Set up counters
    iyear <- catch[ii,1] # actual year
    i <- as.integer(year_lookup[which(year_lookup[,1]==iyear),2]) # year index
    k <- catch[ii,2] # gear
    m <- catch[ii,3] # type: 1=catch in weight; 2=catch in numbers
    d_ct <- catch[ii,4] # observed catch

    # Baranov catch equation
    if(m==1) {
      # catch in weight
      ct[ii] = (ft[k,i]/(Mt[i] + ft[k,i]))*(1-exp(-Mt[i]-ft[k,i]))*biomass[i]

    }
    if(m==2){
      # catch in numbers
      ct[ii] = (ft[k,i]/(Mt[i] + ft[k,i]))*(1-exp(-Mt[i]-ft[k,i]))*numbers[i]
    }
    if(!m %in% 1:2){
      stop("Catch type must be 1 (weight) or 2 (numbers). Set in column 3 of dat$catch")
    }

    # catch residual
    eta[ii] = log(d_ct+TINY) - log(ct[ii]+TINY)

    # NOTE: Catches are not exactly as in rep file - could be rounding in the reported
    # log_ft_pars from the par file. GO BACK AND CHECK ALL CALCS AND VALUES
  } # end for ii

  # End calcFisheryObservations_deldiff
#|---------------------------------------------------------------------|


#|---------------------------------------------------------------------|
  # 5. calcSurveyObservations_deldiff()
  # Purpose: This function calculates predicted survey observations for each year

  # Needed to determin if q is random walk
  q_prior <- prior_settings_q[1,]

  # Set up vector for mle qs (per Walters&Ludwig 1993 https://doi.org/10.1139/f94-07)
  q <- vector(length=nit) # vector of q for each survey

  # set up lists for storing residuals and predicted indices
  # these are ragged arrays in iscam
  epsilon <- list() # residuals
  it_hat <- list() # predicted indices
  qt <- list() # time varying q (if ctl$surv.q prior type == 2)

  #loop over surveys to create list
  for(kk in 1:nit){
    epsilon[[kk]] <- vector(length=nitnobs[kk])
    it_hat[[kk]]  <- vector(length=nitnobs[kk])
    qt[[kk]]  <- vector(length=nitnobs[kk])
    epsilon[[kk]][1:nitnobs[kk]] <- 0.
    it_hat[[kk]][1:nitnobs[kk]] <- 0.
    qt[[kk]][1:nitnobs[kk]] <- 0.
  }

  # Now loop through gears and index obs to predict survey obs
  for(kk in 1:nit){
    # Vulnerable numbers
    V <- vector(length=nitnobs[kk])
    V[1:nitnobs[kk]] <- 0.

    # not sure why we need these
    #nz <- 0 # counter for number of observations
    #iz <- 1  # index for first year of data for prospective analysis

    for(ii in 1:nitnobs[kk]){
      iyear    <- indices[[kk]][ii,1] # actual year (not used)
      i <- as.integer(year_lookup[which(year_lookup[,1]==iyear),2]) # year index - need this to match biomass to survey obs
      k    <- indices[[kk]][ii,3] # gear
      di   <- indices[[kk]][ii,5] # timing

      #nz <- nz+1  # counter for number of observations

      # TODO: check these correctly matched up with years - yes, checked against rep file
      z = ft[k,i]+Mt[i]
      # Adjust numbers and biomass for survey timing. If di=0, no adjustment
      Np = numbers[i] * exp( -z * di)
      Bp = biomass[i] * exp( -z * di)

      # Two different survey types: 1=prop to numbers; 2=prop to biomass
      if(survtype[kk]==1){
        V[ii] <- Np
      }
      if(survtype[kk]==2){
        V[ii] <- Bp
      }
      if(!survtype[kk] %in% 1:2){
        stop(("Survey type must be 1 (survey proporional to numbers) or 2 (survey proporional to biomass). Set in dat$survtype"))
      }
    } #end of ii loop

    it 	<- t(indices[[kk]])[2,1:nitnobs[kk]] # index
    wt 	<- t(indices[[kk]])[4,1:nitnobs[kk]] # index weight
    wt 	<- wt/sum(wt) # normalized weight (sum to 1) # Q: here the weights are normalized by the sum but in the likelihood, iscam normalises it_wt by the mean (in the data section)

    # get mle q (same as for ASM, from Walters&Ludwig 1993 https://doi.org/10.1139/f94-07)
    zt 	= log(it) - log(V[1:nitnobs[kk]])
    zbar = sum(zt*wt)
    q[kk] = exp(zbar)

    # survey residuals - checked against rep file
    epsilon[[kk]][1:nitnobs[kk]] = zt - zbar
    it_hat[[kk]][1:nitnobs[kk]] = q[kk] * V[1:nitnobs[kk]]

    # SPECIAL CASE: penalized random walk in q.
    # !!!NOT TESTED!!! This is dimensioned correctly but have not checked calcs
    if(q_prior[kk]==2 ){
        epsilon[[kk]][1:nitnobs[kk]] <- 0 # initialize epsilon

        # iscam ADMB code:
         # fd_zt <- first_difference(zt)
        # From the admb source code, looks like first_difference returns a vector of
        # i+1 - i:
        # i.e., differences.elem(i) = values.elem(i + 1) - values.elem(i);
        fd_zt <- diff(zt)
        zw_bar <- sum(fd_zt*wt[1:(nitnobs[kk]-1)])
        epsilon[[kk]][1:(nitnobs[kk]-1)] <- fd_zt - zw_bar
        qt[[kk]][1] = exp(zt[1])

        for(ii in 2:nitnobs[k]){
          qt[[kk]][ii] = qt[[kk]][ii-1] * exp(fd_zt[ii-1])
        }
        it_hat[[kk]][1:nitnobs[kk]] = qt[[kk]][1:nitnobs[kk]]*V[1:nitnobs[kk]]
      }
  } #end kk loop

  # End calcSurveyObservations_deldiff
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
  #rt     <- vector(length=nyrs-sage)  # estimated recruits from calcNumbersBiomass_deldiff()
  #delta  <- vector(length=nyrs-sage)  # residuals between estimated R and R from S-R curve (process err)
  tmp_rt <- vector(length=nyrs) # recruits derived from stock-recruit model
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

#|---------------------------------------------------------------------|
  # 7. calcAnnualMeanWeight_deldiff() //RF added this for P cod - only gets added to objective function if cntrl(15)==1
  # Purpose: This function calculates the mean weight of the catch for each year, gear by dividing the total
  #          biomass by the total numbers

  # Set up lists for obs and predicted annual mean weights
  # these are ragged arrays in iscam
  annual_mean_weight     <- list() # observed
  obs_annual_mean_weight <- list() # predicted
  epsilon_mean_weight <- list() # residuals

  #loop over series to create list
  for(kk in 1:nmeanwt){
    annual_mean_weight[[kk]]     <- vector(length=nmeanwtobs[kk])
    obs_annual_mean_weight[[kk]] <- vector(length=nmeanwtobs[kk])
    epsilon_mean_weight[[kk]]    <- vector(length=nmeanwtobs[kk])
    annual_mean_weight[[kk]][1:nmeanwtobs[kk]] <- 0.
    obs_annual_mean_weight[[kk]][1:nmeanwtobs[kk]] <- 0.
    epsilon_mean_weight[[kk]][1:nmeanwtobs[kk]] <- 0.
  }

  # loop through series with empirical annual mean weight data
  for(kk in 1:nmeanwt){

    Vn <- vector(length = nmeanwtobs[kk])	      # Vulnerable numbers to gear
    Vb <- vector(length = nmeanwtobs[kk])	      # Vulnerable biomass to gear
    Vn[1:nmeanwtobs[kk]] <- Vb[1:nmeanwtobs[kk]] <- 0

    # Loop through observations
    for(ii in 1:nmeanwtobs[kk]){

      iyear    <- meanwtdata[[kk]][ii,1]  # actual year - not used
      i <- as.integer(year_lookup[which(year_lookup[,1]==iyear),2]) # year index - need this to match meanwt to obs
      k    <- meanwtdata[[kk]][ii,3]  # gear
      di   <- meanwtdata[[kk]][ii,4]  # timing

      ws  = exp(-Zt[i]*di)   # Total mortality that accounts for timing
      wN  = numbers[i]*ws
      wB  = biomass[i]*ws
      Vn[ii] <- wN
      Vb[ii] <- wB

      # RF Checked against rep file
      annual_mean_weight[[kk]][ii] <- Vb[ii]/Vn[ii]
      obs_annual_mean_weight[[kk]][ii] <-  meanwtdata[[kk]][ii,2]	  # fill a list of vectors with observed annual mean weights
      #  residual
      epsilon_mean_weight[[kk]][ii] <-  log(annual_mean_weight[[kk]][ii]) - log(obs_annual_mean_weight[[kk]][ii])
    }	# end ii loop
  } # end kk loop

  # End calcAnnualMeanWeight_deldiff
#|---------------------------------------------------------------------|

#|---------------------------------------------------------------------|
# calcObjectiveFunction();

  # Likelihood for catch
  for(ii in 1:nctobs){
    jnll <- jnll - dnorm(eta[ii], 0.0, sig_c)
  }

  # Likelihood for relative abundance indices
  # loop over surveys
  for(kk in 1:nit){
    sig_it <- rep(0,nitnobs[kk]) # vector for weights for each obs in survey k

    # Loop over observations in the survey
    for(ii in 1:nitnobs[kk])
    {
      # Get the weightings. Note that iscam normalizes it_wt in the datasection
      # by dividing by the mean. But above, where q is calculated, the weights (called wt)
      # are normalized by dividing by sum
      # for now, follow iscam
      it_wt <- dat$indices[[k]][ii,4]
      tmp_mean <- mean(it_wt)
      it_wt <- it_wt/tmp_mean # Normalise. This happens on L466 of devs/iscam.tpl
      sig_it[ii] <- sig/it_wt # divide global sig by individual weights (which are inverted, so divide)

      # update likelihood
      jnll <- jnll - dnorm(epsilon[[kk]][ii], 0.0, sig_it[ii])
    }
  }

  # Likelihood for recruitment
  for(ii in 1:length(rt)){
    jnll <- jnll - dnorm(delta[ii], 0.0, tau)
  }

  # Likelihood for mean weight
  # We are entering the likelihood in log space here - do we need a Jacobian transformation?
  for(kk in 1:nmeanwt){
    for(ii in 1:nmeanwtobs[kk]){
        jnll <- jnll - dnorm(epsilon_mean_weight[[kk]][ii], 0.0, sig_w)
    }
  }

 #==============================================================================================
 # ~PRIORS~
 #==============================================================================================

 for(i in 1:ctl$num.params){
  ptype <- theta_control[i,5]





 }





 jnll <- jnll - pll
 # End calcObjectiveFunction
#|---------------------------------------------------------------------|


  jnll # return joint neg log likelihood
}

model(pars)

## MakeADFun builds the graph, basically "compiles" the model with random effects identified
## from TMB help: map = List defining how to optionally collect and fix parameters
## Means you can fix some instances of a vector/matrix of parameters and estimate ones with the same factor id to be the same
## Might be good for q or selectivity for example when you want all the same value for a given age

# from babysam
#obj <- MakeADFun(f, par, random=c("logN", "logF", "missing"), map=list(logsdF=as.factor(rep(0,length(par$logsdF)))), silent=FALSE)
obj <- MakeADFun(f, pars, silent=FALSE)
# The optimization step - gets passed the parameters, likelihood function and gradients Makeby ADFun
opt <- nlminb(obj$par, obj$fn, obj$gr, control=list(eval.max=1000, iter.max=1000))
opt$objective

