# This is a development script for a delay-difference model in rtmb
# The first draft attempts to reproduce the iscam delay difference model from
#   the 2018 Pacific Cod stock assessment for Haida Gwaii/Queen Charlotte
#   Sound stock in Pacific Canada (Area 5ABCD) [published in 2020]

# **For this first version, this is a ONE group, ONE area, ONE sex model
# Therefore, we can strip out a lot of the indexing in the iscam inputs**
# M is fixed, not time-varying - adapt later for tvm

# Package name: DDRTMB

# Authors: Robyn Forrest (RF) and Catarina Wor (CW) (Pacific Biological Station)

# Date created:  May 15, 2024
# Last Modified: May 27, 2024

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
#  - Can probably relax requirement of setting nmeanwt to 1 when no mean weight data
#  - Probably don't need all the counters
#  - Check that pars list is complete
#  - Move some of model sections into separate functions
#  - Implement MCMC

# 2. Potential model changes:
#  - Tidy up the three recruitment parameters - currently set to all be the same as per Paul Starr's request
#  - Transform normal space parameters to log
#  - Need for Jacobian transformations?
#  - Look at bias correction (see Thorson and Kristensen paper)
#  - Look at alternate settings for Errors in Variables (e.g., weights additive instead of multiplicative)
#  - Look at state-space implementation
#  - Add time-varying M
#  - Think about how to make this a multi-species, multi-area model? (MICE)
#     - Some of the architecture is already in iscam


# 3. Graphic outputs and diagnostics
#  - Work with Sean, Nick, Catarina, others, ... for standardized set of visualizations of outputs and diagnostics

# Document and build package (these are also buttons in Rstudio)
#    this will incorporate new functions saved in the R and data folders

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# IMPORTANT: BEFORE RUNNING THIS CODE, PLEASE SOURCE devs/filter-inputs.R

# This will load the raw 2020 5ABCD Pacific Cod inputs into the package
# Need to do this while package is in development
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load documentation and inputs
devtools::document()
devtools::load_all()

library(here)
library(tidyverse)
library(RTMB)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set up data and parameter controls
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Data inputs and controls
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Input lists (run filter-inputs.R first)
dat <- pcod2020dat # Data inputs. Use ?pcod2020dat to see definitions
ctl <- pcod2020ctl # Control inputs. Use ?pcod2020ctl to see definitions
pfc <- pcod2020pfc # Control inputs for projections. Use ?pcod2020pfc to see definitions
nyrs <- dat$nyr-dat$syr+1
yrs <-  dat$syr:dat$nyr

# get the number of and index for commercial (fishery) fleets
nfleet <- 0 # number of fishing fleets (not surveys)
for(i in 1:dat$ngear){
  nfleet <- ifelse(dat$alloc[i]>0, nfleet+1, nfleet)
}

# Index to identfy which gears are fishing fleets
fleetindex <- which(dat$alloc>0)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Parameters
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par <- list()
# Begin by using initial values specified in pcod2020ctl
# Check all estimated pars are here
# Leading parameters
par$log_ro <- pcod2020ctl$params[1,1] # log unfished recruitment (syr)
par$h <- pcod2020ctl$params[2,1] # steepness
par$log_m <- pcod2020ctl$params[3,1] # log natural mortality
# In the current version of the P cod model, these are set the same as ro so do not estimate
#   par$log_avgrec <- pcod2020ctl$params[4,1] # log average recruitment (syr+1 to nyr)
#   par$log_recinit <- pcod2020ctl$params[5,1] #l og average of initial recruitments to fill first year if population if population is not unfished at syr
par$rho <- pcod2020ctl$params[6,1] # Errors in Variables: fraction of the total variance associated with observation error
par$kappa <- pcod2020ctl$params[7,1] # Errors in Variables: total precision (inverse of variance) of the total error.
# TODO: Check these are dimensioned correctly
par$log_ft <- matrix(0, nrow=nyrs, ncol=nfleet)
par$init_log_rec_devs <- rep(0,(dat$nage-dat$sage + 1)) #vector(length=dat$nage - dat$sage + 1) # I think this is length nage-sage+1 (i.e., length 2:9)
par$log_rec_devs <- rep(0,nyrs) #vector(length=nyrs)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. Model
# Eventually move functions to separate R scripts
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

model <- function(par){

  getAll(par,dat, pfc)
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

  # MOVE SOME OF THESE SECTIONS INTO SEPARATE FUNCTIONS

#|---------------------------------------------------------------------|
  # 1. initParameters
  ro        <- exp(par$log_ro)
  steepness <- par$h
  m         <- exp(par$log_m)
  rho       <- par$rho
  varphi    <- sqrt(1.0/par$log_kappa)
  sig       <- sqrt(rho)*varphi  # elem_prod(sqrt(rho) , varphi)
  tau       <- sqrt(1.0-rho)*varphi # elem_prod(sqrt(1.0-rho) , varphi)
  Ft        <- exp(log_ft)

  # A decision was made in 2018 to fix these two parameters to be the same
  #   as ro (P. Starr). May revisit this later
  log_avgrec <- par$log_ro
  #log_recinit <- par$log_ro # RF: I don't think it makes sense to use this one


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~ TESTING VALUES FROM data-raw/iscam.rep and data-raw/iscam.par (2020 PCOD RESULTS) ~
  # DELETE THIS ONCE MODEL EQUATIONS ARE TESTED
  ro        <- 3376.55
  steepness <- 0.809901
  m         <- 0.306306
  rho       <- 0.058824
  varphi    <- 0.824621
  log_avgrec <- log(ro)
  Ft <- c(0.100901, 0.162811, 0.218133, 0.172062, 0.135691, 0.0919969, 0.139881, 0.0772674, 0.129473, 0.191035, 0.236091, 0.209665, 0.18531, 0.149684, 0.0791808, 0.134725, 0.135795, 0.114886, 0.152902, 0.178127, 0.167726, 0.138626, 0.128951, 0.25108, 0.210596, 0.145312, 0.137206, 0.120306, 0.104469, 0.0661442, 0.16161, 0.360395, 0.24175, 0.141831, 0.152064, 0.339467, 0.306043, 0.273561, 0.121854, 0.100894, 0.111752, 0.12076, 0.1013, 0.0790893, 0.0660031, 0.0333779, 0.037986, 0.0554713, 0.0657476, 0.0801039, 0.0748695, 0.0352134, 0.0291054, 0.0524231, 0.095605, 0.0769179, 0.0558788, 0.0562869, 0.0645536, 0.0709053, 0.0428528, 0.0288908, 0.0189564, 0.0296574, 0.0258504)
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
  #          specificed by user.
  Mt <- rep(m,nyrs) # natural mortality
  Zt <- Ft+Mt # total mortality
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
  log_rt <- vector(length=nage-sage + nyrs) # log recruits for entire time series including init period

  # Initial Conditions
  snat <- exp(-m) # natural survival rate
  # Eqm mean weight
  wbar <- (snat*alpha.g + wk*(1-snat))/(1-snat*rho.g) # RF: calculation checked against rep file

  # Unfished numbers and biomass
  no <- ro/(1-snat)
  bo <- no*wbar # RF: calculation checked against rep file

  # Parameters of Stock-Recruit relationship
  # Maximum juvenile survival rate (same for BH and Ricker)
  alpha.sr <- (CR*ro)/bo # NOTE: *Steve called this so*. RF thinks this is confusing notation because the o implies unfished, whereas this is max juv survival rate when the biomass is close to zero

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
    # Then decay each age according to M
    nAge <- vector(length=length(sage:nage))
    for(j in 1:length(nAge)){
      nAge[j] <- exp(log_avgrec + init_log_rec_devs[j])*exp(-Mt[1]*(j-1))
    }
    nAge[length(nAge)] <- nAge[length(nAge)]/(1 - exp(-Mt[1])) # plus group

    ## TEST (from iscam calcNumbersAtAge)
    #tr <- vector(length(sage:nage))
    #tr(sage)        = ( log_avgrec(ih)+log_rec_devs(ih)(syr));
    #tr(sage+1,nage) = (log_recinit(ih)+init_log_rec_devs(ih));
    #tr(sage+1,nage) = tr(sage+1,nage)+log(lx(sage+1,nage));



    # Add up the numbers at age to get numbers in year 1
    numbers[1] <- sum(nAge)

    log_rt[1] <- log_avgrec+log_rec_devs[1] # Check: this should be same as log(nAge[1]) but only if


    # RFUpdate Correction: below biomass is the sum of weight at age x numbers at age not wbar
    biomass[1] = sum(elem_prod(tmp_N,d3_wt_avg[1]))
    annual_mean_wt[1] = biomass(ig,syr)/numbers[1]

  }
  if(ctl$misc[5]==1){
    # Unfished equilibrium

  }
  if(ctl$misc[5]==2){
    # Fished at equilibrium with Ft in first year

  }
  if(!ctl$misc[5] %in% 0:2){
    stop("Starting conditions not specified! ctl$misc[5] must be set to 0, 1 or 2 (see ?pcod2020ctl for definitions).")
  }



  # Time dynamics



  # Add recruitment for projection year ... assume it is average
  rnplus=exp(log_avgrec)

 # End calcNumbersBiomass_deldiff
#|---------------------------------------------------------------------|


#|---------------------------------------------------------------------|
  # 4. calcFisheryObservations_deldiff()


  # End calcFisheryObservations_deldiff
#|---------------------------------------------------------------------|


#|---------------------------------------------------------------------|
  # 5. calcSurveyObservations_deldiff()


  # End calcSurveyObservations_deldiff
#|---------------------------------------------------------------------|


#|---------------------------------------------------------------------|
  # 6. calcStockRecruitment_deldiff()

  # End calcStockRecruitment_deldiff
#|---------------------------------------------------------------------|

#|---------------------------------------------------------------------|
  # 7. calcAnnualMeanWeight_deldiff() //RF added this for P cod - only gets added to objective function if cntrl(15)==1


  # End calcAnnualMeanWeight_deldiff
#|---------------------------------------------------------------------|

#|---------------------------------------------------------------------|
# calcObjectiveFunction();



 jnll <- jnll + pll
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. Parameters for priors
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


