# This is a development script for a delay-difference model in rtmb
# The first draft attempts to reproduce the iscam delay difference model from
#   the 2018 Pacific Cod stock assessment for Haida Gwaii/Queen Charlotte
#   Sound stock in Pacific Canada (Area 5ABCD) [published in 2020]

# **For this first version, this is a ONE group, ONE area, ONE sex model
# Therefore, we can strip out a lot of the indexing in the iscam inputs**
# M is fixed, not time-varying - adapt later for tvm

# Package name: DDRTMB

# Authors: Robyn Forrest and Catarina Wor (Pacific Biological Station)

# Date created: May 15, 2024
# Last Modified: May 15, 2024

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
par$init_log_rec_devs <- rep(0,(dat$nage-dat$sage)) #vector(length=dat$nage - dat$sage) # I think this is length nage-sage
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
  # 1. Initiate parameters
  ro        <- exp(par$log_ro)
  steepness <- par$steepness
  m         <- exp(par$log_m)
  rho       <- par$rho
  varphi    <- sqrt(1.0/par$log_kappa)
  sig       <- sqrt(rho)*varphi  # elem_prod(sqrt(rho) , varphi)
  tau       <- sqrt(1.0-rho)*varphi # elem_prod(sqrt(1.0-rho) , varphi)
  Ft        <- exp(log_ft)

  # A decision was made in 2018 to fix these two parameters to be the same
  #   as ro (P. Starr). May revisit this later
  log_avrec <- par$log_ro
  log_recinit <- par$log_ro
#|---------------------------------------------------------------------|

#|---------------------------------------------------------------------|
  # 2. calcTotalMortality_deldiff();
  # Purpose: This function calculates fishing mortality, total mortality and annual
  #          survival rates S=exp(-Z) for each year.
  #          Z also is updated with time-varying natural mortality rates if
  #          specificed by user.
  Mt <- rep(m,nyrs) # natural mortality
  Zt <- Ft+Mt # total mortality
  surv <- exp(-Zt) # survival rate

#|---------------------------------------------------------------------|

#|---------------------------------------------------------------------|
  # 3. calcNumbersBiomass_deldiff();


#|---------------------------------------------------------------------|


#|---------------------------------------------------------------------|
  # 4. calcFisheryObservations_deldiff();


#|---------------------------------------------------------------------|


#|---------------------------------------------------------------------|
# 5. calcSurveyObservations_deldiff()


#|---------------------------------------------------------------------|


#|---------------------------------------------------------------------|
# 6. calcStockRecruitment_deldiff();


#|---------------------------------------------------------------------|

#|---------------------------------------------------------------------|
# 7. calcAnnualMeanWeight_deldiff(); //RF added this for P cod - only gets added to objective function if cntrl(15)==1


#|---------------------------------------------------------------------|

#|---------------------------------------------------------------------|
# calcObjectiveFunction();





 jnll <- jnll + pll
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


