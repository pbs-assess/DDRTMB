# This is a development script for a delay-difference model in rtmb
# The first draft attempts to reproduce the iscam delay difference model from
#   the 2020 Pacific Cod stock assessment for Haida Gwaii/Queen Charlotte
#   Sound stock in Pacific Canada (Area 5ABCD) [published in 2021]

# **For this first version, this is a ONE group, ONE area, ONE sex model**
# M is fixed, not time-varying - adapt later for tvm

# Package name: DDRTMB (not a package yet!)

# Authors: Robyn Forrest (RF), Catarina Wor (CW), Sean Anderson (SA) (Pacific Biological Station, Nanaimo, Canada)

# Date created:  May 15, 2024
# Last Modified: July 25, 2024

# Notes:
# - The iscam input files are already loaded into the package:
#   1. pcod2020dat (the input data)
#   2. pcod2020ctl (controls for running the model,
#      including prior descriptions and misc settings)
#   3. pcod2020pfc (controls for projections and decision tables)
#
# - Iscam uses an errors-in-variables approach to partition observation error,
#       with multiplicative weighting of each index observation using annual CVs
# - Possibly want to change this to additive weightings (as per SS3), then also explore
#       state space options. But first try to reproduce the iscam results.

# TODO (move these to Issues on gitHub repo):
# 1. Translation of iscam to RTMB:
#  - Move some of model sections into separate functions
#  - Implement MCMC - tmbstan
#  - Check source of fished equilibrium equation and check code - or remove it (in calcNumbersBiomass_deldiff)

# 2. Potential model changes:
#  - Tidy up the three recruitment parameters - currently set to all be the same as per Paul Starr's request
#  - Need for Jacobian transformations?
#  - Look at bias correction (see Thorson and Kristensen paper - but not needed for Bayesian)
#  - Look at alternate settings for Errors in Variables (e.g., weights additive instead of multiplicative)
#  - Look at state-space implementation
#  - Add time-varying M
#  - Think about how to make this a multi-species, multi-area model? (MICE)
#  - Some of the architecture is already in iscam
#  - consider a stan version if this doesn't work well
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

# eventually move to standard R statistical functions
# Currently using facsimiles of the needed functions from ADMB statsLib.h
source(here("R/likelihood_funcs.R"))

# The model function is in a separate file
# There is a bunch of stuff in the global space that it needs
source(here("R/model.R"))

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
ctl <- pcod2020ctl # Control inputs. Use ?pcod2020ctl to see definitions. Not all used in d-d model
pfc <- pcod2020pfc # Control inputs for projections. Use ?pcod2020pfc to see definitions
nyrs <- dat$nyr-dat$syr+1
yrs  <-  dat$syr:dat$nyr
ages <-  dat$sage:dat$nage

# **Maybe add all of these objects to dat**
# get the number of and index for commercial (fishery) fleets
nfleet <- 0 # number of fishing fleets (not surveys)
for(i in 1:dat$ngear){
  nfleet <- ifelse(dat$alloc[i]>0, nfleet+1, nfleet) # not sure if this is needed
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
d3_wt_avg <- wa # just to be consistent with iscam rep file

# Get settings for priors
# Leading parameters
num_params <- ctl$num.params
theta_control <- ctl$params %>%  as.data.frame() # settings for initialization and priors
q_control <- t(ctl$surv.q) %>%  as.data.frame() # settings for q priors (transposed from ctl file)

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
par$log_m <- theta_control[3,1] #-1.18317 # (rep file estimate) theta_control[3,1] # log natural mortality
# In the current version of the P cod model, these are set the same as ro so do not estimate
#par$log_avgrec <- theta_control[4,1] # log average recruitment (syr+1 to nyr)
#par$log_recinit <- ptheta_control[5,1] #l og average of initial recruitments to fill first year if population if population is not unfished at syr

# Variance parameters are fixed for P cod so should not be in par
par$rho <- theta_control[6,1] # Errors in Variables: fraction of the total variance associated with observation error
par$kappa <- theta_control[7,1] # Errors in Variables: total precision (inverse of variance) of the total error.
# TODO: Check these are dimensioned and initialized correctly
par$log_ft_pars <- numeric(nyrs) # estimated log fishing mortalities (total across fleets)
par$init_log_rec_devs <- numeric(length=length((dat$sage+1):dat$nage)) # I think this is length nage-sage+1 (i.e., length 2:9)
par$log_rec_devs <- numeric(nyrs)
par

# Test obj function: iscam has 195.806
# Current test with pars fixed from iscam.rep: 195.806  :-)
# Current RTMB value: 195.806  :-) :-) :-)
#model(par)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. Estimate parameters with RTMB
# PROCEDURE_SECTION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## MakeADFun builds the graph, basically "compiles" the model with random effects identified
## from TMB help: map = List defining how to optionally collect and fix parameters
## Means you can fix some instances of a vector/matrix of parameters and estimate ones with the same factor id to be the same
# Fixing rho and kappa  log_m=factor(NA) h=factor(NA),
# Note the model is in model.R and was sourced above
obj <- MakeADFun(model, par, silent=FALSE,
               map=list(rho=factor(NA), kappa=factor(NA)))
# The optimization step - gets passed the parameters, likelihood function and gradients Makeby ADFun
opt <- nlminb(obj$par, obj$fn, obj$gr, control=list(eval.max=1000, iter.max=1000))
opt$objective
#sdr <- summary(sdreport(obj))

# Estimated MPD parameters
pl <- as.list(sdreport(obj),"Est")
plsd <- as.list(sdreport(obj),"Std")
# Derived quantities with standard errors (from ADREPORT)
plrad <- as.list(sdreport(obj),"Est", report=TRUE)
plradsd <- as.list(sdreport(obj),"Std", report=TRUE)
# Derived quantities without standard errors (from REPORT)
# get them out like this: obj$report()$X
plr <- as.list(obj$report())

# Write out MPD results for plotting
saveRDS(pl, here("outputs","ParameterEstimates.rda"))
saveRDS(plsd, here("outputs","ParameterSDs.rda"))
saveRDS(plrad, here("outputs","DerivedEstimates.rda"))
saveRDS(plradsd, here("outputs","DerivedSDs.rda"))

# Plot comparisons with iscam
# Delete this for package
source(here("devs","plot_iscam_compare_mpd.r"))
source(here("devs","plot_rtmb_results_mpd.r"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. MCMCs - Posterior parameter estimates
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Use tmbstan to run MCMC
# See example code https://github.com/kaskr/tmbstan
library(tmbstan)
# SA suggested 1000 samples with 4 chains as a start
fitmcmc <- tmbstan(obj, chains=2,
                   iter=1000,
                   init=list(opt$par),
                   lower=c(1.,0.2,-2.3,
                           rep(-5,length(par$log_ft_pars)),
                           rep(-5, length(par$log_rec_devs)),
                           rep(-5, length(par$init_log_rec_devs))),
                   upper=c(12.,1.,0.,
                           rep(5,length(par$log_ft_pars)),
                           rep(5, length(par$log_rec_devs)),
                           rep(5, length(par$init_log_rec_devs))))

mc <- extract(fitmcmc, pars=names(obj$par),
              inc_warmup=TRUE, permuted=FALSE)

## Can also get ESS and Rhat from rstan::monitor
# https://github.com/kaskr/tmbstan
mon <- monitor(mc)
max(mon$Rhat)
min(mon$Tail_ESS)

mc.df <- as.data.frame(mc[,1,])
hist(mc.df$log_ro)
hist(mc.df$h)
hist(exp(mc.df$log_m))

saveRDS(mc.df, here("outputs","MCMCParameterEstimates.rda"))
saveRDS(mon, here("outputs","MCMCDiagnostics.rda"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5. Posterior derived parameters and model variables
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Rerun model with posterior parameter estimates,
# Extract REPORT objects

# https://github.com/kaskr/tmbstan
## What if you want a posterior for derived quantities in the report? Just
## loop through each posterior sample (row) and call the report function
## which returns a list. The last column is the log-posterior density (lp__)
## and needs to be dropped
#post <- as.matrix(mc) # this just returns one massive column of the posterior samples
post <- as.matrix(mc[,1,]) # I think this is what I want
#obj$report(post[1,-ncol(post)])
obj$report(post[1,]) # the last column is the final year of rec devs so don't drop it

posteriors <- list()
posteriors$biomass <- matrix(NA, ncol=nyrs+1, nrow=nrow(post))
posteriors$recruits <- matrix(NA, ncol=nyrs+1, nrow=nrow(post))

for(i in 1:nrow(post)){
  r <- obj$report(post[i,])
  post_biomass[i,] <- r$biomass
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5. Diagnostics
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



