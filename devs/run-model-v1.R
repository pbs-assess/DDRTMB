# This is a development script for a delay-difference model in rtmb
# The first draft attempts to reproduce the iscam delay difference model from
#   the 2020 Pacific Cod stock assessment for Haida Gwaii/Queen Charlotte
#   Sound stock in Pacific Canada (Area 5ABCD) [published in 2021]

# **For this first version, this is a ONE group, ONE area, ONE sex model**
# M is fixed, not time-varying - adapt later for tvm

# Package name: DDRTMB (not a package yet!)

# Authors: Robyn Forrest (RF), Catarina Wor (CW), Sean Anderson (SA) (Pacific Biological Station, Nanaimo, Canada)

# Date created:  May 15, 2024
# Last Modified: August 1, 2024

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
#  - Move some of model sections (R\model.R) into separate functions:
#      - Typically, move components where there is a choice to be made
#         into separate functions (e.g., recruitment)
#      - Possibly move the major model steps into separate functions,
#         but be careful that the model steps are preserved in order
#         i.e., will still need the model function to be called (this is what
#          is optimised), rather than letting users build their own
#  - Implement MCMC - tmbstan
#  - Check source of fished equilibrium equation and check code - or remove it (in calcNumbersBiomass_deldiff)
# Think about a consistent way to move objects around.
#  Some are in global space and some come in through the getAll func

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
library(purrr)
library(RTMB)
library(tmbstan)

# eventually move to standard R statistical functions
# Currently using facsimiles of the needed functions from ADMB statsLib.h
source(here("R/likelihood_funcs.R"))

# The model function is in a separate file
# There is a bunch of stuff in the global space that it needs
source(here("R/model.R"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  ~ SETTINGS ~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nsample <- 5000 # number of posterior samples for the mcmc
nchain <- # number of chains for the mcmc (outputs only look at one chain for now)
proj_years <- 1 # How many projection years for decision table
# Set test mode for testing model with MPD estimates from pcod2020 test file
# Delete this eventually. The test code is currently commented out anyway.
test <- FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

# TODO: Everything being put in the global space here should be
#  inside the model function and should come from dat (via the getAll function)
#  i.e., the model function needs a DATA_SECTION
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
saveRDS(pl, here("outputs","parameter_estimates.rda"))
saveRDS(plsd, here("outputs","parameter_sds.rda"))
saveRDS(plrad, here("outputs","derived_estimates.rda"))
saveRDS(plradsd, here("outputs","derived_sds.rda"))

# Plot results and comparisons with iscam
# Delete this for package
source(here("devs","plot_iscam_compare_mpd.r"))
source(here("devs","plot_rtmb_results_mpd.r"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. MCMCs - Posterior parameter estimates
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Use tmbstan to run MCMC
# See example code https://github.com/kaskr/tmbstan

# SA suggested 1000 samples with 4 chains as a start
# parallel causing errors - I'm not setting it up correctly
# cores <- parallel::detectCores()-1
# options(mc.cores = cores)
fitmcmc <- tmbstan(obj, chains=nchain,
                   iter=nsample,
                   init=list(opt$par),
                   lower=c(1.,0.2,-2.3,
                           rep(-5,length(par$log_ft_pars)),
                           rep(-5, length(par$log_rec_devs)),
                           rep(-5, length(par$init_log_rec_devs))),
                   upper=c(12.,1.,0.,
                           rep(5,length(par$log_ft_pars)),
                           rep(5, length(par$log_rec_devs)),
                           rep(5, length(par$init_log_rec_devs))))

# Remove burn in (warmup)
mc <- extract(fitmcmc, pars=names(obj$par),
              inc_warmup=FALSE, permuted=FALSE)

## Can also get ESS and Rhat from rstan::monitor
# https://github.com/kaskr/tmbstan
mon <- monitor(mc)
max(mon$Rhat)
min(mon$Tail_ESS)

# Save results as a data frame
# (if running more than one chain, then would do this for each chain)
mc.df <- as.data.frame(mc[,1,])
saveRDS(mc.df, here("outputs","MCMC_parameter_estimates.rda"))
saveRDS(mon, here("outputs","MCMC_diagnostics.rda"))

# # Four chains
# fitmcmc_4ch <- tmbstan(obj, chains=4,
#                    iter=Iter,
#                    init=list(opt$par),
#                    lower=c(1.,0.2,-2.3,
#                            rep(-5,length(par$log_ft_pars)),
#                            rep(-5, length(par$log_rec_devs)),
#                            rep(-5, length(par$init_log_rec_devs))),
#                    upper=c(12.,1.,0.,
#                            rep(5,length(par$log_ft_pars)),
#                            rep(5, length(par$log_rec_devs)),
#                            rep(5, length(par$init_log_rec_devs))))
# mc4ch <- extract(fitmcmc_4ch, pars=names(obj$par),
#               inc_warmup=FALSE, permuted=FALSE)
#
# ## Can also get ESS and Rhat from rstan::monitor
# # https://github.com/kaskr/tmbstan
# mon4ch <- monitor(mc4ch)
#
# mc.df.4ch <- list()
# for(i in 1:4){
#   mc.df.4ch[[i]] <- as.data.frame(mc4ch[,i,])
# }
# saveRDS(mc.df.4ch, here("outputs","MCMCParameterEstimates_4chain_list.rda"))
# saveRDS(mon4ch, here("outputs","MCMCDiagnostics_4chain_list.rda"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5. Get posteriors for derived variables (REPORT objects)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Rerun model with posterior parameter estimates (extract REPORT objects)
#  See readme at https://github.com/kaskr/tmbstan:
#    "What if you want a posterior for derived quantities in the report? Just
#    loop through each posterior sample (row) and call the report function"

# Version 1: for reporting, graphs etc
# This is a list where each element is a variable of interest
# Make a list for putting posterior estimates of biomass,
#    recruits and q
# Note that parameters, rec devs and logf are
#    already reported in the mc object
#   (but have added Ft to REPORT to simplify projection model)
posteriors_by_variable <- list()
posteriors_by_variable$biomass  <- matrix(NA, ncol=nyrs+1, nrow=nrow(post))
posteriors_by_variable$numbers <- matrix(NA, ncol=nyrs+1, nrow=nrow(post))
posteriors_by_variable$recruits <- matrix(NA, ncol=nyrs-dat$sage, nrow=nrow(post))
posteriors_by_variable$S <- matrix(NA, ncol=nyrs, nrow=nrow(post))
posteriors_by_variable$Ft <- matrix(NA, ncol=nyrs, nrow=nrow(post))
posteriors_by_variable$q <- matrix(NA, ncol=dat$nit, nrow=nrow(post))

# Version 2: for projections
# This is a list where each element is a posterior sample
#  containing all the variables needed for the proj model
#  Gets passed to projection model using purrr::map2_df()
posteriors_by_sample <- list()

# Get the posterior output from tmbstan, as matrix
post <- as.matrix(mc.df)

for(i in 1:nrow(post)){
  r <- obj$report(post[i,])

  # Posteriors by variable
  posteriors_by_variable$biomass[i,]  <- r$biomass
  posteriors_by_variable$numbers[i,]  <- r$numbers
  posteriors_by_variable$recruits[i,] <- r$rt
  posteriors_by_variable$surv[i,]  <- r$surv
  posteriors_by_variable$Ft[i,] <- r$Ft
  posteriors_by_variable$q[i,]  <- r$q

  # Posteriors by sample (for project_model)
  posteriors_by_sample[[i]] <- r
  # Append 3 leading parameters and proj_years
  posteriors_by_sample[[i]]$log_ro <- mc.df$log_ro[i]
  posteriors_by_sample[[i]]$h      <- mc.df$h[i]
  posteriors_by_sample[[i]]$log_m  <- mc.df$log_m[i]
  posteriors_by_sample[[i]]$proj_years  <- proj_years # for now add projection years here
}

saveRDS(posteriors_by_variable, here("outputs","MCMC_derived_estimates.rda"))
saveRDS(posteriors_by_sample, here("outputs","MCMC_outputs_bysample.rda"))

# Plot results and comparisons with iscam
# Delete this for package
source(here("devs","plot_iscam_compare_mcmc.r"))
source(here("devs","plot_rtmb_results_mcmc.r"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5. Projections - post MCMC step
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# TODO: THIS SHOULD BE A FUNCTION

# NEED:
# 1. Posterior samples from the historical period (syr:nyr),
#     plus the leading parameters log_ro, h and log_m
#   i) Biomass, Numbers, Fishing mortality, Survival rate
#        and Recruits from the historical period
#        (from posteriors object)
# 2. A vector of future TACs (from pfc file)
# 3. Number of projection years

# try putting outputs in a list:
#   one list item for each tac (from pfc$tac.vec)
# each list item is a matrix with dim
#    nrow=number of posterior samples,
#    ncol=however many variables we want for the decision tables
# NOTE: iscam had a giant matrix. It was slightly challenging to
#       get the results out for the decision tables.
#       Having one list object per tac will be easier for getting
#       out the probabilities

# List object for projection outputs
# This will be the inputs for decision tables
proj_out <- list()

# Read in posterior samples (in case just doing projections)
posteriors_by_sample <- readRDS(here("outputs","MCMC_outputs_bysample.rda"))

# Need to loop over future TACs but do not need to loop
#  over posterior samples. Let purrr do that.
for(i in 1:5){
  tac <- pfc$tac.vec[i]

  # Run the projection model for tac[i]
  proj_out[[i]] <- purrr::map2_df(posteriors_by_sample, tac, project_model)
}
names(proj_out) <- pfc$tac.vec[1:5]

# TODO:
# Add decision table code
# Produce decision tables and graphics

