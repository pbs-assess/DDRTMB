# Read in Pacific Cod 2020 assessment inputs (Stock Area 5ABCD)
# Strip out any inputs that only relate to the iscam age-structured model
# For this first version, this is a ONE group, ONE area, ONE sex model
# Therefore, we can also strip out a lot of the indexing in the iscam inputs

# ONLY DO THIS ONCE. THE LAST LINE ADDS THE FILTERED INPUTS TO THE PACKAGE #

# Authors: Robyn Forrest and Catarina Wor (Pacific Biological Station)

# Date created: May 15, 2024
# Last Modified: May 21, 2024

# Document and build package (these are also buttons in Rstudio)
#    this will incorporate new functions saved in the R and data folders
devtools::document()
devtools::load_all()

library(tidyverse)
source("devs/load-models.R")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set up data and parameter controls. Better to make this a function, or even just
#  make it part of the package.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in raw input files
rawpcod2020dat<-read.data.file("data-raw/pcod.dat")
rawpcod2020ctl<-read.control.file("data-raw/pcod.ctl",
                               num.gears =6,
                               num.age.gears = 1,)
rawpcod2020pfc<-read.projection.file("data-raw/pcod.pfc")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Strip out anything not needed for this d-d model
pcod2020dat <- list()
pcod2020dat$gearNames <- rawpcod2020dat$gearNames # names of the fleets
# Index markers
pcod2020dat$syr     <- rawpcod2020dat$syr # Start year
pcod2020dat$nyr     <- rawpcod2020dat$nyr # Final year
pcod2020dat$sage    <- rawpcod2020dat$sage # Start age
pcod2020dat$nage    <- rawpcod2020dat$nage # Max age
pcod2020dat$ngear   <- rawpcod2020dat$ngear # Number of gears (includes commercial and survey)
pcod2020dat$alloc   <- rawpcod2020dat$alloc # CHECK: Catch allocation: 0 < alloc <= 1 if commercial, 0 if survey
# Fixed parameters
pcod2020dat$linf    <- rawpcod2020dat$linf # Growth linfinity (not sure if used in d-d, has its own growth pars)
pcod2020dat$k       <- rawpcod2020dat$k # Growth vonB K (not sure if used in d-d, has its own growth pars)
pcod2020dat$to      <- rawpcod2020dat$to # Growth t_0 (not sure if used in d-d, has its own growth pars)
pcod2020dat$lwscal  <- rawpcod2020dat$lwscal # Growth L-W a (not sure if used in d-d)
pcod2020dat$lwpow   <- rawpcod2020dat$lwpow # Growth L-W b (not sure if used in d-d)
pcod2020dat$kage    <- rawpcod2020dat$dd.kage # Knife-edge age at recruitment
pcod2020dat$alpha.g <- rawpcod2020dat$dd.alpha.g # Growth Ford-Walford alpha
pcod2020dat$rho.g   <- rawpcod2020dat$dd.rho.g # Growth Ford-Walford rho
pcod2020dat$wk      <- rawpcod2020dat$wk # Weight at age of recruitment
# Catch and Index observations
# Catch
# Unneeded columns removed below. See below for column descriptions.
pcod2020dat$nctobs  <- rawpcod2020dat$nctobs # Number of catch observations
pcod2020dat$catch   <- rawpcod2020dat$catch # Catch observations (tonnes).
# Survey indices
# Unneeded columns removed below. See below for column descriptions.
pcod2020dat$nit     <- rawpcod2020dat$nit # Number of surveys
pcod2020dat$nitnobs <- rawpcod2020dat$nitnobs # Number of survey obs (a vector of length pcod2020dat$nit)
pcod2020dat$survtype <- rawpcod2020dat$survtype # Survey type (a vector of length pcod2020dat$nit) - see below
pcod2020dat$indices <- rawpcod2020dat$indices # a list of dataframes (n list elements = length pcod2020dat$nit)
# Annual commercial mean weight data
# Unneeded columns removed below. See below for column descriptions.
pcod2020dat$nmeanwt <- rawpcod2020dat$nmeanwt # Number of mean weight series (if none, must be set to 1 - CAN NOW CHANGE THIS REQUIREMENT)
pcod2020dat$nmeanwtobs <- rawpcod2020dat$nmeanwtobs # Number of mean weight observations
pcod2020dat$meanwtdata <- rawpcod2020dat$meanwtdata # Annual commercial mean weights (kg)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Clean up unnecessary data columns
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Catch data: a matrix of pcod2020dat$nctobs x 4
# Type: 1 = catch in weight; 2 = catch in numbers
# gear is gear index from pcod2020dat$ngear, which includes both catch and surveys
# If no data, just skip the row for that year
pcod2020dat$catch  <- pcod2020dat$catch %>%
  as.data.frame() %>%
  select(year, gear, type, value)

# Index data: a list with pcod2020dat$nit elements, each is a matrix of 1:pcod2020dat$nitnobs[i] x 5
# it: Index value. Set type in pcod2020dat$survtype:
#   1 = survey is proportional to vulnerable numbers
#   2 = survey is proportional to vulnerable biomass
#   3 = survey is proportional to spawning biomass (e.g., a spawn survey)
# gear: gear index from pcod2020dat$ngear, which includes both catch and surveys.
#      So if there is one commercial gear then all the survey gears will be 1:pcod2020dat$nit + 1
# wt: 1/CV (used as precision to multiplicatively weight observations).
#     Relative weights for each relative abundance normalized to have a
#     mean = 1 so rho = sig^2/(sig^2+tau^2) holds true in variance pars.
# timing: The fraction of total mortality that has occurred prior to survey. Usually zero.
# If no data, just skip the row for that year
for(i in 1:pcod2020dat$nit){
  pcod2020dat$indices[[i]] <- pcod2020dat$indices[[i]] %>%
    as.data.frame() %>%
    select(iyr, it, gear, wt, timing)
}

# Annual mean weight data: a matrix of pcod2020dat$nmeanwtobs x 4
# Timing should match survey (usually 0)
# If no data, just skip the row for that year
pcod2020dat$meanwtdata  <- pcod2020dat$meanwtdata %>%
  as.data.frame() %>%
  select(year, meanwt, gear, timing)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Parameter controls
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pcod2020ctl <- list()
# Strip out anything not needed for this d-d model
pcod2020ctl$num.params <- rawpcod2020ctl$num.params # Number of leading parameters to be estimated with priors
#                                  (does not include nuisance pars like q)
pcod2020ctl$params <- rawpcod2020ctl$params # A matrix rawpcod2020ctl$num.params x 7 (with row names)
# Some of these settings will be redundant in tmb (e.g., lb, ub, phz)?
# NEED to convert some of these to log space
# ROWS
# 1. log_ro = log unfished recruitment
# 2. steepness = steepness
# 3. log_m = log natural Mortality
# 4. log_avrec = log average recruitment (set the same as ro)
# 5. log_recinit = log average of initial recruitments to fill first year (set the same as ro)
# 6. rho = EIV: fraction of the total variance associated with observation error
# 7. kappa (varphi) = EIV: total precision (inverse of variance) of the total error
# NOTE: For P. cod VARPHI AND RHO ARE FIXED TO make average ratio of sds from the two surveys the same as 2005

# ival = starting value
# lb = lower bound (not sure if this can be used in tmb)
# ub = upper bound (not sure if this can be used in tmb)
# phz = phase of estimation (in ADMB). Won't need this in tmb but need to check how to fix a par
# prior = prior type:
#   -0 uniform      (0,0)
#   -1 normal       (p1=mu,p2=sig)
#   -2 lognormal    (p1=log(mu),p2=sig)
#   -3 beta         (p1=alpha,p2=beta)
#   -4 gamma        (p1=alpha,p2=beta)
# p1, p2 = prior distribution parameters, see above

# Priors for q
pcod2020ctl$num.indices <- pcod2020dat$nit # this was a separate parameter in iscam in the ctl file - don't use this one, but keep for now
pcod2020ctl$surv.q      <- rawpcod2020ctl$surv.q # parameters for priors on q. a matrix 3 x pcod2020dat$nit with row names
# ROWS:
# priortype - see below
# priormeanlog = mean of prior in log space
# priorsd = sd of prior
# Prior type:
#       0 - uninformative prior
#       1 - normal prior density for log(q)
#       2 - random walk in q
## Need one column for each survey.

# Controls for fitting mean weight data
pcod2020ctl$fit.mean.weight <- rawpcod2020ctl$fit.mean.weight # 1 = fit to annual mean weights, 0 = do not fit to annual mean weights
pcod2020ctl$num.mean.weight <- rawpcod2020ctl$num.mean.weight.cv # Number of annual mean weight series
pcod2020ctl$weight.sig <- rawpcod2020ctl$weight.sig # SD for likelihood for fitting to annual mean weight (one for each series)

# Miscellaneous controls
pcod2020ctl$misc      <- rawpcod2020ctl$misc[1:13] # A matrix 13 x 1 with row names
# 1  -verbose ADMB output (0=off, 1=on)
# 2  -recruitment model (1=beverton-holt, 2=ricker)
# 3  -std in observed catches in first phase.
# 4  -std in observed catches in last phase.
# 5  -Assume unfished equilibrium in first year (0=FALSE, 1=TRUE, 2 = AT EQUILIBRIUM WITH FISHING MORTALITY IN SYR - IMPLEMENTED ONLY IN DELAY DIFF MODEL)
# 6  -Maternal effects multiplier
# 7  -Mean fishing mortality for regularizing the estimates of Ft
# 8  -std in mean fishing mortality in first phase
# 9  -std in mean fishing mortality in last phase
# 10 -phase for estimating m_deviations (use -1 to turn off mdevs)
# 11 -std in deviations for natural mortality
# 12 -number of estimated nodes for deviations in natural mortality
# 13 -fraction of total mortality that takes place prior to spawning

# Projection control file
pcod2020pfc <- list()
pcod2020pfc$num.tac <- rawpcod2020pfc$num.tac # Number of TAC options for decision table projections
pcod2020pfc$tac.vec <- rawpcod2020pfc$tac.vec # TAC options for decision table projections
pcod2020pfc$num.pcod2020ctl.options <- rawpcod2020pfc$num.pcod2020ctl.options # Number of options in pcod2020ctl.options
pcod2020pfc$pcod2020ctl.options <- rawpcod2020pfc$pcod2020ctl.options # options for projections: Matrix 1 x 9
## - 1) Start year for mean natural mortality rate
## - 2)  Last year for mean natural mortality rate

## - 3) Start year for average fecundity/weight-at-age in projections. CHECK IF USED IN D-D
## - 4)  Last year for average fecundity/weight-at-age in projections. CHECK IF USED IN D-D

## - 5) Start year for average recruitment period in projections
## - 6) End year for average recruitment period in projections

## - 7)Short time series: "historical" control points based on biomass and F reconstruction
## - 8) Long time series:"historical" control points based on biomass and F reconstruction

## - 9) bmin for "minimum biomass from which the stock recovered to above average" for "historical" control points based on biomass and F reconstruction

#save these to a data folder
usethis::use_data(pcod2020dat, overwrite = TRUE)
usethis::use_data(pcod2020ctl, overwrite = TRUE)
usethis::use_data(pcod2020pfc, overwrite = TRUE)

