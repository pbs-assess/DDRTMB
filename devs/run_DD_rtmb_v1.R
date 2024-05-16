# This is a development script for a delay-difference model in rtmb
# The first draft attempts to reproduce the iscam delay difference model from
#   the 2018 Pacific Cod stock assessment for Haida Gwaii/Queen Charlotte
#   Sound stock in Pacific Canada (Area 5ABCD) [published in 2020]

# **For this first version, this is a ONE group, ONE area, ONE sex model
# Therefore, we can strip out a lot of the indexing in the iscam inputs**

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

# TODO:
# CHECK how to set dat$alloc if more than one commercial gear
# Can probably relax requirement of setting nmeanwt to 1 when no mean weight data
# Put the dat, ctl and pfc stripping below into functions

# document and build package (these are also buttons in Rstudio)
# this will incorporate new functions saved in the R and data folders
devtools::document()
devtools::load_all()

library(tidyverse)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set up data and parameter controls. Maybe make this a function.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rawdat<- pcod2020dat # this is a list already loaded into the package
# Strip out anything not needed for this d-d model
dat <- list()
dat$gearNames <- rawdat$gearNames # names of the fleets
# Index markers
dat$syr     <- rawdat$syr # Start year
dat$nyr     <- rawdat$nyr # Final year
dat$sage    <- rawdat$sage # Start age
dat$nage    <- rawdat$nage # Final age
dat$ngear   <- rawdat$ngear # Number of gears (includes commercial and survey)
dat$alloc   <- rawdat$alloc # CHECK: Catch allocation: 0 < alloc <= 1 if commercial, 0 if survey
# Fixed parameters
dat$linf    <- rawdat$linf # Growth linfinity (not sure if used in d-d, has its own growth pars)
dat$k       <- rawdat$k # Growth vonB K (not sure if used in d-d, has its own growth pars)
dat$to      <- rawdat$to # Growth t_0 (not sure if used in d-d, has its own growth pars)
dat$lwscal  <- rawdat$lwscal # Growth L-W a (not sure if used in d-d)
dat$lwpow   <- rawdat$lwpow # Growth L-W b (not sure if used in d-d)
dat$kage    <- rawdat$dd.kage # Knife-edge age at recruitment
dat$alpha.g <- rawdat$dd.alpha.g # Growth Ford-Walford alpha
dat$rho.g   <- rawdat$dd.rho.g # Growth Ford-Walford rho
dat$wk      <- rawdat$wk # Weight at age of recruitment
# Catch and Index observations
# Catch
# Unneeded columns removed below. See below for column descriptions.
dat$nctobs  <- rawdat$nctobs # Number of catch observations
dat$catch   <- rawdat$catch # Catch observations (tonnes).
# Survey indices
# Unneeded columns removed below. See below for column descriptions.
dat$nit     <- rawdat$nit # Number of surveys
dat$nitnobs <- rawdat$nitnobs # Number of survey obs (a vector of length dat$nit)
dat$survtype <- rawdat$survtype # Survey type (a vector of length dat$nit) - see below
dat$indices <- rawdat$indices # a list of dataframes (n list elements = length dat$nit)
# Annual commercial mean weight data
# Unneeded columns removed below. See below for column descriptions.
dat$nmeanwt <- rawdat$nmeanwt # Number of mean weight series (if none, must be set to 1 - CAN NOW CHANGE THIS REQUIREMENT)
dat$nmeanwtobs <- rawdat$nmeanwtobs # Number of mean weight observations
dat$meanwtdata <- rawdat$meanwtdata # Annual commercial mean weights (kg)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Clean up unnecessary data columns
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Catch data: a matrix of dat$nctobs x 4
# Type: 1 = catch in numbers; 2 = catch in weight
# gear is gear index from dat$ngear, which includes both catch and surveys
# If no data, just skip the row for that year
dat$catch  <- dat$catch %>%
  as.data.frame() %>%
  select(year, gear, type, value)

# Index data: a list with dat$nit elements, each is a matrix of 1:dat$nitnobs[i] x 5
# it: Index value. Set type in dat$survtype:
#   1 = survey is proportional to vulnerable numbers
#   2 = survey is proportional to vulnerable biomass
#   3 = survey is proportional to spawning biomass (e.g., a spawn survey)
# gear: gear index from dat$ngear, which includes both catch and surveys.
#      So if there is one commercial gear then all the survey gears will be 1:dat$nit + 1
# wt: 1/CV (used as precision to multiplicatively weight observations).
#     Relative weights for each relative abundance normalized to have a
#     mean = 1 so rho = sig^2/(sig^2+tau^2) holds true in variance pars.
# timing: The fraction of total mortality that has occurred prior to survey. Usually zero.
# If no data, just skip the row for that year
for(i in 1:dat$nit){
  dat$indices[[i]] <- dat$indices[[i]] %>%
    as.data.frame() %>%
    select(iyr, it, gear, wt, timing)
}

# Annual mean weight data: a matrix of dat$nmeanwtobs x 4
# Timing should match survey (usually 0)
# If no data, just skip the row for that year
dat$meanwtdata  <- dat$meanwtdata %>%
  as.data.frame() %>%
  select(year, meanwt, gear, timing)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Parameter controls
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rawctl<- pcod2020ctl # this is a list already loaded into the package
# Strip out anything not needed for this d-d model

