# This is a development script for a delay-difference model in rtmb
# The first draft attempts to reproduce the iscam delay difference model from
#   the 2018 Pacific Cod stock assessment for Haida Gwaii/Queen Charlotte
#   Sound stock in Pacific Canada (Area 5ABCD) [published in 2020]

# **For this first version, this is a ONE group, ONE area, ONE sex model
# Therefore, we can strip out a lot of the indexing in the iscam inputs**

# Package name: DDRTMB

# Authors: Robyn Forrest and Catarina Wor (Pacific Biological Station)

# Date this script started: May 15, 2024
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
#
#
#

# document and build package (these are also buttons in Rstudio)
# this will incorporate new functions saved in the R and data folders
devtools::document()
devtools::load_all()

library(tidyverse)

# Set up data and parameters
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
dat$alloc   <- rawdat$alloc # Catch allocation: 1 if commercial, 0 if survey
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
dat$nctobs  <- rawdat$nctobs # Number of catch observations
dat$catch  <- rawdat$catch # Catch observations (tonnes)




# Clean up unnecessary data columns
dat$catch  <- dat$catch %>%
  select(year, gear, value)
