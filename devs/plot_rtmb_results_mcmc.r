# Simple plots to look at RTMB MCMC results
# Robyn Forrest
# Date created:  July 26, 2024
# Last Modified: July 26, 2024
#
# FIRST run the model in run-model-v1.R
# Load documentation and inputs
# Uncomment these 2 lines if running standalone
devtools::document()
devtools::load_all()

library(tidyverse)
library(reshape2)
library(cowplot)
library(here)
source("devs/load-models.R")

# Settings
# Colours for data, iscam and RTMB
datcol <- 1
rtmbcol <- "darkblue"

if(!file.exists(here("outputs"))) dir.create(here("outputs"), recursive = TRUE)
if(!file.exists(here("outputs","figs"))) dir.create(here("outputs","figs"), recursive = TRUE)

# 1. Read in MCMC outputs
mcmcpars <- readRDS(here("outputs","MCMCParameterEstimates.rda"))
mcmcderived <- readRDS( here("outputs","MCMCDerivedEstimates.rda"))
mcmcdiagnostics <- readRDS(here("outputs","MCMCDiagnostics.rda"))

# Pairs and trace plots


# Parameter histograms (without and with priors)


# Biomass, recruits and fishing mortality time series



