# Simple plots to compare MPD estimates of results
# Robyn Forrest
# Date created:  July 23, 2024
# Last Modified: July 23, 2024
#
# FIRST run the model in run-model-v1.R
# Load documentation and inputs
devtools::document()
devtools::load_all()

library(tidyverse)
library(reshape2)
library(cowplot)
library(here)
source("devs/load-models.R")

# 1. Read in output files from RTMB and iscam
dat <- pcod2020dat
ctl <- pcod2020ctl
pl <- readRDS(here("outputs","ParameterEstimates.rda"))
plsd <- readRDS(here("outputs","ParameterSDs.rda"))
plr <- readRDS(here("outputs","DerivedEstimates.rda"))
plrsd <- readRDS(here("outputs","DerivedSDs.rda"))
pcod2020rep <- read.report.file("data-raw/iscam.rep")

yrs <-  dat$syr:dat$nyr
pyrs <- dat$syr:(dat$nyr+1) # includes projection year

# Basic plots
# Biomass
maxX <- max(c(pcod2020rep$biomass,(plr$biomass+2*plrsd$biomass)))
plot(pyrs, plr$biomass, type="l", lwd=3, col=1, ylim=c(0,maxX), xlab="Years", ylab="Biomass (t)")
lines(pyrs, plr$biomass-2*plrsd$biomass, lwd=3, col=1, lty="dotted")
lines(pyrs, plr$biomass+2*plrsd$biomass, lwd=3, col=1, lty="dotted")
lines(pyrs, pcod2020rep$biomass, lwd=3, col="darkred")
legend("topright", legend=c("RTMB", "iscam MPD"), lwd=3, col=c(1,"darkred") ,bty="n")


