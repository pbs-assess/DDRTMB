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

# Settings
datcol <- 1
iscamcol <- "red"
rtmbcol <- "darkblue"

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
# Catch
maxY <- max(dat$catch$value)
plot(yrs,pcod2020rep$ct, col=iscamcol, lwd=3,type="l", ylim=c(0,maxY), xlab="Years", ylab="Catch (t)")
lines(yrs,plr$ct, col=rtmbcol, lwd=3, lty=2)
points(yrs, dat$catch$value, type="p", pch=19, cex=1.2, col=datcol)
arrows(x0=yrs, y0=dat$catch$value-0.05*dat$catch$value, x1 = yrs, y1=dat$catch$value+0.05*dat$catch$value, code = 0)
legend("topright", legend=c("Data", "RTMB", "iscam MPD"), lwd=3, col=c(datcol, rtmbcol, iscamcol) ,bty="n")

# Indices
par(mfrow=c(2,3))
for(i in 1:dat$nit){
  obsindex <- dat$indices[[i]]$it
  iwt <- 1/(dat$indices[[i]]$wt)
  iscamindex <- pcod2020rep$it_hat[i,][which(!is.na(pcod2020rep$it_hat[1,]))]
  rtmbindex <- plr$`it_hat[[1]]`
  maxY <- max(index)
  plot(yrs,pcod2020rep$ct, col=iscamcol, lwd=3,type="l", ylim=c(0,maxY), xlab="Years", ylab="Catch (t)")
  lines(yrs,plr$ct, col=rtmbcol, lwd=3, lty=2)
  points(yrs, dat$catch$value, type="p", pch=19, cex=1.2, col=datcol)
  arrows(x0=yrs, y0=dat$catch$value-0.05*dat$catch$value, x1 = yrs, y1=dat$catch$value+0.05*dat$catch$value, code = 0)
  legend("topright", legend=c("Data", "RTMB", "iscam MPD"), lwd=3, col=c(datcol, rtmbcol, iscamcol) ,bty="n")
}

# Biomass
maxY <- max(c(pcod2020rep$biomass,(plr$biomass+2*plrsd$biomass)))
plot(pyrs, plr$biomass, type="l", lwd=3, col=rtmbcol, ylim=c(0,maxY), xlab="Years", ylab="Biomass (t)")
lines(pyrs, plr$biomass-2*plrsd$biomass, lwd=3, col=rtmbcol, lty="dotted")
lines(pyrs, plr$biomass+2*plrsd$biomass, lwd=3, col=rtmbcol, lty="dotted")
lines(pyrs, pcod2020rep$biomass, lwd=3, col=iscamcol)
legend("topright", legend=c("RTMB", "iscam MPD"), lwd=3, col=c(rtmbcol, iscamcol) ,bty="n")

# Rec devs

