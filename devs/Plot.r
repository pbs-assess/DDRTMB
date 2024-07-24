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
# Colours for data, iscam and RTMB
datcol <- 1
iscamcol <- "red"
rtmbcol <- "darkblue"

# 1. Read in input files and output files from RTMB and iscam
# inputs
dat <- pcod2020dat
ctl <- pcod2020ctl
# RTMB outputs
pl <- readRDS(here("outputs","ParameterEstimates.rda"))
plsd <- readRDS(here("outputs","ParameterSDs.rda"))
plr <- readRDS(here("outputs","DerivedEstimates.rda"))
plrsd <- readRDS(here("outputs","DerivedSDs.rda"))
# iscam outputs
pcod2020rep <- read.report.file("data-raw/iscam.rep")
pcod2020par <- read.par.file("data-raw/iscam.par")

yrs <-  dat$syr:dat$nyr
pyrs <- dat$syr:(dat$nyr+1) # includes projection year

# Basic plots
# Catch
maxY <- max(dat$catch$value)
plot(yrs, pcod2020rep$ct, col=iscamcol, lwd=3,type="l", ylim=c(0,maxY), xlab="Years", ylab="Catch (t)")
lines(yrs, plr$ct, col=rtmbcol, lwd=3, lty=2)
points(yrs, dat$catch$value, type="p", pch=19, cex=1.2, col=datcol)
arrows(x0=yrs, y0=dat$catch$value-0.05*dat$catch$value, x1 = yrs, y1=dat$catch$value+0.05*dat$catch$value, code = 0)
legend("topright", legend=c("Data", "RTMB", "iscam MPD"), lwd=3, col=c(datcol, rtmbcol, iscamcol) ,bty="n")

# Indices
par(mfrow=c(2,3))
for(i in 1:dat$nit){
  obsindex <- dat$indices[[i]]$it
  iwt <- 1/(dat$indices[[i]]$wt) # Equivalent to CV (I think)
  iyrs <- dat$indices[[i]]$iyr
  iscamindex <- pcod2020rep$it_hat[i,][which(!is.na(pcod2020rep$it_hat[i,]))]
  # right now, rtmb is reporting all the it_hats in one vector.
  # Need to get them out
  if(i==1){
    n1 <- 1
    n2 <- dat$nitnobs[i]
    rtmbindex <- plr$it_hat_all[n1:n2]
  }else{
    n1 <- n2+1
    n2 <- n1+dat$nitnobs[i]-1
    rtmbindex <- plr$it_hat_all[n1:n2]
  }

  maxY <- max(c((obsindex+2*iwt*obsindex),iscamindex,rtmbindex))
  plot(iyrs,iscamindex, col=iscamcol, lwd=3,type="l", ylim=c(0,maxY), xlab="Years", ylab="Catch (t)")
  lines(iyrs,rtmbindex, col=rtmbcol, lwd=3, lty=1)
  points(iyrs, obsindex, type="p", pch=19, cex=1.2, col=datcol)
  arrows(x0=iyrs, y0=obsindex-2*iwt*obsindex, x1 = iyrs, y1=obsindex+2*iwt*obsindex, code = 0)
  legend("topright", legend=c("Data", "RTMB", "iscam MPD"), lwd=3, col=c(datcol, rtmbcol, iscamcol) ,bty="n")
}

# Biomass
par(mfrow=(c(1,1)))
maxY <- max(c(pcod2020rep$biomass,(plr$biomass+2*plrsd$biomass)))
plot(pyrs, plr$biomass, type="l", lwd=3, col=rtmbcol, ylim=c(0,maxY), xlab="Years", ylab="Biomass (t)")
lines(pyrs, plr$biomass-2*plrsd$biomass, lwd=3, col=rtmbcol, lty="dotted")
lines(pyrs, plr$biomass+2*plrsd$biomass, lwd=3, col=rtmbcol, lty="dotted")
lines(pyrs, pcod2020rep$biomass, lwd=3, col=iscamcol)
legend("topright", legend=c("RTMB", "iscam MPD"), lwd=3, col=c(rtmbcol, iscamcol) ,bty="n")

# Log rec devs
maxY <- max(c(pcod2020par$log_rec_devs,(pl$log_rec_devs+2*plsd$log_rec_devs)))
minY <- min(c(pcod2020par$log_rec_devs,(pl$log_rec_devs-2*plsd$log_rec_devs)))
plot(yrs, pl$log_rec_devs, pch=19, cex=1.2, col=rtmbcol, ylim=c(minY,maxY), xlab="Years", ylab="log recruit devs")
abline(h=0, lty=2, lwd=0.5)
arrows(x0=yrs, y0=(pl$log_rec_devs-2*plsd$log_rec_devs), x1 = yrs, y1=(pl$log_rec_devs+2*plsd$log_rec_devs), code = 0, col=rtmbcol)
points(yrs, pcod2020par$log_rec_devs, pch=1, col=iscamcol)
legend("topright", legend=c("RTMB", "iscam MPD"), pch=c(19,1), col=c(rtmbcol, iscamcol) ,bty="n")

# Init log rec devs
inityrs <- (yrs[1]-length(pcod2020par$init_log_rec_devs)):(yrs[1]-1)
maxY <- max(c(pcod2020par$init_log_rec_devs,(pl$init_log_rec_devs+2*plsd$init_log_rec_devs)))
minY <- min(c(pcod2020par$init_log_rec_devs,(pl$init_log_rec_devs-2*plsd$init_log_rec_devs)))
plot(inityrs, pl$init_log_rec_devs, pch=19, cex=1.2, col=rtmbcol, ylim=c(minY,maxY), xlab="Years", ylab="Init log recruit devs")
abline(h=0, lty=2, lwd=0.5)
arrows(x0=inityrs, y0=(pl$init_log_rec_devs-2*plsd$init_log_rec_devs), x1 = inityrs, y1=(pl$init_log_rec_devs+2*plsd$init_log_rec_devs), code = 0, col=rtmbcol)
points(inityrs, pcod2020par$init_log_rec_devs, pch=1, col=iscamcol)
legend("topright", legend=c("RTMB", "iscam MPD"), pch=c(19,1), col=c(rtmbcol, iscamcol) ,bty="n")

# Derived recruits
recyrs <- (dat$syr+dat$sage):dat$nyr
maxY <- max(c(pcod2020rep$rt,(plr$rt+2*plrsd$rt)))
plot(recyrs, plr$rt, pch=19, col=rtmbcol, ylim=c(0,maxY), xlab="Years", ylab="Recruits")
arrows(x0=recyrs, y0=(plr$rt-2*plrsd$rt), x1 = recyrs, y1=(plr$rt+2*plrsd$rt), code = 0, col=rtmbcol)
points((recyrs+0.5), pcod2020rep$rt, pch=19, col=iscamcol)
legend("topright", legend=c("RTMB", "iscam MPD"), pch=19, col=c(rtmbcol, iscamcol) ,bty="n")

# Log ft pars
maxY <- max(c(pcod2020par$log_ft_pars,(pl$log_ft_pars+2*plsd$log_ft_pars)))
minY <- min(c(pcod2020par$log_ft_pars,(pl$log_ft_pars-2*plsd$log_ft_pars)))
plot(yrs, pl$log_ft_pars, pch=19, cex=1.2, col=rtmbcol, ylim=c(minY,maxY), xlab="Years", ylab="Log ft pars")
abline(h=0, lty=2, lwd=0.5)
arrows(x0=yrs, y0=(pl$log_ft_pars-2*plsd$log_ft_pars), x1 = yrs, y1=(pl$log_ft_pars+2*plsd$log_ft_pars), code = 0, col=rtmbcol)
points(yrs, pcod2020par$log_ft_pars, pch=1, col=iscamcol)
legend("topright", legend=c("RTMB", "iscam MPD"), pch=c(19,1), col=c(rtmbcol, iscamcol) ,bty="n")

# Derived Ft pars
maxY <- max(c(exp(pcod2020par$log_ft_pars),(exp(pl$log_ft_pars))))
plot(yrs, exp(pl$log_ft_pars), type="l", lwd=3,col=rtmbcol, ylim=c(0,maxY), xlab="Years", ylab="Derived ft")
lines(yrs, exp(pcod2020par$log_ft_pars), lwd=3, col=iscamcol)
legend("topright", legend=c("RTMB", "iscam MPD"), lwd=3, col=c(rtmbcol, iscamcol) ,bty="n")

# Parameters
Pars <- matrix(nrow=2, ncol=5)
Pars[1,] <- c(pl$log_ro, pl$h, pl$log_m, pl$rho, pl$kappa)
Pars[2,] <- c(pcod2020par$theta1,pcod2020par$theta2,pcod2020par$theta3,pcod2020par$theta6,pcod2020par$theta7)
colnames(Pars) <- c("log_ro", "h", "log_m", "rho", "kappa")
barplot(Pars, beside=T, col=c(rtmbcol, iscamcol), ylim=c(min(Pars),max(Pars)))
legend("topright", legend=c("RTMB", "iscam MPD"), pch=15, col=c(rtmbcol, iscamcol) ,bty="n")
