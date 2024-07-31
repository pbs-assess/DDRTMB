# Simple plots to look at RTMB results
# Robyn Forrest
# Date created:  July 23, 2024
# Last Modified: July 23, 2024
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

# 1. Read in input files and output files from RTMB and iscam
# inputs
dat <- pcod2020dat
ctl <- pcod2020ctl
# RTMB outputs
pl <- readRDS(here("outputs","parameter_estimates.rda"))
plsd <- readRDS(here("outputs","parameter_sds.rda"))
plrad <- readRDS(here("outputs","derived_estimates.rda"))
plradsd <- readRDS(here("outputs","derived_sds.rda"))

yrs <-  dat$syr:dat$nyr
pyrs <- dat$syr:(dat$nyr+1) # includes projection year

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  MPDs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Basic plots - just use base R for now, for quickness
# Catch
png(here("outputs","figs","RTMB_MPD_CatchFit.png"), width=8, height=6, units="in", res=300)
  maxY <- max(dat$catch$value)
  plot(yrs, plrad$ct, col=rtmbcol, lwd=3,type="l", ylim=c(0,maxY), xlab="Years", ylab="Catch (t)")
  points(yrs, dat$catch$value, type="p", pch=19, cex=1.2, col=datcol)
  arrows(x0=yrs, y0=dat$catch$value-0.05*dat$catch$value, x1 = yrs, y1=dat$catch$value+0.05*dat$catch$value, code = 0)
  legend("topright", legend=c("Data", "RTMB MPD"), lwd=3, col=c(datcol, rtmbcol) ,bty="n")
dev.off()

# Indices
png(here("outputs","figs","RTMB_MPD_IndicesFit.png"), width=8, height=6, units="in", res=300)
  par(mfrow=c(2,3))
  for(i in 1:dat$nit){
    obsindex <- dat$indices[[i]]$it
    iwt <- 1/(dat$indices[[i]]$wt) # Equivalent to CV (I think)
    iyrs <- dat$indices[[i]]$iyr
    # right now, rtmb is reporting all the it_hats in one vector.
    # Need to get them out
    if(i==1){
      n1 <- 1
      n2 <- dat$nitnobs[i]
      rtmbindex <- plrad$it_hat_all[n1:n2]
    }else{
      n1 <- n2+1
      n2 <- n1+dat$nitnobs[i]-1
      rtmbindex <- plrad$it_hat_all[n1:n2]
    }

    maxY <- max(c((obsindex+2*iwt*obsindex),rtmbindex))
    plot(iyrs,rtmbindex, col=rtmbcol, lwd=3,type="l", ylim=c(0,maxY), xlab="Years", ylab="Catch (t)")
    points(iyrs, obsindex, type="p", pch=19, cex=1.2, col=datcol)
    arrows(x0=iyrs, y0=obsindex-2*iwt*obsindex, x1 = iyrs, y1=obsindex+2*iwt*obsindex, code = 0)
    legend("topright", legend=c("Data", "RTMB MPD"), lwd=3, col=c(datcol, rtmbcol) ,bty="n")
  }
dev.off()

# Annual mean weight
# This is not robust to there being more than one mean weight series but good for initial testing
obsindex <- dat$meanwtdata[[1]]$it
iwt <- ctl$weight.sig
iyrs <- dat$meanwtdata[[1]]$iyr
rtmbindex <- plrad$annual_mean_weight_all
png(here("outputs","figs","RTMB_MPD_MeanWeightFit.png"), width=8, height=6, units="in", res=300)
  maxY <- max(c((obsindex+2*iwt*obsindex),rtmbindex))
  plot(iyrs,rtmbindex, col=rtmbcol, lwd=3,type="l", ylim=c(0,maxY), xlab="Years", ylab="Annual mean weight (kg)")
  points(iyrs, obsindex, type="p", pch=19, cex=1.2, col=datcol)
  arrows(x0=iyrs, y0=obsindex-2*iwt*obsindex, x1 = iyrs, y1=obsindex+2*iwt*obsindex, code = 0)
  legend("topright", legend=c("Data", "RTMB MPD"), lwd=3, col=c(datcol, rtmbcol) ,bty="n")
dev.off()

# Biomass
png(here("outputs","figs","RTMB_MPD_Biomass.png"), width=8, height=6, units="in", res=300)
  par(mfrow=(c(1,1)))
  maxY <- max(c(plrad$biomass+2*plradsd$biomass))
  plot(pyrs, plrad$biomass, type="l", lwd=3, col=rtmbcol, ylim=c(0,maxY), xlab="Years", ylab="Biomass (t)")
  lines(pyrs, plrad$biomass-2*plradsd$biomass, lwd=1, col=rtmbcol, lty="dotted")
  lines(pyrs, plrad$biomass+2*plradsd$biomass, lwd=1, col=rtmbcol, lty="dotted")
dev.off()

# Log rec devs
png(here("outputs","figs","RTMB_MPD_Logrecdevs.png"), width=8, height=6, units="in", res=300)
  maxY <- max(c(pl$log_rec_devs+2*plsd$log_rec_devs))
  minY <- min(c(pl$log_rec_devs-2*plsd$log_rec_devs))
  plot(yrs, pl$log_rec_devs, pch=19, cex=1.2, col=rtmbcol, ylim=c(minY,maxY), xlab="Years", ylab="log recruit devs")
  abline(h=0, lty=2, lwd=0.5)
  arrows(x0=yrs, y0=(pl$log_rec_devs-2*plsd$log_rec_devs), x1 = yrs, y1=(pl$log_rec_devs+2*plsd$log_rec_devs), code = 0, col=rtmbcol)
dev.off()

# Init log rec devs
png(here("outputs","figs","RTMB_MPD_Initlogrecdevs.png"), width=8, height=6, units="in", res=300)
  inityrs <- (yrs[1]-length(pl$init_log_rec_devs)):(yrs[1]-1)
  maxY <- max(c(pl$init_log_rec_devs+2*plsd$init_log_rec_devs))
  minY <- min(c(pl$init_log_rec_devs-2*plsd$init_log_rec_devs))
  plot(inityrs, pl$init_log_rec_devs, pch=19, cex=1.2, col=rtmbcol, ylim=c(minY,maxY), xlab="Years", ylab="Init log recruit devs")
  abline(h=0, lty=2, lwd=0.5)
  arrows(x0=inityrs, y0=(pl$init_log_rec_devs-2*plsd$init_log_rec_devs), x1 = inityrs, y1=(pl$init_log_rec_devs+2*plsd$init_log_rec_devs), code = 0, col=rtmbcol)
dev.off()

# Derived recruits
png(here("outputs","figs","RTMB_MPD_Recruits.png"), width=8, height=6, units="in", res=300)
  recyrs <- (dat$syr+dat$sage):dat$nyr
  maxY <- max(plrad$rt+2*plradsd$rt)
  plot(recyrs, plrad$rt, pch=19, col=rtmbcol, ylim=c(0,maxY), xlab="Years", ylab="Recruits")
  arrows(x0=recyrs, y0=(plrad$rt-2*plradsd$rt), x1 = recyrs, y1=(plrad$rt+2*plradsd$rt), code = 0, col=rtmbcol)
dev.off()

# Log ft pars
png(here("outputs","figs","RTMB_MPD_Logft.png"), width=8, height=6, units="in", res=300)
  maxY <- max(c(pl$log_ft_pars+2*plsd$log_ft_pars))
  minY <- min(c(pl$log_ft_pars-2*plsd$log_ft_pars))
  plot(yrs, pl$log_ft_pars, pch=19, cex=1.2, col=rtmbcol, ylim=c(minY,maxY), xlab="Years", ylab="Log ft pars")
  abline(h=0, lty=2, lwd=0.5)
  arrows(x0=yrs, y0=(pl$log_ft_pars-2*plsd$log_ft_pars), x1 = yrs, y1=(pl$log_ft_pars+2*plsd$log_ft_pars), code = 0, col=rtmbcol)
dev.off()

# Derived Ft pars
png(here("outputs","figs","RTMB_MPD_Ft.png"), width=8, height=6, units="in", res=300)
  maxY <- max(exp(pl$log_ft_pars)+exp(2*pl$log_ft_pars))
  plot(yrs, exp(pl$log_ft_pars), type="l", lwd=3,col=rtmbcol, ylim=c(-0.1,maxY), xlab="Years", ylab="Derived ft")
  lines(yrs, exp(pl$log_ft_pars)-exp(2*pl$log_ft_pars), lwd=1, col=rtmbcol, lty="dotted")
  lines(yrs, exp(pl$log_ft_pars)+exp(2*pl$log_ft_pars), lwd=1, col=rtmbcol, lty="dotted")
dev.off()

# Parameters
# Leading
png(here("outputs","figs","RTMB_MPD_LeadingParams.png"), width=8, height=6, units="in", res=300)
  Pars <- matrix(nrow=3, ncol=5)
  Pars[1,] <- c(pl$log_ro+2*plsd$log_ro, pl$h +2*plsd$h, pl$log_m+2*plsd$log_m, pl$rho, pl$kappa)
  Pars[2,] <- c(pl$log_ro, pl$h, pl$log_m, pl$rho, pl$kappa)
  Pars[3,] <- c(pl$log_ro-2*plsd$log_ro, pl$h -2*plsd$h, pl$log_m-2*plsd$log_m, pl$rho, pl$kappa)
  minY <- min(Pars[3,])
  maxY <- max(Pars[1,])
  boxplot(Pars, ylim=c(1.5*minY,maxY), names=c("log_ro", "h", "log_m", "rho", "kappa"))
  abline(h=0, lwd=1, lty=2)
dev.off()

png(here("outputs","figs","RTMB_MPD_hm.png"), width=8, height=6, units="in", res=300)
  Pars <- matrix(nrow=3, ncol=2)
  Pars[1,] <- c(pl$h +2*plsd$h, pl$log_m+2*plsd$log_m)
  Pars[2,] <- c(pl$h, pl$log_m)
  Pars[3,] <- c(pl$h -2*plsd$h, pl$log_m-2*plsd$log_m)
  minY <- min(Pars[3,])
  maxY <- max(Pars[1,])
  boxplot(Pars, ylim=c(minY,maxY), names=c("h", "log_m"))
  abline(h=0, lwd=1, lty=2)
dev.off()

# Catchability, q
qnames <- "q1"
for(i in 2:dat$nit){
  qnames <- c(qnames, paste0("q",i))
}
png(here("outputs","figs","RTMB_MPD_q.png"), width=8, height=6, units="in", res=300)
  Pars <- matrix(nrow=3, ncol=dat$nit)
  Pars[1,] <- c(plrad$q+2*plradsd$q)
  Pars[2,] <- c(plrad$q)
  Pars[3,] <- c(plrad$q-2*plradsd$q)
  maxY <- max(Pars[1,])
  boxplot(Pars, ylim=c(0,maxY), names=qnames)
dev.off()
