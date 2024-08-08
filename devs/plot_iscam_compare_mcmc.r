# Simple plots to compare MCMC estimates of results (iscam vs RTMB)
# Robyn Forrest
# Date created:  July 25, 2024
# Last Modified: July 25, 2024

# Write this code once we have developed plots of RTMB MCMC results
# Simple plots to look at RTMB MCMC results
# Robyn Forrest
# Date created:  July 26, 2024
# Last Modified: July 26, 2024
#
# FIRST run the model in run-model-v1.R
# Load documentation and inputs
# # Uncomment these 4 lines if running standalone
# devtools::document()
# devtools::load_all()
# Also the dat and control file for model dimensions and priors settings
# dat <- pcod2020dat
# ctl <- pcod2020ctl
# pfc <- pcod2020pfc

library(tidyverse)
library(reshape2)
library(cowplot)
library(here)
library(gfplot)
source("devs/load-models.R")
source("R/make_decision_table.R")
source("R/make_decision_table_iscam.R")

# Settings
# Colours for data, iscam and RTMB
datcol <- 1
rtmbcol <- "darkblue"

if(!file.exists(here("outputs"))) dir.create(here("outputs"), recursive = TRUE)
if(!file.exists(here("outputs","figs"))) dir.create(here("outputs","figs"), recursive = TRUE)

# 1. RTMB Read in MCMC outputs - burnin already removed
mcmcpars <- readRDS(here("outputs","MCMC_parameter_estimates.rda"))
mcmcderived <- readRDS(here("outputs","MCMC_outputs_byvariable.rda"))
mcmcdiagnostics <- readRDS(here("outputs","MCMC_diagnostics.rda"))
# Projection_model results including ref points
projoutput <- readRDS(here("outputs","Projections_output.rda")) # contains reference points

# 2. iscam Read in MCMC outputs
iscampars <- read_csv(here("data-raw/iscam_mcmc.csv"))
iscambiomass <- read_csv(here("data-raw/iscam_sbt_mcmc.csv"))
iscamrecruits <- read_csv(here("data-raw/iscam_rt_mcmc.csv"))
iscamrecdevs <- read_csv(here("data-raw/iscam_rdev_mcmc.csv"))
iscamfmort <- read_csv(here("data-raw/iscam_ft_mcmc.csv"))
# Projection_model results including ref points
projoutput_iscam <- read_csv(here("data-raw/iscammcmc_proj_Gear1.csv"))

# Remove burn in samples for iscam ouputs
iscampars <- iscampars[1001:2000,]
iscamrecruits <- iscamrecruits[1001:2000,]
iscamrecdevs <- iscamrecdevs[1001:2000,]
iscamfmort <- iscamfmort[1001:2000,1:(dat$nyr-dat$syr+1)] # take  only gear 1, remove survey gears
projoutput_iscam <- projoutput_iscam[30971:62001,] # results in the desired 1000 samples for each TAC

# Parameter histograms (without and with priors)
# 1. Simple individual histograms
#RTMB
png(here("outputs","figs","CompareLeadingParams_MCMC.png"), width=8, height=6, units="in", res=300)
  par(mfrow=c(2,3),mai=rep(0.3,4),oma=c(0.5,0.5,0.1,0.1))
  hist(mcmcpars$log_ro, col="darkblue", main="RTMB log R0", xlab="")
  hist(mcmcpars$h, col="darkblue", main="RTMB h", xlab="")
  hist(mcmcpars$log_m, col="darkblue", main="RTMB log M", xlab="")
  hist(log(iscampars$ro_gr1), col="red", main="iscam log R0", xlab="")
  hist(iscampars$h_gr1, col="red", main="iscam h", xlab="")
  hist(log(iscampars$m_gs1), col="red", main="iscam log M", xlab="")
dev.off()

png(here("outputs","figs","Compareq_MCMC.png"), width=8, height=6, units="in", res=300)
  par(mfrow=c(2,dat$nit),mai=rep(0.3,4),oma=c(0.5,0.5,0.1,0.1))
  for(i in 1:dat$nit) {
    qit <- paste0("q",i)
    hist(mcmcderived$q[,i],col="darkblue", main=qit, xlab="")
  }
  for(i in 1:dat$nit) {
    qit <- paste0("q",i)
    hist(as.matrix(iscampars)[,qit], col="red", main=qit, xlab="")
  }
dev.off()

# Biomass, recruits and fishing mortality time series
# Biomass
post_biomass_rtmb <- mcmcderived$biomass %>%
  apply(2,quantile,probs=c(0.025,0.5,0.975))%>%
  t() %>%
  as.data.frame() %>%
  rename(lwr=`2.5%`, med=`50%`, upr=`97.5%`) %>%
  mutate(Year=dat$syr:(dat$nyr+1)) %>%
  ggplot() +
  geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr), fill="blue", alpha = 0.3) +
  geom_line(aes(x=Year,y=med),color="blue")+
  xlab("Year") + ylab("Posterior biomass (t)")+ggtitle("rtmb")+
  theme_pbs()

post_biomass_iscam <- iscambiomass %>%
  apply(2,quantile,probs=c(0.025,0.5,0.975)) %>%
  t() %>%
  as.data.frame() %>%
  rename(lwr=`2.5%`, med=`50%`, upr=`97.5%`) %>%
  mutate(Year=dat$syr:(dat$nyr+1)) %>%
  ggplot() +
  geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr), fill="red", alpha = 0.3) +
  geom_line(aes(x=Year,y=med),color="red")+
  xlab("Year") + ylab("Posterior biomass (t)")+ggtitle("iscam")+
  theme_pbs()
#print(post_biomass_iscam)
cowplot::plot_grid(post_biomass_rtmb,post_biomass_iscam, ncol=1)
ggsave(here("outputs","figs","CompareBiomass_MCMC.png"), width=8, height=6, units="in")

# Fishing mortality
post_ft_rtmb <- mcmcderived$Ft %>%
  apply(2,quantile,probs=c(0.025,0.5,0.975))%>%
  t() %>%
  as.data.frame() %>%
  rename(lwr=`2.5%`, med=`50%`, upr=`97.5%`) %>%
  mutate(Year=dat$syr:dat$nyr) %>%
  ggplot() +
  geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr), fill="blue", alpha = 0.3) +
  geom_line(aes(x=Year,y=med),color="blue")+
  xlab("Year") + ylab("Posterior fishing mortality")+ggtitle("rtmb")+
  theme_pbs()
#print(post_ft_rtmb)

# Only take gear 1 in (ft for surveys is all 0)
post_ft_iscam <- iscamfmort %>%
  apply(2,quantile,probs=c(0.025,0.5,0.975)) %>%
  t() %>%
  as.data.frame() %>%
  rename(lwr=`2.5%`, med=`50%`, upr=`97.5%`) %>%
  mutate(Year=dat$syr:dat$nyr) %>%
  ggplot() +
  geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr), fill="red", alpha = 0.3) +
  geom_line(aes(x=Year,y=med),color="red")+
  xlab("Year") + ylab("Posterior fishing mortality")+ggtitle("iscam")+
  theme_pbs()
#print(post_ft_iscam)
cowplot::plot_grid(post_ft_rtmb,post_ft_iscam, ncol=1)
ggsave(here("outputs","figs","CompareFt_MCMC.png"), width=8, height=6, units="in")

# Recruits ... small differences, look at devs
post_recruits_rtmb <- mcmcderived$recruits %>%
  apply(2,quantile,probs=c(0.025,0.5,0.975))%>%
  t() %>%
  as.data.frame() %>%
  rename(lwr=`2.5%`, med=`50%`, upr=`97.5%`) %>%
  mutate(Year=(dat$syr+dat$sage):dat$nyr) %>%
  ggplot() +
  geom_pointrange(aes(x=Year,y=med, ymin=lwr, ymax=upr), col="blue", alpha = 0.5) +
  xlab("Year") + ylab("Posterior recruits")+ggtitle("rtmb")+
  theme_pbs()
#print(post_recruits_rtmb)

post_recruits_iscam <- iscamrecruits %>%
  apply(2,quantile,probs=c(0.025,0.5,0.975)) %>%
  t() %>%
  as.data.frame() %>%
  rename(lwr=`2.5%`, med=`50%`, upr=`97.5%`) %>%
  mutate(Year=(dat$syr+dat$sage):dat$nyr) %>%
  ggplot() +
  geom_pointrange(aes(x=Year,y=med, ymin=lwr, ymax=upr), col="red", alpha = 0.5) +
  xlab("Year") + ylab("Posterior recruits")+ggtitle("iscam")+
  theme_pbs()
#print(post_recruits_iscam)
cowplot::plot_grid(post_recruits_rtmb,post_recruits_iscam, ncol=1)
ggsave(here("outputs","figs","CompareRecruits_MCMC.png"), width=8, height=6, units="in")

# Log rec devs - offset. Iscam mcmc_output only writes out 3:nyrs
post_logrecdevs_rtmb <- mcmcpars[,79:141] %>%
  apply(2,quantile,probs=c(0.025,0.5,0.975))%>%
  t() %>%
  as.data.frame() %>%
  rename(lwr=`2.5%`, med=`50%`, upr=`97.5%`) %>%
  mutate(Year=(dat$syr+dat$sage):dat$nyr) %>%
  ggplot() +
  geom_pointrange(aes(x=Year,y=med, ymin=lwr, ymax=upr), col="blue", alpha = 0.5) +
  geom_hline(yintercept=0, lty=2)+
  xlab("Year") + ylab("Posterior log recruit devs")+ggtitle("rtmb")+
  theme_pbs()
#print(post_logrecdevs_rtmb)

post_logrecdevs_iscam <- iscamrecdevs %>%
  apply(2,quantile,probs=c(0.025,0.5,0.975)) %>%
  t() %>%
  as.data.frame() %>%
  rename(lwr=`2.5%`, med=`50%`, upr=`97.5%`) %>%
  mutate(Year=(dat$syr+dat$sage):dat$nyr) %>%
  ggplot() +
  geom_pointrange(aes(x=Year,y=med, ymin=lwr, ymax=upr), col="red", alpha = 0.5) +
  geom_hline(yintercept=0, lty=2)+
  xlab("Year") + ylab("Posterior log recruit devs")+ggtitle("iscam")+
  theme_pbs()
#print(post_logrecdevs_iscam)
cowplot::plot_grid(post_logrecdevs_rtmb,post_logrecdevs_iscam, ncol=1)
ggsave(here("outputs","figs","CompareLogrecdevs_MCMC.png"), width=8, height=6, units="in")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~ Compare projection outputs including ref points ~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get reference points quantiles (same for every TAC so just take first)
refpt_quants_rtmb <- projoutput[[1]] %>%
  apply(2,quantile,probs=c(0.025,0.5,0.975))%>%
  t() %>%
  as.data.frame() %>%
  rename(lwr=`2.5%`, med=`50%`, upr=`97.5%`) %>%
  tibble::rownames_to_column(var <- "refpt")

refpt_quants_iscam <- projoutput_iscam %>%
  filter(TAC==0) %>%
  apply(2,quantile,probs=c(0.025,0.5,0.975))%>%
  t() %>%
  as.data.frame() %>%
  rename(lwr=`2.5%`, med=`50%`, upr=`97.5%`)%>%
  tibble::rownames_to_column(var <- "refpt")

# Table of reference point quantiles
write_csv(refpt_quants_rtmb, here("outputs","Posterior_reference_points_quants_rtmb.csv"))
write_csv(refpt_quants_iscam, here("outputs","Posterior_reference_points_quants_iscam.csv"))

# Biomass
#1. Historical reference points
USRrtmb <- refpt_quants_rtmb %>%
  filter(refpt=="bavg")
LRPrtmb <- refpt_quants_rtmb %>%
  filter(refpt=="bmin")

post_biomass_rp_rtmb <- mcmcderived$biomass %>%
  apply(2,quantile,probs=c(0.025,0.5,0.975))%>%
  t() %>%
  as.data.frame() %>%
  rename(lwr=`2.5%`, med=`50%`, upr=`97.5%`) %>%
  mutate(Year=dat$syr:(dat$nyr+1), USRlwr=USRrtmb$lwr, USRmed=USRrtmb$med, USRupr=USRrtmb$upr,LRPlwr=LRPrtmb$lwr, LRPmed=LRPrtmb$med, LRPupr=LRPrtmb$upr) %>%
  ggplot()+
  geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr), fill="darkgrey", alpha = 0.5) +
  geom_line(aes(x=Year,y=med),color="black")+
  geom_ribbon(aes(x=Year,ymin=USRlwr,ymax=USRupr), fill="green",alpha=0.2)+
  geom_ribbon(aes(x=Year,ymin=LRPlwr,ymax=LRPupr), fill="red",alpha=0.2)+
  geom_line(aes(x=Year,y=USRmed),color="green", lty=2)+
  geom_line(aes(x=Year,y=LRPmed),color="red", lty=2)+
  xlab("Year") + ylab("Posterior biomass (t)")+ggtitle("rtmb")+
  theme_pbs()

USRiscam <- refpt_quants_iscam %>%
  filter(refpt=="BAvgS")
LRPiscam <- refpt_quants_iscam %>%
  filter(refpt=="Bmin")

post_biomass_rp_iscam <- iscambiomass %>%
  apply(2,quantile,probs=c(0.025,0.5,0.975))%>%
  t() %>%
  as.data.frame() %>%
  rename(lwr=`2.5%`, med=`50%`, upr=`97.5%`) %>%
  mutate(Year=dat$syr:(dat$nyr+1), USRlwr=USRiscam$lwr, USRmed=USRiscam$med, USRupr=USRiscam$upr,LRPlwr=LRPiscam$lwr, LRPmed=LRPiscam$med, LRPupr=LRPiscam$upr) %>%
  ggplot()+
  geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr), fill="darkgrey", alpha = 0.5) +
  geom_line(aes(x=Year,y=med),color="black")+
  geom_ribbon(aes(x=Year,ymin=USRlwr,ymax=USRupr), fill="green",alpha=0.2)+
  geom_ribbon(aes(x=Year,ymin=LRPlwr,ymax=LRPupr), fill="red",alpha=0.2)+
  geom_line(aes(x=Year,y=USRmed),color="green", lty=2)+
  geom_line(aes(x=Year,y=LRPmed),color="red", lty=2)+
  xlab("Year") + ylab("Posterior biomass (t)")+ggtitle("iscam")+
  theme_pbs()
cowplot::plot_grid(post_biomass_rp_rtmb,post_biomass_rp_iscam, ncol=1)
ggsave(here("outputs","figs","CompareBiomass_MCMC_HistRefpts.png"), width=8, height=6, units="in")

# Fishing mortality
#1. Historical reference points
LRRrtmb <- refpt_quants_rtmb %>%
  filter(refpt=="favg")

post_ft_rp_rtmb <- mcmcderived$Ft %>%
  apply(2,quantile,probs=c(0.025,0.5,0.975))%>%
  t() %>%
  as.data.frame() %>%
  rename(lwr=`2.5%`, med=`50%`, upr=`97.5%`) %>%
  mutate(Year=dat$syr:(dat$nyr), LRRlwr=LRRrtmb$lwr, LRRmed=LRRrtmb$med, LRRupr=LRRrtmb$upr) %>%
  ggplot()+
  geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr), fill="darkgrey", alpha = 0.3) +
  geom_line(aes(x=Year,y=med),color="black")+
  geom_ribbon(aes(x=Year,ymin=LRRlwr,ymax=LRRupr), fill="darkgray",alpha=0.2)+
  geom_line(aes(x=Year,y=LRRmed),color="black", lty=2)+
  xlab("Year") + ylab("Posterior biomass (t)")+ggtitle("rtmb")+
  theme_pbs()

LRRiscam <- refpt_quants_iscam %>%
  filter(refpt=="FAvgS")

post_ft_rp_iscam <- iscamfmort %>%
  apply(2,quantile,probs=c(0.025,0.5,0.975))%>%
  t() %>%
  as.data.frame() %>%
  rename(lwr=`2.5%`, med=`50%`, upr=`97.5%`) %>%
  mutate(Year=dat$syr:(dat$nyr), LRRlwr=LRRiscam$lwr, LRRmed=LRRiscam$med, LRRupr=LRRiscam$upr) %>%
  ggplot()+
  geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr), fill="darkgrey", alpha = 0.3) +
  geom_line(aes(x=Year,y=med),color="black")+
  geom_ribbon(aes(x=Year,ymin=LRRlwr,ymax=LRRupr), fill="darkgray",alpha=0.2)+
  geom_line(aes(x=Year,y=LRRmed),color="black", lty=2)+
  xlab("Year") + ylab("Posterior biomass (t)")+ggtitle("iscam")+
  theme_pbs()
cowplot::plot_grid(post_ft_rp_rtmb,post_ft_rp_iscam, ncol=1)
ggsave(here("outputs","figs","CompareFt_MCMC_HistRefpts.png"), width=8, height=6, units="in")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot densities of reference points and other projection benchmarks and stock status
rtmb_proj  <- projoutput[[1]]%>%
  as.data.frame()
iscam_proj <- projoutput_iscam %>%
  rename(bavg=BAvgS, favg=FAvgS, bmin=Bmin,
         bo=B0,fmsy=FMSY,bmsy=BMSY,B2022Bavg=B2022BAvgS,
         B202208Bmsy=B20220.8BMSY,B202204Bmsy=B20220.4BMSY,F2021Favg=F2021FAvgS) %>%
  filter(TAC==0) %>%
  as.data.frame()

# Bavg
rtmbden <- density(rtmb_proj$bavg)
iscamden <- density(iscam_proj$bavg)
png(here("outputs","figs","CompareRefPt_HistUSR_MCMC.png"), width=8, height=6, units="in", res=300)
  plot(rtmbden, main="USR (Historical average biomass)", col="blue",lwd=2, ylim=c(0,max(c(rtmbden$y, iscamden$y))), xlim=c(0,max(c(rtmbden$x, iscamden$x))))
  lines(iscamden, col="red",lwd=2)
  polygon(rtmbden, col=adjustcolor("blue", alpha.f = 0.2))
  polygon(iscamden, col=adjustcolor("red", alpha.f = 0.2))
  legend("topright",legend=c("RTMB", "iscam"), lwd=2, col=c("blue","red"), bty="n")
dev.off()

# FAvg
rtmbden <- density(rtmb_proj$favg)
iscamden <- density(iscam_proj$favg)
png(here("outputs","figs","CompareRefPt_HistLRR_MCMC.png"), width=8, height=6, units="in", res=300)
  plot(rtmbden, main="LRR (Historical average F)", col="blue",lwd=2, ylim=c(0,max(c(rtmbden$y, iscamden$y))), xlim=c(0,max(c(rtmbden$x, iscamden$x))))
  lines(iscamden, col="red",lwd=2)
  polygon(rtmbden, col=adjustcolor("blue", alpha.f = 0.2))
  polygon(iscamden, col=adjustcolor("red", alpha.f = 0.2))
  legend("topright",legend=c("RTMB", "iscam"), lwd=2, col=c("blue","red"), bty="n")
dev.off()

rtmbden <- density(rtmb_proj$bmin)
iscamden <- density(iscam_proj$bmin)
png(here("outputs","figs","CompareRefPt_HistLRP_MCMC.png"), width=8, height=6, units="in", res=300)
plot(rtmbden, main="LRP (Historical Bmin)", col="blue",lwd=2, ylim=c(0,max(c(rtmbden$y, iscamden$y))), xlim=c(0,max(c(rtmbden$x, iscamden$x))))
  lines(iscamden, col="red",lwd=2)
  polygon(rtmbden, col=adjustcolor("blue", alpha.f = 0.2))
  polygon(iscamden, col=adjustcolor("red", alpha.f = 0.2))
  legend("topright",legend=c("RTMB", "iscam"), lwd=2, col=c("blue","red"), bty="n")
dev.off()

# Bmsy
rtmbden <- density(rtmb_proj$bmsy)
iscamden <- density(iscam_proj$bmsy)
png(here("outputs","figs","CompareRefPt_Bmsy_MCMC.png"), width=8, height=6, units="in", res=300)
  plot(rtmbden, main="Bmsy", col="blue",lwd=2, ylim=c(0,max(c(rtmbden$y, iscamden$y))), xlim=c(0,max(c(rtmbden$x, iscamden$x))))
  lines(iscamden, col="red",lwd=2)
  polygon(rtmbden, col=adjustcolor("blue", alpha.f = 0.2))
  polygon(iscamden, col=adjustcolor("red", alpha.f = 0.2))
  legend("topright",legend=c("RTMB", "iscam"), lwd=2, col=c("blue","red"), bty="n")
dev.off()

# Fmsy
rtmbden <- density(rtmb_proj$fmsy)
iscamden <- density(iscam_proj$fmsy)
png(here("outputs","figs","CompareRefPt_Fmsy_MCMC.png"), width=8, height=6, units="in", res=300)
  plot(rtmbden, main="Fmsy", col="blue",lwd=2, ylim=c(0,max(c(rtmbden$y, iscamden$y))), xlim=c(0,max(c(rtmbden$x, iscamden$x))))
  lines(iscamden, col="red",lwd=2)
  polygon(rtmbden, col=adjustcolor("blue", alpha.f = 0.2))
  polygon(iscamden, col=adjustcolor("red", alpha.f = 0.2))
  legend("topright",legend=c("RTMB", "iscam"), lwd=2, col=c("blue","red"), bty="n")
dev.off()

# B0
rtmbden <- density(rtmb_proj$bo)
iscamden <- density(iscam_proj$bo)
png(here("outputs","figs","CompareRefPt_B0_MCMC.png"), width=8, height=6, units="in", res=300)
  plot(rtmbden, main="B0", col="blue",lwd=2, ylim=c(0,max(c(rtmbden$y, iscamden$y))), xlim=c(0,max(c(rtmbden$x, iscamden$x))))
  lines(iscamden, col="red",lwd=2)
  polygon(rtmbden, col=adjustcolor("blue", alpha.f = 0.2))
  polygon(iscamden, col=adjustcolor("red", alpha.f = 0.2))
  legend("topright",legend=c("RTMB", "iscam"), lwd=2, col=c("blue","red"), bty="n")
dev.off()

# Look at some stock status indicators
rtmbden <- density(rtmb_proj$B2021)
iscamden <- density(iscam_proj$B2021)
png(here("outputs","figs","CompareRefPt_B2021_MCMC.png"), width=8, height=6, units="in", res=300)
  plot(rtmbden, main="B2021", col="blue",lwd=2, ylim=c(0,max(c(rtmbden$y, iscamden$y))), xlim=c(0,max(c(rtmbden$x, iscamden$x))))
  lines(iscamden, col="red",lwd=2)
  polygon(rtmbden, col=adjustcolor("blue", alpha.f = 0.2))
  polygon(iscamden, col=adjustcolor("red", alpha.f = 0.2))
  legend("topright",legend=c("RTMB", "iscam"), lwd=2, col=c("blue","red"), bty="n")
dev.off()

rtmbden <- density(rtmb_proj$B2022)
iscamden <- density(iscam_proj$B2022)
png(here("outputs","figs","CompareRefPt_B2022_MCMC.png"), width=8, height=6, units="in", res=300)
  plot(rtmbden, main="B2022", col="blue",lwd=2, ylim=c(0,max(c(rtmbden$y, iscamden$y))), xlim=c(0,max(c(rtmbden$x, iscamden$x))))
  lines(iscamden, col="red",lwd=2)
  polygon(rtmbden, col=adjustcolor("blue", alpha.f = 0.2))
  polygon(iscamden, col=adjustcolor("red", alpha.f = 0.2))
  legend("topright",legend=c("RTMB", "iscam"), lwd=2, col=c("blue","red"), bty="n")
dev.off()

rtmbden <- density(rtmb_proj$B2022B2021)
iscamden <- density(iscam_proj$B2022B2021)
png(here("outputs","figs","CompareRefPt_B2022relB2021_MCMC.png"), width=8, height=6, units="in", res=300)
  plot(rtmbden, main="B2022 relative to B2021", col="blue",lwd=2, ylim=c(0,max(c(rtmbden$y, iscamden$y))), xlim=c(0,max(c(rtmbden$x, iscamden$x))))
  lines(iscamden, col="red",lwd=2)
  polygon(rtmbden, col=adjustcolor("blue", alpha.f = 0.2))
  polygon(iscamden, col=adjustcolor("red", alpha.f = 0.2))
  legend("topright",legend=c("RTMB", "iscam"), lwd=2, col=c("blue","red"), bty="n")
dev.off()

rtmbden <- density(rtmb_proj$B2022Bavg)
iscamden <- density(iscam_proj$B2022Bavg)
png(here("outputs","figs","CompareRefPt_B2022relHistUSR_MCMC.png"), width=8, height=6, units="in", res=300)
  plot(rtmbden, main="B2022 relative to Historical USR", col="blue",lwd=2, ylim=c(0,max(c(rtmbden$y, iscamden$y))), xlim=c(0,max(c(rtmbden$x, iscamden$x))))
  lines(iscamden, col="red",lwd=2)
  polygon(rtmbden, col=adjustcolor("blue", alpha.f = 0.2))
  polygon(iscamden, col=adjustcolor("red", alpha.f = 0.2))
  legend("topright",legend=c("RTMB", "iscam"), lwd=2, col=c("blue","red"), bty="n")
dev.off()

rtmbden <- density(rtmb_proj$B2022Bmin)
iscamden <- density(iscam_proj$B2022Bmin)
png(here("outputs","figs","CompareRefPt_B2022relHistLRP_MCMC.png"), width=8, height=6, units="in", res=300)
  plot(rtmbden, main="B2022 relative to Historical LRP", col="blue",lwd=2, ylim=c(0,max(c(rtmbden$y, iscamden$y))), xlim=c(0,max(c(rtmbden$x, iscamden$x))))
  lines(iscamden, col="red",lwd=2)
  polygon(rtmbden, col=adjustcolor("blue", alpha.f = 0.2))
  polygon(iscamden, col=adjustcolor("red", alpha.f = 0.2))
  legend("topright",legend=c("RTMB", "iscam"), lwd=2, col=c("blue","red"), bty="n")
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Compare decision tables and stock status densities for some other TACs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Whole list
rtmb_proj  <- projoutput
# convert dataframe to list
iscam_proj <- projoutput_iscam %>%
  rename(bavg=BAvgS, favg=FAvgS, bmin=Bmin,
         bo=B0,fmsy=FMSY,bmsy=BMSY,B2022Bavg=B2022BAvgS,
         B202208Bmsy=B20220.8BMSY,B202204Bmsy=B20220.4BMSY,
         F2021Favg=F2021FAvgS, F2021Fmsy=F2021FMSY)
iscam_proj <- split(iscam_proj, list(iscam_proj$TAC))

tactest <- c(2,5,10,20,25,30)

# plot B2022/Hist ref points and F2021/Historical ref point for 6 TACs
png(here("outputs","figs","CompareRefPt_F2021relHistLRR_sixtac.png"), width=8, height=6, units="in", res=300)
par(mfrow=c(2,3))
for(i in tactest){
  rtmb_proj_tac <- rtmb_proj[[i]]
  iscam_proj_tac <- iscam_proj[[i]]

  rtmbden <- density(rtmb_proj_tac$F2021Favg)
  iscamden <- density(iscam_proj_tac$F2021Favg)
  plot(rtmbden, main=paste("TAC=",pfc$tac.vec[i]), col="blue",lwd=2,
       ylim=c(0,max(c(rtmbden$y, iscamden$y))),
       xlim=c(0,max(c(rtmbden$x, iscamden$x))),
       xlab="F2021/HistLRR")
  lines(iscamden, col="red",lwd=2)
  polygon(rtmbden, col=adjustcolor("blue", alpha.f = 0.2))
  polygon(iscamden, col=adjustcolor("red", alpha.f = 0.2))
  legend("topright",legend=c("RTMB", "iscam"), lwd=2, col=c("blue","red"), bty="n")
 }
dev.off()

png(here("outputs","figs","CompareRefPt_B2022relHistUSR_sixtac.png"), width=8, height=6, units="in", res=300)
par(mfrow=c(2,3))
for(i in tactest){
  rtmb_proj_tac <- rtmb_proj[[i]]
  iscam_proj_tac <- iscam_proj[[i]]

  rtmbden <- density(rtmb_proj_tac$B2022Bavg)
  iscamden <- density(iscam_proj_tac$B2022Bavg)
  plot(rtmbden, main=paste("TAC=",pfc$tac.vec[i]), col="blue",lwd=2,
       ylim=c(0,max(c(rtmbden$y, iscamden$y))),
       xlim=c(0,max(c(rtmbden$x, iscamden$x))),
       xlab="B2022/HistUSR")
  lines(iscamden, col="red",lwd=2)
  polygon(rtmbden, col=adjustcolor("blue", alpha.f = 0.2))
  polygon(iscamden, col=adjustcolor("red", alpha.f = 0.2))
  legend("topright",legend=c("RTMB", "iscam"), lwd=2, col=c("blue","red"), bty="n")
}
dev.off()

png(here("outputs","figs","CompareRefPt_B2022relHistLRP_sixtac.png"), width=8, height=6, units="in", res=300)
par(mfrow=c(2,3))
for(i in tactest){
  rtmb_proj_tac <- rtmb_proj[[i]]
  iscam_proj_tac <- iscam_proj[[i]]

  rtmbden <- density(rtmb_proj_tac$B2022Bmin)
  iscamden <- density(iscam_proj_tac$B2022Bmin)
  plot(rtmbden, main=paste("TAC=",pfc$tac.vec[i]), col="blue",lwd=2,
       ylim=c(0,max(c(rtmbden$y, iscamden$y))),
       xlim=c(0,max(c(rtmbden$x, iscamden$x))),
       xlab="B2022/HistLRP")
  lines(iscamden, col="red",lwd=2)
  polygon(rtmbden, col=adjustcolor("blue", alpha.f = 0.2))
  polygon(iscamden, col=adjustcolor("red", alpha.f = 0.2))
  legend("topright",legend=c("RTMB", "iscam"), lwd=2, col=c("blue","red"), bty="n")
}
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get decision tables, just for one year projection
decision_table_rtmb  <- make_decision_table(rtmb_proj,npyr=1)
decision_table_iscam <- make_decision_table_iscam(iscam_proj,npyr=1) # has its own function bc no B20xx/B0

# remove column P(B2022<B0) because iscam doesn't have it
decision_table_rtmb <- decision_table_rtmb %>%
  select(-"P(B2022<B0)")

decision_table_rtmb <- rbind(rep("RTMB",ncol(decision_table_rtmb)),decision_table_rtmb)
decision_table_iscam <- rbind(rep("iscam",ncol(decision_table_iscam)),decision_table_iscam)
decision_table_compare <- rbind(decision_table_rtmb,decision_table_iscam)
write_csv(decision_table_compare,here("outputs","Compare_decision_table.csv"))

