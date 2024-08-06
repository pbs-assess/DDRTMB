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

library(tidyverse)
library(reshape2)
library(cowplot)
library(here)
library(gfplot)
source("devs/load-models.R")

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
projoutput_iscam <- projoutput_iscam[31002:62001,]

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
# Get reference points (same for every TAC so just take first)
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
write_csv(refpt_quants_rtmb, here("outputs","Posterior_reference_points_quants.csv"))
write_csv(refpt_quants_iscam, here("outputs","Posterior_reference_points_quants.csv"))

# Biomass
#1. Historical reference points
USRrtmb <- refpt_quants_rtmb %>%
  filter(refpt=="bavg")
LRPrtmb <- refpt_quants_rtmb %>%
  filter(refpt=="bmin")
post_biomass_rp <- mcmcderived$biomass %>%
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
  xlab("Year") + ylab("Posterior biomass (t)")+
  theme_pbs()
#print(post_biomass_rp)
ggsave(here("outputs","figs","RTMB_MCMC_Biomass_HistRefPts.png"), width=8, height=6, units="in")

# 2. MSY reference points
BMSY <- refpt_quants_rtmb %>%
  filter(refpt=="bmsy") %>%
  select(lwr,med,upr)
USRrtmb <- 0.8*BMSY
LRPrtmb <- 0.4*BMSY
post_biomass_rp <- mcmcderived$biomass %>%
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
  xlab("Year") + ylab("Posterior biomass (t)")+
  theme_pbs()
#print(post_biomass_rp)
ggsave(here("outputs","figs","RTMB_MCMC_Biomass_MSYRefPts.png"), width=8, height=6, units="in")

# Fishing mortality
#1. Historical reference points
LRRrtmb <- refpt_quants_rtmb %>%
  filter(refpt=="favg")
post_ft_rp <- mcmcderived$Ft %>%
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
  xlab("Year") + ylab("Posterior biomass (t)")+
  theme_pbs()
#print(post_ft_rp)
ggsave(here("outputs","figs","RTMB_MCMC_Ft_HistRefPts.png"), width=8, height=6, units="in")

# 2. MSY reference points
LRRrtmb <- refpt_quants_rtmb %>%
  filter(refpt=="fmsy")
post_ft_rp <- mcmcderived$Ft %>%
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
  xlab("Year") + ylab("Posterior biomass (t)")+
  theme_pbs()
#print(post_ft_rp)
ggsave(here("outputs","figs","RTMB_MCMC_Ft_MSYRefPts.png"), width=8, height=6, units="in")



