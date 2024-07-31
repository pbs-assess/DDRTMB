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
# Uncomment these 2 lines if running standalone
devtools::document()
devtools::load_all()

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
mcmcderived <- readRDS( here("outputs","MCMC_derived_estimates.rda"))
mcmcdiagnostics <- readRDS(here("outputs","MCMC_diagnostics.rda"))

# Also the dat and control file for model dimensions and priors settings
dat <- pcod2020dat
ctl <- pcod2020ctl

# 2. iscam Read in MCMC outputs
iscampars <- read_csv(here("data-raw/iscam_mcmc.csv"))
iscambiomass <- read_csv(here("data-raw/iscam_sbt_mcmc.csv"))
iscamrecruits <- read_csv(here("data-raw/iscam_rt_mcmc.csv"))
iscamrecdevs <- read_csv(here("data-raw/iscam_rdev_mcmc.csv"))
iscamfmort <- read_csv(here("data-raw/iscam_ft_mcmc.csv"))

# Remove burn in samples
iscampars <- iscampars[1001:2000,]
iscamrecruits <- iscamrecruits[1001:2000,]
iscamrecdevs <- iscamrecdevs[1001:2000,]
iscamfmort <- iscamfmort[1001:2000,1:(dat$nyr-dat$syr+1)] # take  only gear 1, remove survey gears

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


