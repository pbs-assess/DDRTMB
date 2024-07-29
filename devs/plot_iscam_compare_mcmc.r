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
mcmcpars <- readRDS(here("outputs","MCMCParameterEstimates.rda"))
mcmcderived <- readRDS( here("outputs","MCMCDerivedEstimates.rda"))
mcmcdiagnostics <- readRDS(here("outputs","MCMCDiagnostics.rda"))
# Also the dat and control file for model dimensions and priors settings
dat <- pcod2020dat
ctl <- pcod2020ctl

# 2. iscam Read in MCMC outputs - **need to remove burnin**
iscampars <- read_csv(here("data-raw/iscam_mcmc.csv"))
iscambiomass <- read_csv(here("data-raw/iscam_sbt_mcmc.csv"))
iscamrecruits <- read_csv(here("data-raw/iscam_rt_mcmc.csv"))
iscamrecdevs <- read_csv(here("data-raw/iscam_rdev_mcmc.csv"))
iscamfmort <- read_csv(here("data-raw/iscam_ft_mcmc.csv"))

# Parameter histograms (without and with priors)
# 1. Simple individual histograms
hist(mcmcpars$log_ro)
hist(mcmcpars$h)
hist(exp(mcmcpars$log_m))
par(mfrow=c(2,3))
for(i in 1:dat$nit) hist(mcmcderived$q[,i])

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

post_biomass_iscam <- iscambiomass[1001:2000,] %>%
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

# Fishing mortality
post_ft_rtmb <- exp(mcmcpars[,4:68]) %>%
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
post_ft_iscam <- iscamfmort[1001:2000,1:(dat$nyr-dat$syr+1)] %>%
  apply(2,quantile,probs=c(0.025,0.5,0.975)) %>%
  t() %>%
  as.data.frame() %>%
  rename(lwr=`2.5%`, med=`50%`, upr=`97.5%`) %>%
  mutate(Year=dat$syr:dat$nyr) %>%
  ggplot() +
  geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr), fill="red", alpha = 0.3) +
  geom_line(aes(x=Year,y=med),color="red")+
  xlab("Year") + ylab("Posterior biomass (t)")+ggtitle("iscam")+
  theme_pbs()
#print(post_ft_iscam)
cowplot::plot_grid(post_ft_rtmb,post_ft_iscam, ncol=1)
