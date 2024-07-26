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

# 1. Read in MCMC outputs - burnin already removed
mcmcpars <- readRDS(here("outputs","MCMCParameterEstimates.rda"))
mcmcderived <- readRDS( here("outputs","MCMCDerivedEstimates.rda"))
mcmcdiagnostics <- readRDS(here("outputs","MCMCDiagnostics.rda"))
# Also the dat and control file for model dimensions and priors settings
dat <- pcod2020dat
ctl <- pcod2020ctl

# Pairs and trace plots
# Trace plots


# Parameter histograms (without and with priors)
# 1. Simple individual histograms
hist(mcmcpars$log_ro)
hist(mcmcpars$h)
hist(exp(mcmcpars$log_m))
par(mfrow=c(2,3))
for(i in 1:dat$nit) hist(mcmcderived$q[,i])

# Biomass, recruits and fishing mortality time series
post_biomass_rtmb <- mcmcderived$biomass %>%
  apply(2,quantile,probs=c(0.025,0.5,0.975))%>%
  t() %>%
  as.data.frame() %>%
  rename(lwr=`2.5%`, med=`50%`, upr=`97.5%`) %>%
  mutate(Year=dat$syr:(dat$nyr+1)) %>%
  ggplot() +
  geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr), fill="blue", alpha = 0.3) +
  geom_line(aes(x=Year,y=med),color="blue")+
  xlab("Year") + ylab("Posterior biomass (t)")+
  theme_pbs()
print(post_biomass)

post_biomass_iscam <- mcmcderived$biomass %>%
  apply(2,quantile,probs=c(0.025,0.5,0.975))%>%
  t() %>%
  as.data.frame() %>%
  rename(lwr=`2.5%`, med=`50%`, upr=`97.5%`) %>%
  mutate(Year=dat$syr:(dat$nyr+1)) %>%
  ggplot() +
  geom_ribbon(aes(x=Year, ymin=lwr, ymax=upr), fill="blue", alpha = 0.3) +
  geom_line(aes(x=Year,y=med),color="blue")+
  xlab("Year") + ylab("Posterior biomass (t)")+
  theme_pbs()
print(post_biomass)
