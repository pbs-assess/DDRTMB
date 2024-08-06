# Simple plots to look at RTMB MCMC results
# Robyn Forrest
# Date created:  July 26, 2024
# Last Modified: July 26, 2024
#
# FIRST run the model in run-model-v1.R
# Load documentation and inputs
# Uncomment these 4 lines if running standalone
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
library(coda)
source("devs/load-models.R")
source("devs/mcmc_plots.R")

# Settings
# Colours for data, iscam and RTMB
datcol <- 1
rtmbcol <- "darkblue"

if(!file.exists(here("outputs"))) dir.create(here("outputs"), recursive = TRUE)
if(!file.exists(here("outputs","figs"))) dir.create(here("outputs","figs"), recursive = TRUE)

# 1. Read in MCMC outputs - burnin already removed
mcmcpars <- readRDS(here("outputs","MCMC_parameter_estimates.rda"))
mcmcderived <- readRDS(here("outputs","MCMC_outputs_byvariable.rda"))
mcmcdiagnostics <- readRDS(here("outputs","MCMC_diagnostics.rda"))
projoutput <- readRDS(here("outputs","Projections_output.rda")) # contains reference points

# Put all estimated scalar parameters together
# Normally EIV pars would be here too but they are fixed for pcod
qit <- paste0("q",1)
for(i in 2:dat$nit){
  qit <- c(qit,paste0("q",i))
}
post_pars <- cbind(mcmcpars[,1:3], mcmcderived$q)
colnames(post_pars) <- c("log_ro","h","log_m",qit)

# Pairs and trace plots
# Plots from iscam
png(here("outputs","figs","RTMB_MCMC_Priors_Posts.png"), width=8, height=6, units="in", res=300)
  make.priors.posts.plot(post_pars,
                       ctl=ctl,
                      priors.only = FALSE)
dev.off()

png(here("outputs","figs","RTMB_MCMC_Pairs.png"), width=8, height=6, units="in", res=300)
  make.pairs.plot(post_pars)
dev.off()

png(here("outputs","figs","RTMB_MCMC_Trace.png"), width=8, height=6, units="in", res=300)
  make.traces.plot(post_pars,axis.lab.freq = 200)
dev.off()

png(here("outputs","figs","RTMB_MCMC_Autocor.png"), width=8, height=6, units="in", res=300)
  make.autocor.plot(post_pars,ylim = c(-1,1))
dev.off()

# Biomass, recruits and fishing mortality time series
post_biomass <- mcmcderived$biomass %>%
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
ggsave(here("outputs","figs","RTMB_MCMC_Biomass.png"), width=8, height=6, units="in")

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
ggsave(here("outputs","figs","RTMB_MCMC_Ft.png"), width=8, height=6, units="in")

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
ggsave(here("outputs","figs","RTMB_MCMC_Recruits.png"), width=8, height=6, units="in")

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
ggsave(here("outputs","figs","RTMB_MCMC_Logrecdevs.png"), width=8, height=6, units="in")

# Table of mcmc diagnostics
write_csv(mcmcdiagnostics, here("outputs","RTMB_MCMC_Diagnostic.csv"))

# PLOTS WITH REFERENCE POINTS
# Get reference points quantiles
refpt_quants_rtmb <- projoutput[[1]] %>%
  apply(2,quantile,probs=c(0.025,0.5,0.975))%>%
  t() %>%
  as.data.frame() %>%
  rename(lwr=`2.5%`, med=`50%`, upr=`97.5%`) %>%
  tibble::rownames_to_column(var <- "refpt")

# Table of reference point quantiles
write_csv(refpt_quants_rtmb, here("outputs","Posterior_reference_points_quants.csv"))

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

# Plot densities of reference points and other projection benchmarks and stock status
rtmb_proj  <- projoutput[[1]]%>%
  as.data.frame()

# Bavg
rtmbden <- density(rtmb_proj$bavg)
png(here("outputs","figs","RTMB_RefPt_HistUSR_MCMC.png"), width=8, height=6, units="in", res=300)
  plot(rtmbden, main="USR (Historical average biomass)", col="blue",lwd=2)
  polygon(rtmbden, col=adjustcolor("blue", alpha.f = 0.2))
dev.off()

# FAvg
rtmbden <- density(rtmb_proj$favg)
png(here("outputs","figs","RTMB_RefPt_HistLRR_MCMC.png"), width=8, height=6, units="in", res=300)
  plot(rtmbden, main="LRR (Historical average F)", col="blue",lwd=2)
  polygon(rtmbden, col=adjustcolor("blue", alpha.f = 0.2))
dev.off()

rtmbden <- density(rtmb_proj$bmin)
png(here("outputs","figs","RTMB_RefPt_HistLRP_MCMC.png"), width=8, height=6, units="in", res=300)
  plot(rtmbden, main="LRP (Historical Bmin)", col="blue",lwd=2)
  polygon(rtmbden, col=adjustcolor("blue", alpha.f = 0.2))
dev.off()

# Bmsy
rtmbden <- density(rtmb_proj$bmsy)
png(here("outputs","figs","RTMB_RefPt_Bmsy_MCMC.png"), width=8, height=6, units="in", res=300)
  plot(rtmbden, main="Bmsy", col="blue",lwd=2)
  polygon(rtmbden, col=adjustcolor("blue", alpha.f = 0.2))
dev.off()

# Fmsy
rtmbden <- density(rtmb_proj$fmsy)
png(here("outputs","figs","RTMB_RefPt_Fmsy_MCMC.png"), width=8, height=6, units="in", res=300)
  plot(rtmbden, main="Fmsy", col="blue",lwd=2)
  polygon(rtmbden, col=adjustcolor("blue", alpha.f = 0.2))
dev.off()

# B0
rtmbden <- density(rtmb_proj$bo)
png(here("outputs","figs","RTMB_RefPt_B0_MCMC.png"), width=8, height=6, units="in", res=300)
  plot(rtmbden, main="Fmsy", col="blue",lwd=2)
  polygon(rtmbden, col=adjustcolor("blue", alpha.f = 0.2))
dev.off()

# Look at some stock status indicators
rtmbden <- density(rtmb_proj$B2021)
png(here("outputs","figs","RTMB_RefPt_B2021_MCMC.png"), width=8, height=6, units="in", res=300)
  plot(rtmbden, main="B2021", col="blue",lwd=2)
  polygon(rtmbden, col=adjustcolor("blue", alpha.f = 0.2))
dev.off()

rtmbden <- density(rtmb_proj$B2022)
png(here("outputs","figs","RTMB_RefPt_B2022_MCMC.png"), width=8, height=6, units="in", res=300)
  plot(rtmbden, main="B2022", col="blue",lwd=2)
  polygon(rtmbden, col=adjustcolor("blue", alpha.f = 0.2))
dev.off()

rtmbden <- density(rtmb_proj$B2022B2021)
png(here("outputs","figs","RTMB_RefPt_B2022relB2021_MCMC.png"), width=8, height=6, units="in", res=300)
  plot(rtmbden, main="B2022 relative to B2021", col="blue",lwd=2)
  polygon(rtmbden, col=adjustcolor("blue", alpha.f = 0.2))
dev.off()

rtmbden <- density(rtmb_proj$B2022Bavg)
png(here("outputs","figs","RTMB_RefPt_B2022relHistUSR_MCMC.png"), width=8, height=6, units="in", res=300)
  plot(rtmbden, main="B2022 relative to Historical USR", col="blue",lwd=2)
  polygon(rtmbden, col=adjustcolor("blue", alpha.f = 0.2))
dev.off()

rtmbden <- density(rtmb_proj$B2022Bmin)
png(here("outputs","figs","RTMB_RefPt_B2022relHistLRP_MCMC.png"), width=8, height=6, units="in", res=300)
  plot(rtmbden, main="B2022 relative to Historical LRP", col="blue",lwd=2)
  polygon(rtmbden, col=adjustcolor("blue", alpha.f = 0.2))
dev.off()

