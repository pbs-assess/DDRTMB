library(RTMB)
library(here)

## THIS IS THE BABYSAM CODE FROM THE NANAIMO RTMB COURSE, FEB 26 - MARCH 1 2024
## CODE BY ANDERS NIELSEN
## ADDITIONAL COMMENTS BY ROBYN FORREST (RF)

## Comments with double pound sign are RF's comments
# Comments with single pound sign are from Anders' original file

## TODO
## Understand IGAR settings
## It means Irregular Grid Auto Regression. To do with obs error covariance
## Also in the folder sam-obs (IGAR.R), see sam-obs.pdf
## See also US.R and sam-obs.pdf (unstructured covariance)
## Understand how selectivity is dealt with ... looks like directly through F at age

#library(stockassessment)
#fit<-fitfromweb("NEA_sei_21_v5Reca")
load(here("devs","babysam","neasaithe.RData"))

## Data
dat<-list()
dat$logobs <- fit$data$logobs ## note that Anders puts all the obs in one long data object and has indices to get them out
dat$aux <- fit$data$aux ## looks like an index file for identifying location of fleet and age data (?)
dat$minYear <- min(fit$data$years)
dat$minAge <- min(fit$data$minAgePerFleet)
dat$fleetTypes <- fit$data$fleetTypes ## 3 fleets, not sure yet about types codes
dat$sampleTimes <- fit$data$sampleTimes
dat$year <- fit$data$years
dat$age <-  min(fit$data$minAgePerFleet):max(fit$data$maxAgePerFleet) ## 10 age classes (3:12)
dat$M <- fit$data$natMor ## M (constant in age and time)
dat$SW <- fit$data$stockMeanWeight ## Mean weight at age
dat$MO <- fit$data$propMat ## proportion mature each year
dat$PF <- fit$data$propF ## proportion of Z that is F?? - all zero
dat$PM <- fit$data$propM ## proportion of Z that is natM??  - all zero

dat$srmode <- 2 ## SR model (0=Random Walk, 1=Ricker, 2=BH)
dat$fcormode <- 2 ## I think this is switch for autocor in F (0 or not 0, see pars section)
dat$keyF <- fit$conf$keyLogFsta[1,]+1 ## not sure. This has 10 values, i.e., ages? But what do the numbers mean? The max value (6) is used to set the length of the vectorlogsdF
dat$keyQ <- fit$conf$keyLogFpar+1 ## This is a setting for q (see pars section). Seems to be a q for every gear and age
dat$keySd <- fit$conf$keyVarObs+1 ## I think this is a setting for sd of obs? One for every age and year
dat$keySd[dat$keySd<=0] <- NA ## setting zeros in dat$keySd to NA
dat$fleetDim <- apply(dat$keySd,1,function(x)sum(!is.na(x))) ## setting fleet dimension but not sure exactly how
dat$covType <-  as.integer(fit$conf$obsCorStruct)-1 # c(0,1,2) ## not sure
dat$keyIGAR <- fit$conf$keyCorObs +1 ## not sure but see sam-obs.pdf
dat$keyIGAR[fit$conf$keyCorObs==0] <- NA ## not sure but see sam-obs.pdf
dat$keyIGAR[is.na(fit$conf$keyCorObs)] <- -1 ## not sure but see sam-obs.pdf
#dat$keyIGAR[2, 1:4] <- 1

## Parameters - Anders usually sets starting values to 0
## No priors
par <- list()
par$logsdR <- 0 ## sd in recruitment, one value
par$logsdS <- 0 ## Not sure what this is yet, one value
par$logsdF <- numeric(max(dat$keyF)) # sd in F Vector of length 6, which is max(dat$keyF)
par$rickerpar <- if(dat$srmode==1){c(1,1)}else{numeric(0)} ## if ricker, set starting vals to 1
par$transRhoF <- if(dat$fcormode==0){numeric(0)}else{0.1} ## I think this is the transformed rho term for autocor in F (set to 0.1 if non-zero)
par$bhpar <- if(dat$srmode==2){c(1,1)}else{numeric(0)} ## if bh, set starting vals to 1
par$logQ <- numeric(max(dat$keyQ, na.rm=TRUE))-5 ## there seems to be a q for every age and gear
par$logsdO <- numeric(max(dat$keySd, na.rm=TRUE)) ## Obs errors for each gear?
par$logIGARdist <- if(sum(dat$covType==1)==0){numeric(0)}else{numeric(max(dat$keyIGAR,na.rm=TRUE))} ## not sure but see sam-obs.pdf
par$parUS <- unlist(sapply(1:length(dat$covType), function(f)if(dat$covType[f]==2)unstructured(dat$fleetDim[f])$parms())) ## This is for unstructured covariance (see sam-obs.pdf)
par$missing <- numeric(sum(is.na(dat$logobs))) ## missing is identified in makeADFun as a random effect so this is probably accounting for missing data (none here)
par$logN <- matrix(0, nrow=length(dat$year), ncol=length(dat$age)) ## matrix for Na,t initialize at zero
par$logF <- matrix(0, nrow=length(dat$year), ncol=max(dat$keyF)) ## matrix for Fa,t initialize at zero

## spawning biomass - what are PF and PM for?
## The decay in N is already taken care of in the calc of numbers?
## When would you set non-zero values of PF and PM
ssbFUN <- function(N, F, M, SW, MO, PF, PM){
  rowSums(N*exp(-PF*F-PM*M)*MO*SW)
}

f<-function(par){
  getAll(par,dat) ## puts all the parameters and data into the global environment (like attach)
  logobs <- OBS(logobs) ## the OBS function tells RTMB what is an observation ... not sure if you always need it but do need it for simulating from the data
  logobs[is.na(logobs)] <- missing # I think this sets missing values to zero but the value of missing is zero because there were no missing obs so not sure how this works if the missing calc above is > 0
  nobs <- length(logobs)
  nrow <- nrow(M) ## years in rows
  ncol <- ncol(M) ## ages in columns
  ## exponentiate all the estimated sds, estimated in log space
  sdR <- exp(logsdR) ## sd in recruitment, one value (process error)
  sdS <- exp(logsdS) ## sd in numbers at age, one value (process error)
  sdF <- exp(logsdF) ## sd in F Vector of length 6. M & F for 3 fleets?
  sdO <- exp(logsdO) ## sd in observations for each fleet
  logFF <- logF[,keyF] # expand F ## still not sure why logF has 6 columns
  F <- exp(logFF) ## probably should call F something other than F, maybe Fmort

  ssb <- ssbFUN(exp(logN),exp(logFF),M,SW,MO,PF,PM)
  jnll <- 0 ## initialize joint neg log likelihood

  ## loop over years
  ## note Anders calculates jnll at each time and age step
  ## get recruits from ssb and the stock-recruit curve
  for(y in 2:nrow){
    thisSSB <- ifelse((y-minAge-1)>(-.5),ssb[y-minAge],ssb[1])
    if(srmode==0){ #RW ## not sure where to find the parameters of the RW. In the nll?
      pred <- logN[y-1,1]
    }
    if(srmode==1){ #Ricker
      pred <- rickerpar[1]+log(thisSSB)-exp(rickerpar[2])*thisSSB
    }
    if(srmode==2){ #BH
      pred <- bhpar[1]+log(thisSSB)-log(1.0+exp(bhpar[2])*thisSSB)
    }
    if(!(srmode%in%c(0,1,2))){
      stop(paste("srmode", srmode,"not implemented yet"))
    }
    jnll <- jnll - dnorm(logN[y,1],pred,sdR,log=TRUE) ## Add recruits to nll
  }
  for(y in 2:nrow){
    for(a in 2:ncol){
      pred <- logN[y-1,a-1]-F[y-1,a-1]-M[y-1,a-1] ## predict log numbers at age
      if(a==ncol){
        pred <- log(exp(pred)+exp(logN[y-1,a]-F[y-1,a]-M[y-1,a])) ## plus group
      }
      jnll <- jnll - dnorm(logN[y,a],pred,sdS,log=TRUE) ##!! add process error in N at age to nll - because this is a state space model - do we want process errors on every N at age??
    }
  }

  ## Need to better understand what is going on with covariance here
  ## look at sam-obs.pdf
  SigmaF <- matrix(0, ncol(logF),ncol(logF))
  if(fcormode==0){
    diag(SigmaF) <- sdF*sdF
  }
  if(fcormode==1){
    diag(SigmaF) <- sdF*sdF
    rhoF <- 2*plogis(transRhoF[1])-1 ## I think this is autocorrelation in sigmaF among fleets
    for(i in 2:ncol(logF)){
      for(j in 1:(i-1)){
        SigmaF[i,j] <- rhoF*sdF[i]*sdF[j]
        SigmaF[j,i] <- SigmaF[i,j]
      }
    }
  }
  if(fcormode==2){
    diag(SigmaF) <- sdF*sdF
    rhoF <- 2*plogis(transRhoF[1])-1
    for(i in 2:ncol(logF)){
      for(j in 1:(i-1)){
        SigmaF[i,j] <- sdF[i]*sdF[j]*(rhoF^(i-j))
        SigmaF[j,i] <- SigmaF[i,j]
      }
    }
  }
  for(y in 2:nrow){
    jnll <- jnll - dmvnorm(logF[y,],logF[y-1,],SigmaF,log=TRUE) #modelling F as mvnorm process
  }

  logPred <- numeric(nobs)
  for(i in 1:nobs){
    y <- aux[i,1]-minYear+1 ## year
    f <- aux[i,2] ## fleet
    a <- aux[i,3]-minAge+1 ## age
    Z <- F[y,a]+M[y,a] ## Total mortality Z

    if(fleetTypes[f]==0){
      logPred[i] <- logN[y,a]-log(Z)+log(1-exp(-Z))+log(F[y,a]) ## commercial catch?
    }
    if(fleetTypes[f]==2){
      logPred[i] <- logQ[keyQ[f,a]]+logN[y,a]-Z*sampleTimes[f] ## surveys?
    }
    if(!(fleetTypes[f]%in%c(0,2))){
      stop("This fleet type is has not been implemented yet")
    }
  }
  Slist<-list()
  for(f in unique(aux[,2])){ # each fleet
    if(covType[f]==0){# independent
      S <- diag(sdO[na.omit(keySd[f,])]^2)
    }
    if(covType[f]==1){# IGAR ## irregular grid AR1 - covariance across fleets - see sam-obs.pdf
      S <- diag(sdO[na.omit(keySd[f,])]^2)
      dist <- cumsum(c(0,exp(logIGARdist[na.omit(keyIGAR[f,])])))
      for(i in 2:ncol(S)){
        for(j in 1:(i-1)){
          S[i,j] <- sqrt(S[i,i])*sqrt(S[j,j])*(0.5^(dist[i]-dist[j]))
          S[j,i] <- S[i,j]
        }
      }
    }
    if(covType[f]==2){# US ## unstructured covariance - see sam-obs.pdf
      idx <- which(f==unlist(sapply(1:length(dat$covType), function(f)if(dat$covType[f]==2)unstructured(dat$fleetDim[f])$parms()+f)))
      D <- diag(sdO[na.omit(keySd[f,])])
      R <- unstructured(nrow(D))$corr(parUS[idx])
      S <- D%*%R%*%D
    }
    if(!covType[f]%in%c(0,1,2)){#
      stop("Covariance type not implemented")
    }
    Slist[[length(Slist)+1]]<-S
    for(y in unique(aux[,1])){ # year within fleet
      idx <- which((aux[,2]==f) & (aux[,1]==y))
      if(length(idx)!=0){
        jnll <- jnll - dmvnorm(logobs[idx],logPred[idx],S,log=TRUE) ## add obs err to nll
      }
    }
  }
  ## reporting. See Likelihood.pdf in MaxLikelihood folder (slide Getting results out)
  ## Note the more things you report, the longer it takes to run sdreport
  ## From RTMB help:
  ## ADREPORT(): Can be used inside the objective function to report quantities for which uncertainties will be calculated by sdreport.
  ## REPORT(): Can be used inside the objective function to report quantities via the model object using obj$report()
  REPORT(Slist)
  REPORT(logPred)
  logssb<-log(ssb)
  ADREPORT(logssb)
  jnll ## return obj function value
}

## MakeADFun builds the graph, basically "compiles" the model with random effects identified
## from TMB help: map = List defining how to optionally collect and fix parameters
## Means you can fix some instances of a vector/matrix of parameters and estimate ones with the same factor id to be the same
## Might be good for q or selectivity for example when you want all the same value for a given age
obj <- MakeADFun(f, par, random=c("logN", "logF", "missing"), map=list(logsdF=as.factor(rep(0,length(par$logsdF)))), silent=FALSE)
# The optimization step - gets passed the parameters, likelihood function and gradients Makeby ADFun
opt <- nlminb(obj$par, obj$fn, obj$gr, control=list(eval.max=1000, iter.max=1000))
opt$objective

## looks like there is a stock assessment package for making graphs
stockassessment::ssbplot(fit)
sdr<-sdreport(obj)
plr<-as.list(sdr,report=TRUE, "Est") # param estimates
plrsd<-as.list(sdr,report=TRUE, "Std") # param sds
lines(dat$year, exp(plr$logssb), lwd=3, col="darkred")
lines(dat$year, exp(plr$logssb-2*plrsd$logssb), lwd=3, col="darkred", lty="dotted")
lines(dat$year, exp(plr$logssb+2*plrsd$logssb), lwd=3, col="darkred", lty="dotted")

## Look at parameters
plr
## Look at Numbers at age
