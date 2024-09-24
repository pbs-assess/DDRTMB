#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MODEL
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# model (based on iscam delay difference model (see devs\iscam.tpl))

# Author: Robyn Forrest
# Date created:  July 26, 2024
# Last Modified: August 1, 2024

# This model is called [in run-model-v1]*
#   * This model will be called by a yet-to-be-written function

# This first draft attempts to reproduce the iscam delay difference model from
#   the 2020 Pacific Cod stock assessment for Haida Gwaii/Queen Charlotte
#   Sound stock in Pacific Canada (Area 5ABCD) [published in 2021]

# **For this first version, this is a ONE group, ONE area, ONE sex model**
# M is fixed, not time-varying - adapt later for tvm

# arguments:
# pars= a list of starting parameters

# Returns:
# objective function value

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Model
# Eventually define some components as functions and move to separate R scripts
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
model <- function(par){

  `[<-` <- RTMB::ADoverload("[<-") # Need this to avoid problem of some variables being reassigned from ADvariable

  # Should we put ctl in here too?
  getAll(par,dat, pfc) # RTMB function. Puts arguments into global space

  # Pseudocode from iscam
  # 1. initParameters()
  # 2. calcTotalMortality_deldiff()
  # 3. calcNumbersBiomass_deldiff()
  # 4. calcFisheryObservations_deldiff()
  # 5. calcSurveyObservations_deldiff()
  # 6. calcStockRecruitment_deldiff()
  # 7. calcAnnualMeanWeight_deldiff() //RF added this for P cod - only gets added to objective function if cntrl(15)==1
  # 8. calcObjectiveFunction()

  #|---------------------------------------------------------------------|
  # probably don't need to use par$ anywhere in model bc of the getAll function above
  theta <- c(par$log_ro,
             par$h,
             par$log_m,
             par$rho,
             par$kappa)

  # Initialize objective function components at 0
  nlvec_dd_ct <- numeric(1) # initialize joint neg log likelihood for catch data (nlvec_dd[[1]] in iscam)
  nlvec_dd_it <- numeric(nit) # initialize joint neg log likelihood for survey data (nlvec_dd[[2]] in iscam)
  nlvec_dd_rt <- numeric(1) # initialize joint neg log likelihood for recruitment (nlvec_dd[[3]] in iscam)
  nlvec_dd_wt <- numeric(nmeanwt) # initialize joint neg log likelihood for mean weight (nlvec_dd[[4]] in iscam)
  priors   <- numeric(length(theta)) # initialize priors component of neg log likelihood (priors in iscam)
  qvec     <- numeric(nit) # initialize q priors component of neg log likelihood (qvec in iscam)
  pvec     <- numeric(5) # initialize penalties component of neg log likelihood (pvec in iscam)
  objfun   <- numeric(1) # objective function to be minimized

  # 1. initParameters
  ro        <- exp(theta[1])
  steepness <- theta[2]
  m         <- exp(theta[3])
  rho       <- theta[4]
  kap       <- theta[5] #kappa in par file. Call it kap here bc kappa is an R function

  # don't really need to do this bc getAll function puts the pars into global space
  log_ft_pars <- par$log_ft_pars
  init_log_rec_devs <- par$init_log_rec_devs
  log_rec_devs <- par$log_rec_devs

  # Fixed parameters
  # Variances fixed for P cod
  varphi    <- sqrt(1.0/kap)
  sig       <- sqrt(rho)*varphi # 0.2 for P cod
  tau       <- sqrt(1.0-rho)*varphi # 0.8 for P. cod

  # A decision was made in 2018 to fix these two parameters to be the same
  #   as log_ro (P. Starr). May revisit this later
  log_avgrec  <- log(ro)
  log_recinit <- log(ro)
  sig_c <- ctl$misc[4] # sd in catch likelihood
  sig_w <- ctl$weight.sig # sd in mean weight likelihood (called weight_sig in iscam)
  mean_f <- ctl$misc[7] # mean f constraint for penalty function
  sig_f <- ctl$misc[9] # sd constraint for penalty function

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~ TESTING VALUES FROM data-raw/iscam.rep and data-raw/iscam.par (2020 PCOD RESULTS) ~
  if(test==T){
    # DELETE THIS ONCE MODEL EQUATIONS ARE TESTED
     ro        <- 3376.55
     steepness <- 0.809901
     m         <- 0.306306
     rho       <- 0.058824
     kap       <- 1.470588
     log_avgrec <- log(ro)
     log_recinit <- log(ro)

     # reverse engineer the logs for ro and m, used in the priors calcs
     theta <- c(log(ro),
                steepness,
                log(m),
                rho,
                kap)

    # from iscam.par file
    log_ft_pars <- c(-2.29361, -1.81516, -1.52265, -1.75990, -1.99738, -2.38600, -1.96696, -2.56048, -2.04428, -1.65530, -1.44354, -1.56225, -1.68573, -1.89923, -2.53602, -2.00452, -1.99661, -2.16381, -1.87796, -1.72526, -1.78542, -1.97597, -2.04832, -1.38198, -1.55781, -1.92887, -1.98627, -2.11772, -2.25886, -2.71592, -1.82257, -1.02055, -1.41985, -1.95312, -1.88346, -1.08038, -1.18403, -1.29623, -2.10493, -2.29368, -2.19147, -2.11395, -2.28967, -2.53718, -2.71805, -3.39986, -3.27054, -2.89189, -2.72193, -2.52443, -2.59201, -3.34633, -3.53683, -2.94841, -2.34753, -2.56502, -2.88457, -2.87729, -2.74026, -2.64641, -3.14998, -3.54423, -3.96561, -3.51804, -3.65543)
    init_log_rec_devs <- c(-0.297834, -0.195054, -0.126748, -0.0955628, -0.0934589, -0.108359, -0.130301, 1.04734)
    log_rec_devs <- c(1.05722, 1.10583, -0.139089, -0.165389, -0.298059, -0.336892, -0.173463, 2.84111, 0.284625, 0.163418, -0.0760200, -0.352092, -0.626335, -0.538303, -0.320139, -0.0816409, 2.69634, 0.0765257, 0.524992, 0.510128, 0.356662, 0.953328, 0.574398, 0.840802, 0.173325, 0.402038, 0.278233, -0.103700, 0.166054, 0.213154, 1.49743, 2.13800, -0.221516, -0.0713425, 0.874159, 1.27436, -0.245994, -0.775609, -0.898877, -0.701367, -0.142345, -0.829222, -0.954500, -1.11217, -1.11537, 0.209017, 0.409310, -0.409217, -0.845547, -1.24699, -1.39305, -1.25216, -0.294358, 0.668812, 0.131646, -0.489765, -0.691204, -0.667682, -0.629868, -0.792061, -0.796493, -0.646523, 0.347852, -0.110935, -0.232896)

    # Variances fixed for P cod
    varphi    <- sqrt(1.0/kap)
    sig       <- sqrt(rho)*varphi # 0.2 for P cod
    tau       <- sqrt(1.0-rho)*varphi # 0.8 for P. cod
   } # end if test
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Get the Goodyear Compensation Ratio

  # This is the ratio of juvenile survival at low stock size (alpha) to unfished juvenile survival
  # Calculation differs between BH and Ricker S-R functions (see Forrest et al. 2010 and refs therein https://doi.org/10.1139/F10-077)
  # **IMPORTANT: Steve called this kappa in iscam**
  # **To avoid confusion with the leading variance parameter kappa**
  # **RF has renamed this parameter CR**
  if(ctl$misc[2]==1){
    # Beverton-Holt formulation
    CR <- (4*steepness)/(1-steepness)
  }
  if(ctl$misc[2]==2){
    # Ricker formulation
    CR <- (5*steepness)^1.25
  }
  if(!ctl$misc[2] %in% 1:2){
    stop("S-R relationship not specified! ctl$misc[2] must be set to 1 (BH) or 2 (Ricker).")
  }

  # End initParameters
  #|---------------------------------------------------------------------|

  #|---------------------------------------------------------------------|
  # 2. calcTotalMortality_deldiff();
  # Purpose: This function calculates fishing mortality, total mortality and annual
  #          survival rates S=exp(-Z) for each year.
  #          Z also is updated with time-varying natural mortality rates if
  #          specified by user.

  # Set up F objects
  # In the one fleet case for the delay-diff model, these are the same.
  ft <- matrix(0, nrow=ngear, ncol=nyrs) # a matrix of fishing mortality for each gear and year
  Ft <- numeric(nyrs) # a vector of total F for each year (in the ASM this would account for selectivity but for delay-diff, selectivity is knife-edged)

  # Work down catch obs matrix to associate Fs with catch for each gear
  # The catch matrix in the data file has 4 columns:
  #   year, gear, type (1=catch in weight (tonnes); 2=catch in numbers), catch
  # nctobs is from the dat list and indicates the total number of catch obs
  for(ii in 1:nctobs){
    # Set up counters
    iyear <- catch[ii,1] # actual year
    i <- as.integer(year_lookup[which(year_lookup[,1]==iyear),2]) # year index
    k <- catch[ii,2] # gear

    ftmp    <- exp(log_ft_pars[ii]) # log_ft_pars has length nctobs
    # ft is a matrix with ngear rows and nyr columns
    ft[k,i] <- ftmp # fishing mortality for gear k in year i (in ASM this is modified by selectivity)
    Ft[i]   <- Ft[i] + ftmp # Total fishing mortality in year i
  } # end for ii

  Mt <- rep(m,nyrs) # natural mortality
  Zt <- Ft+Mt       # total mortality
  surv <- exp(-Zt) # fished survival rate

  # End calcTotalMortality_deldiff
  # #|---------------------------------------------------------------------|

  # #|---------------------------------------------------------------------|
  # 3. calcNumbersBiomass_deldiff();
  # Purpose: This function calculates  total biomass and total numbers
  # according to the delay difference model equations from Hilborn and Walters.
  # Quantities are calculated for each year.

  # DD initialization
  # Equilibrium mean weight - obtained from weq = Beq/Neq and solving for weq
  # i.e., weq = [surv(alpha.Neq + rho.Beq + wk.R] / [surv.Neq + R]
  #  with substitutions Neq = Beq/weq and R = Neq(1 - surv)
  # From SJDM, also used by Sinclair in 2005 p cod assessment

  # Initialize numbers and biomass
  numbers <- numeric(length=nyrs+1)
  biomass <- numeric(length=nyrs+1)
  sbt <- numeric(length=nyrs+1) # can eliminate this eventually. It's the same as biomass in the dd model
  # FIXME: I think this should just be from syr:(nyr-sage)
  log_rt <- numeric(length=nyrs) # estimated log recruits
  annual_mean_wt <- numeric(length=nyrs)
  # Add recruitment for projection year ... assume it is average
  rnplus=exp(log_avgrec)

  # Initial Conditions
  snat <- exp(-m) # natural survival rate
  # Eqm mean weight
  wbar <- (snat*alpha.g + wk*(1-snat))/(1-snat*rho.g) # RF: calculation checked against rep file

  # Unfished numbers and biomass
  no <- ro/(1-snat)
  bo <- no*wbar # RF: calculation checked against rep file
  # added sbo for delay diff model - delete later
  sbo <- bo

  # Parameters of Stock-Recruit relationship
  # Maximum juvenile survival rate (same for BH and Ricker)
  alpha.sr <- CR*(ro/bo) # NOTE: *Steve called this so*. RF thinks this is confusing notation because the o implies unfished, whereas this is max juv survival rate when the biomass is close to zero

  # Beta parameter (density-dependence in juvenile survival)
  if(ctl$misc[2]==1){
    # Beverton-Holt
    beta.sr <- (CR -1)/bo
  }
  if(ctl$misc[2]==2){
    # Ricker
    beta.sr <- log(CR)/bo
  }

  # Three options for how first year is initialized (set in ctrl$misc[5]):
  # 0. Unfished and not at equilibrium
  # 1. Unfished equilibrium
  # 2. Fished and at equilibrium with Ft in first year
  tmp_nAge <- numeric(length(sage:nage)) # Numbers at age in first year
  if(ctl$misc[5]==0){
    # Unfished and not at equilibrium - initialize with age structure as in ASM
    # Set up a vector of n-at-age for the first year (length sage:nage) with devs
    # Decay exp(log_recinit + init_log_rec_devs[j-1]) by exp(-M) to fill n-at-age in first year

    # First age is exp(log_avgrec + log_rec_devs[1]) because this is recruits in year 1
    tmp_nAge[1] <- exp(log_avgrec + log_rec_devs[1])

    # Now fill in the other numbers at age for year 1

    for(j in 2:length(tmp_nAge)){
      tmp_nAge[j] <- exp(log_recinit + init_log_rec_devs[j-1]) * exp(-Mt[1]*(j-1))
    }
    tmp_nAge[j] <- (exp(log_recinit + init_log_rec_devs[j-1]) * exp(-Mt[1]*(j-1)))/(1 - exp(-Mt[1])) # plus group

    ## TEST (see iscam calcNumbersAtAge function) - make sure tr is the same as log(tmp_nAge)
    ## DELETE THIS IN LATER VERSIONS
    ## A little tricky  because ADMB indexes actual ages or years rather than indexes
    ## **YES, calculations are identical**
    # tr <- vector(length=length(sage:nage)) # iscam ASM's version of log(tmp_n_Age)
    # lx <- vector(length=length(sage:nage)) # unfished survivorship at age
    # lx[1] <- 1
    # for(j in 2:length(sage:nage)){
    #   lx[j] = lx[j-1] * exp(-Mt[1])
    # }
    # lx[j] <- lx[j]/(1-exp(-Mt[1]))
    #
    # tr[1]  <- log_avgrec+log_rec_devs[1]
    # tr[2:length(tr)] <- log_recinit +init_log_rec_devs
    # tr[2:length(tr)] <- tr[2:length(tr)]+log(lx[2:length(tr)])
    # exp(tr)

    # Add up the numbers at age to get numbers in year 1
    # RF confirmed numbers, biomass and mean wt calculations with rep file
    numbers[1] <- sum(tmp_nAge)
    biomass[1] <- sum(tmp_nAge*d3_wt_avg) # d3_wt_avg calculated in data section from vonB parameters
    annual_mean_wt[1] <- biomass[1]/numbers[1]

    #  Initialise log recruits
    log_rt[1] <- log_avgrec+log_rec_devs[1] # this is just the same as log(tmp_nAge[1])
  }
  if(ctl$misc[5]==1){
    # Unfished equilibrium - **not tested for P. cod**
    numbers[1] <- no
    biomass[1] <- bo
    annual_mean_wt[1] <- biomass[1]/numbers[1]
    log_rt[1] <- log(ro) # Shouldn't this be plus log_rec_devs[1]???
  }
  if(ctl$misc[5]==2){
    # Fished at equilibrium with Ft in first year **not tested for P. cod**
    # NEED TO CHECK THESE CALCS. SOURCE IS 2004 P.COD ASSESSMENT OR H&W 1992 - need to check
    # RF: I think I have some code from SJDM somewhere with these calcs
    sfished = surv[1] # equilibrium survivorship at initial fishing mortality (gear 1 commercial fishery)
    annual_mean_wt[1] = (sfished*alpha.g + wk*(1.-sfished))/(1-rho.g*sfished)

    biomass[1] = -(annual_mean_wt[1]*(wk*alpha.sr-1)
                   +sfished*(alpha.g + rho.g * annual_mean_wt[1]))/
      (beta.sr*(sfished*alpha.g + sfished*rho.g*annual_mean_wt[1]-
                  annual_mean_wt[1]))
    numbers[1] = biomass[1]/annual_mean_wt[1]

    # log rt originally missing from this option
    # chose log_avgrec as placeholder-- dangerous if fishing in first year and before was very high.
    log_rt[1] <- log_avgrec[1] #+ log_rec_devs[1]  # Shouldn't this be plus log_rec_devs[1]???
  }
  if(!ctl$misc[5] %in% 0:2){
    stop("Starting conditions not specified! ctl$misc[5] must be set to 0, 1 or 2 (see ?pcod2020ctl for definitions).")
  }

  # Set sbt[1] <- biomass[1]
  # Eventually can delete this duplication but for now will be easier
  # bc sbt appears elsewhere in code as hangover from ASM
  sbt[1] <- biomass[1]

  # Time dynamics
  for(i in 2:nyrs){

    # Update recruits
    log_rt[i] <- log_avgrec +log_rec_devs[i]

    # Update biomass and numbers
    biomass[i] <- surv[i-1]*(rho.g*biomass[i-1]+alpha.g*numbers[i-1]) +wk*exp(log_rt[i]) # eq. 9.2.5 in HW
    numbers[i] <- surv[i-1]*numbers[i-1]+exp(log_rt[i])
    annual_mean_wt[i] <- biomass[i]/numbers[i]		# calculate predicted weight in dynamics - possible option to fit to it
    sbt[i] <- biomass[i]
  }

  # RF confirmed numbers, biomass and mean wt calculations with rep file
  biomass[nyrs+1]  <- (surv[nyrs]*(rho.g*biomass[nyrs]+alpha.g*numbers[nyrs]) + wk*rnplus)
  numbers[nyrs+1]  <- surv[nyrs]*numbers[nyrs]+rnplus
  sbt[nyrs+1]      <- biomass[nyrs+1] # set spawning biomass to biomass

  # End calcNumbersBiomass_deldiff
  #|---------------------------------------------------------------------|

  #|---------------------------------------------------------------------|
  # 4. calcFisheryObservations_deldiff()
  # Purpose: This function calculates commercial catches for each year and gear
  # The catch matrix in the data file has 4 columns:
  #   year, gear, type (1=catch in weight (tonnes); 2=catch in numbers), catch
  # nctobs is from the dat list and indicates the total number of catch obs

  # Set up predicted catch by gear and year
  # like log_ft_pars, this is one long vector of all catches from all gears
  ct <- numeric(length=nctobs) # predicted catch
  eta <- numeric(length=nctobs) # catch residuals
  tmp <- numeric(length=nctobs) # holder for likelihood components

  for(ii in 1:nctobs){
    # Set up counters
    iyear <- catch[ii,1] # actual year
    i <- as.integer(year_lookup[which(year_lookup[,1]==iyear),2]) # year index
    k <- catch[ii,2] # gear
    mm <- catch[ii,3] # type: 1=catch in weight; 2=catch in numbers
    d_ct <- catch[ii,4] # observed catch

    # Baranov catch equation
    if(mm==1) {
      # catch in weight
      ct[ii] = (ft[k,i]/(Mt[i] + ft[k,i]))*(1-exp(-Mt[i]-ft[k,i]))*biomass[i]

    }
    if(mm==2){
      # catch in numbers
      ct[ii] = (ft[k,i]/(Mt[i] + ft[k,i]))*(1-exp(-Mt[i]-ft[k,i]))*numbers[i]
    }
    if(!mm %in% 1:2){
      stop("Catch type must be 1 (weight) or 2 (numbers). Set in column 3 of dat$catch")
    }

    # catch residual
    eta[ii] = log(d_ct+TINY) - log(ct[ii]+TINY)

    # LIKELIHOOD COMPONENT
    tmp[i] <- dnorm(log(d_ct+TINY), log(ct[ii]+TINY), sig_c, log=T) # Yes. -134.716 Matches nlvec_dd in rep file

    # NOTE: Catches are not exactly as in rep file - could be rounding in the reported
    # log_ft_pars from the par file. GO BACK AND CHECK ALL CALCS AND VALUES
  } # end for ii

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
  # LIKELIHOOD FOR CATCH OBS
  nlvec_dd_ct <- sum(tmp)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|

  # End calcFisheryObservations_deldiff
  #|---------------------------------------------------------------------|

  #|---------------------------------------------------------------------|
  # 5. calcSurveyObservations_deldiff()
  # Purpose: This function calculates predicted survey observations for each year

  # Needed to determine if q is random walk
  q_prior <- q_control$priortype

  # Set up vector for mle qs (per Walters&Ludwig 1993 https://doi.org/10.1139/f94-07)
  q <- numeric(length=nit) # vector of q for each survey

  # Get weights for weighting individual observations for likelihood
  # iscam makes a ragged matrix of the it weights (see L445 of devs/iscam.tpl)
  # then weights by a global mean
  # We can just make a vector since we just need it for the mean
  it_wt <- dat$indices[[1]][,4]
  for(kk in 2:nit){
    tmp <- dat$indices[[kk]][,4]
    it_wt <- c(it_wt, tmp)
  }
  # global mean of it_wts
  it_wt_mean <- mean(it_wt) #devs/iscam.tpl L463

  # set up lists for storing residuals and predicted indices
  # these are ragged arrays in iscam
  epsilon <- list() # residuals
  it_hat <- list() # predicted indices
  qt <- list() # time varying q (if ctl$surv.q prior type == 2)

  #loop over surveys to create list
  for(kk in 1:nit){
    epsilon[[kk]] <- numeric(length=nitnobs[kk])
    it_hat[[kk]]  <- numeric(length=nitnobs[kk])
    qt[[kk]]  <- numeric(length=nitnobs[kk])
    epsilon[[kk]][1:nitnobs[kk]] <- 0.
    it_hat[[kk]][1:nitnobs[kk]] <- 0.
    qt[[kk]][1:nitnobs[kk]] <- 0.
  }

  # Now loop through gears and index obs to predict survey obs
  # Then calculate likelihood component
  for(kk in 1:nit){
    # Vulnerable numbers or biomass
    V <- numeric(length=nitnobs[kk])
    V[1:nitnobs[kk]] <- 0.

    # not sure why we need these
    #nz <- 0 # counter for number of observations
    #iz <- 1  # index for first year of data for prospective analysis

    for(ii in 1:nitnobs[kk]){
      iyear    <- indices[[kk]][ii,1] # actual year (not used)
      i <- as.integer(year_lookup[which(year_lookup[,1]==iyear),2]) # year index - need this to match biomass to survey obs
      k    <- indices[[kk]][ii,3] # gear
      di   <- indices[[kk]][ii,5] # timing

      #nz <- nz+1  # counter for number of observations

      # TODO: check these correctly matched up with years - yes, checked against rep file
      z = ft[k,i]+Mt[i]
      # Adjust numbers and biomass for survey timing. If di=0, no adjustment
      Np = numbers[i] * exp( -z * di)
      Bp = biomass[i] * exp( -z * di)

      # Two different survey types: 1=prop to numbers; 2=prop to biomass
      if(survtype[kk]==1){
        V[ii] <- Np
      }
      if(survtype[kk]==2){
        V[ii] <- Bp
      }
      if(!survtype[kk] %in% 1:2){
        stop(("Survey type must be 1 (survey proporional to numbers) or 2 (survey proporional to biomass). Set in dat$survtype"))
      }
    } #end of ii loop

    it 	<- t(indices[[kk]])[2,1:nitnobs[kk]] # index
    wt 	<- t(indices[[kk]])[4,1:nitnobs[kk]] # index weight
    wt 	<- wt/sum(wt) # normalized weight (sum to 1) # Q: here the weights are normalized by the sum but in the likelihood, iscam normalises it_wt by the mean (in the data section)

    # get mle q (same as for ASM, from Walters&Ludwig 1993 https://doi.org/10.1139/f94-07)
    zt 	  <- log(it) - log(V[1:nitnobs[kk]])
    zbar  <- sum(zt*wt)
    q[kk] <- exp(zbar)

    # survey residuals - checked against rep file
    epsilon[[kk]][1:nitnobs[kk]] <- zt - zbar
    it_hat[[kk]][1:nitnobs[kk]] <- q[kk] * V[1:nitnobs[kk]]

    # SPECIAL CASE: penalized random walk in q.
    # !!!NOT TESTED!!! This is dimensioned correctly but have not checked calcs
    # if(q_prior[kk]==2 ){
    #   epsilon[[kk]][1:nitnobs[kk]] <- 0 # initialize epsilon
    #
    #   # iscam ADMB code:
    #   # fd_zt <- first_difference(zt)
    #   # From the admb source code, looks like first_difference returns a vector of
    #   # i+1 - i:
    #   # i.e., differences.elem(i) = values.elem(i + 1) - values.elem(i);
    #   fd_zt <- diff(zt)
    #   zw_bar <- sum(fd_zt*wt[1:(nitnobs[kk]-1)])
    #   epsilon[[kk]][1:(nitnobs[kk]-1)] <- fd_zt - zw_bar
    #   qt[[kk]][1] = exp(zt[1])
    #
    #   for(ii in 2:nitnobs[k]){
    #     qt[[kk]][ii] = qt[[kk]][ii-1] * exp(fd_zt[ii-1])
    #   }
    #   it_hat[[kk]][1:nitnobs[kk]] = qt[[kk]][1:nitnobs[kk]]*V[1:nitnobs[kk]]
    # }

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
    # LIKELIHOOD FOR RELATIVE ABUNDANCE INDICES
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
    # Normalise the weights for survey kk
    # Side note: iscam normalizes it_wt in the data section
    #   by dividing by the global mean. But above, where q is calculated,
    #   the weights (called wt) are normalized by dividing by sum
    it_wt <- dat$indices[[kk]][,4]/it_wt_mean #devs/iscam.tpl L466
    # vector for weights for each obs in survey k
    sig_it <- sig/it_wt

    # TO DELETE: Old version using admb statslib function
    # NOTE: this gives positive results! dnorm below returns a negative
    # tmp <- admb_dnorm_vector_vector(epsilon[[kk]], sig_it)
    # nlvec_dd_it[kk] <- tmp

    nlvec_dd_it[kk] <- sum(dnorm(zt,zbar,sig_it,log=T))

  } #end kk loop

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
  # End calcSurveyObservations_deldiff
  #|---------------------------------------------------------------------|

  #|---------------------------------------------------------------------|
  # 6. calcStockRecruitment_deldiff()
  # Purpose:
  # This function is used to derive the underlying stock-recruitment
  # relationship that is ultimately used in determining MSY-based reference
  # points.  The objective of this function is to determine the appropriate
  # Ro, Bo and steepness values of either the Beverton-Holt or Ricker  Stock-
  #   Recruitment Model:
  #  Beverton-Holt Model
  #  Rt=k*Ro*St/(Bo+(k-1)*St)*exp(delta-0.5*tau*tau) \f$
  #
  #   Ricker Model
  #  Rt=so*St*exp(-beta*St)*exp(delta-0.5*tau*tau) \f$

  # Set up vectors
  # Don't need to declare first two
  #rt     <- numeric(length=nyrs-sage)  # estimated recruits from calcNumbersBiomass_deldiff()
  #delta  <- numeric(length=nyrs-sage)  # residuals between estimated R and R from S-R curve (process err)
  tmp_rt <- numeric(length=nyrs) # recruits derived from stock-recruit model
  tmp_rt[1:nyrs] <- 0

  # get the process error term from the errors in variables parameters
  # [For the P. cod assessment, rho and varphi are set to give tau=0.8 and obs error=0.2]
  tau <- sqrt(1-rho)*varphi

  # counter
  iicount <- 0

  for(i in 1:nyrs){
    iicount <- iicount+1

    if(ctl$misc[2]==1){
      # Beverton-Holt recruitment
      if(iicount <= kage){
        tmp_rt[i] =  alpha.sr*sbt[1]/(1.+beta.sr*sbt[1])
      }else{
        tmp_rt[i] = alpha.sr*sbt[i-kage]/(1.+beta.sr*sbt[i-kage])
      }
    }
    if(ctl$misc[2]==2){
      # Ricker recruitment
      if(iicount <= kage){
        tmp_rt[i] =  alpha.sr*sbt[1]*exp(-beta.sr*sbt[1])
      }else{
        tmp_rt[i] = alpha.sr*sbt[i-kage]*exp(-beta.sr*sbt[i-kage])
      }
    }
    if(!ctl$misc[2] %in% 1:2){
      stop("Recruitment model must be 1 (Beverton-Holt) or 2 (Ricker). Set in ctl$misc[2]")
    }
  } # end year loop i

  # estimated recruits from calcNumbersBiomass_deldiff()
  # RF Checked against rep file
  rt <- exp(log_rt[(sage+1):nyrs])

  # Calculate delta: process errors (deviations from S-R function)
  # RF: ADD A NOTE WHY BIAS CORRECTION ADDED HERE
  delta <- log(rt)-log(tmp_rt[(sage+1):nyrs])+0.5*tau*tau # RF Checked against rep file

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
  # LIKELIHOOD FOR RECRUITMENT
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
  nlvec_dd_rt <- sum(dnorm(log(rt), log(tmp_rt[(sage+1):nyrs])-0.5*tau*tau,tau, log=T))

  # End calcStockRecruitment_deldiff
  #|---------------------------------------------------------------------|

  #|---------------------------------------------------------------------|
  # 7. calcAnnualMeanWeight_deldiff() //RF added this for P cod - only gets added to objective function if cntrl(15)==1
  # Purpose: This function calculates the mean weight of the catch for each year, gear by dividing the total
  #          biomass by the total numbers

  # Set up lists for obs and predicted annual mean weights
  # these are ragged arrays in iscam
  annual_mean_weight     <- list() # observed
  obs_annual_mean_weight <- list() # predicted
  epsilon_mean_weight <- list() # residuals

  #loop over series to create list
  for(kk in 1:nmeanwt){
    annual_mean_weight[[kk]]     <- numeric(length=nmeanwtobs[kk])
    obs_annual_mean_weight[[kk]] <- numeric(length=nmeanwtobs[kk])
    epsilon_mean_weight[[kk]]    <- numeric(length=nmeanwtobs[kk])
    annual_mean_weight[[kk]][1:nmeanwtobs[kk]] <- 0.
    obs_annual_mean_weight[[kk]][1:nmeanwtobs[kk]] <- 0.
    epsilon_mean_weight[[kk]][1:nmeanwtobs[kk]] <- 0.
  }

  # loop through series with empirical annual mean weight data
  for(kk in 1:nmeanwt){

    Vn <- numeric(length = nmeanwtobs[kk])	      # Vulnerable numbers to gear
    Vb <- numeric(length = nmeanwtobs[kk])	      # Vulnerable biomass to gear
    Vn[1:nmeanwtobs[kk]] <- Vb[1:nmeanwtobs[kk]] <- 0

    # Loop through observations
    for(ii in 1:nmeanwtobs[kk]){

      iyear    <- meanwtdata[[kk]][ii,1]  # actual year - not used
      i <- as.integer(year_lookup[which(year_lookup[,1]==iyear),2]) # year index - need this to match meanwt to obs
      k    <- meanwtdata[[kk]][ii,3]  # gear
      di   <- meanwtdata[[kk]][ii,4]  # timing

      ws  = exp(-Zt[i]*di)   # Total mortality that accounts for timing
      wN  = numbers[i]*ws
      wB  = biomass[i]*ws
      Vn[ii] <- wN
      Vb[ii] <- wB

      # RF Checked against rep file
      annual_mean_weight[[kk]][ii] <- Vb[ii]/Vn[ii]
      obs_annual_mean_weight[[kk]][ii] <-  meanwtdata[[kk]][ii,2]	  # fill a list of vectors with observed annual mean weights
      #  residual
      epsilon_mean_weight[[kk]][ii] <-  log(annual_mean_weight[[kk]][ii]) - log(obs_annual_mean_weight[[kk]][ii])
    }	# end ii loop

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
    # LIKELIHOOD FOR ANNUAL MEAN WEIGHT
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
    nlvec_dd_wt <- sum(dnorm(log(annual_mean_weight[[kk]]),log(obs_annual_mean_weight[[kk]]), sig_w, log=T)) # 3.407 Yes. Matches nlvec_dd in rep file.
  } # end kk loop

  # End calcAnnualMeanWeight_deldiff
  #|---------------------------------------------------------------------|

  #|---------------------------------------------------------------------|
  # calcObjectiveFunction();

  # !!Distributions that match those in ADMB are
  #   provided in the file R\likelihood_funcs.R!! (sourced above)

  #==============================================================================================
  # ~Data and recruitment~ (nlvec_dd)
  #==============================================================================================

  ## ~~ September 20 2024:
  ## ~~ MOVED LIKELIHOOD EQNS INTO RELEVANT MAIN MODEL SECTIONS ~~
  # Likelihood for catch
  # TO DELETE: Old version using admb statslib function
  #tmp <- admb_dnorm_vector_const(eta, sig_c) # Yes. -134.716 Matches nlvec_dd in rep file.

  # Moved likelihood for survey indices. Too long to leave here.

  # Likelihood for recruitment
  # tmp <- admb_dnorm_vector_const(delta, tau)
  # nlvec_dd_rt <- tmp # 82.9127 Yes. Matches nlvec_dd in rep file.

  # Likelihood for mean weight
  # We are entering the likelihood in log space here - do we need a Jacobian transformation?
  # for(kk in 1:nmeanwt){
  #   tmp <- admb_dnorm_vector_const(epsilon_mean_weight[[kk]], sig_w)
  #   nlvec_dd_wt <- tmp # 3.407 Yes. Matches nlvec_dd in rep file.
  # }

  #==============================================================================================
  # ~PRIORS~ (priors)
  #==============================================================================================
  # Leading parameters
  # !The uniform is returning a negative value!
  for(i in 1:length(theta)){
    ptype <- theta_control$prior[i] # prior type

    if(theta_control$phz[i]>=1){
      # Uniform
      if(ptype==0){
        #ptmp <- log(1./(theta_control$p2[i]-theta_control$p1[i])) # Note, iscam used the bounds not p1 and p2
        # For testing use the same as iscam. Why is the uniform set up like this?
        ptmp <- log(1./(theta_control$ub[i]-theta_control$lb[i])) # Note, iscam used the bounds not p1 and p2
      }
      # Normal
      if(ptype==1){
        ptmp <- admb_dnorm_const_const(theta[i],theta_control$p1[i],theta_control$p2[i])
      }
      # Lognormal
      if(ptype==2){
        ptmp <- admb_dlnorm_const_const(theta[i],theta_control$p1[i],theta_control$p2[i])
      }
      # Beta
      if(ptype==3){
        # Constrain between 0.2 and 1
        trans <- (theta[i]-theta_control$lb[i])/(theta_control$ub[i]-theta_control$lb[i])
        ptmp <- admb_dbeta_const_const(trans, theta_control$p1[i],theta_control$p2[i])
      }
      # Gamma
      if(ptype==4){
        ptmp <- admb_dgamma_const_const(theta[i],theta_control$p1[i],theta_control$p2[i]);
      }
      priors[i] <- ptmp
    } # end if
  }# end i

  # Catchability coefficients q
  for(kk in 1:nit){
    if(q_control$priortype[kk] == 1){
      qtmp <- admb_dnorm_const_const(log(q[kk]), q_control$priormeanlog[kk], q_control$priorsd[kk])
      qvec[kk] <- qtmp
    }
  }

  # #==============================================================================================
  # # ~PENALTIES~ pvec
  # # ---------------------------------------------------------------------------------|
  # #  LIKELIHOOD PENALTIES TO REGULARIZE SOLUTION
  # # ---------------------------------------------------------------------------------|
  # #  NOTE: iscam includes another element pvec(2) for m deviations and
  # # also has an empty spot in pvec(3)
  # #  pvec(1)  -> penalty on mean fishing mortality rate.
  # #  pvec(2)  -> penalty on recruitment deviations.
  # #  pvec(3)  -> penalty on initial recruitment vector.
  # #  pvec(4)  -> constraint to ensure sum(log_rec_dev) = 0
  # #  pvec(5)  -> constraint to ensure sum(init_log_rec_dev) = 0
  # #==============================================================================================
  #
  # Fishing mortality
  nft <- length(log_ft_pars)
  mean_log_ft_pars <- sum(log_ft_pars)/nft # getting the mean manually prevents lost class attribute error
  pvec[1] <- admb_dnorm_const_const(mean_log_ft_pars,log(mean_f),sig_f) # Note, there are no phases in rtmb so use last phase settings - might mess up estimation

  # Penalty for log_rec_devs and init_log_rec_devs (large variance here)
  bigsd <- 2. # possibly put this in the data

  pvec[2] <- admb_dnorm_vector_const(log_rec_devs, bigsd)
  pvec[3] <- admb_dnorm_vector_const(init_log_rec_devs, bigsd)

  #constrain so that sum of log_rec_dev and sum of init_log_rec_dev = 0
  ndev <- length(log_rec_devs) # getting the mean manually prevents lost class attribute error
  meandev <- sum(log_rec_devs)/ndev # this was s in iscam
  pvec[4] <- 1.e5*meandev*meandev #mean(log_rec_devs)*mean(log_rec_devs)

  nidev <- length(init_log_rec_devs) # getting the mean manually prevents lost class attribute error
  meanidev <- sum(init_log_rec_devs)/nidev # this was s in iscam
  pvec[5] <- 1.e5*meanidev*meanidev

  # joint likelihood, priors and penalties
  objfun <- nlvec_dd_ct +
    sum(nlvec_dd_it) +
    nlvec_dd_rt +
    nlvec_dd_wt +
    sum(priors) +
    sum(qvec) +
    sum(pvec)

  # if(test==T){
  #   # just for testing likelihood coded correctly. Delete after testing
  #   objfunlist <- list()
  #   objfunlist$objfun <- objfun
  #   objfunlist$nlvec_dd_ct <- nlvec_dd_ct
  #   objfunlist$nlvec_dd_it <- nlvec_dd_it
  #   objfunlist$nlvec_dd_rt <- nlvec_dd_rt
  #   objfunlist$nlvec_dd_wt <- nlvec_dd_wt
  #   objfunlist$priors <- priors
  #   objfunlist$qvec <- qvec
  #   objfunlist$pvec <- pvec
  #   print(objfunlist)
  # } # end if test
  # End calcObjectiveFunction
  #|---------------------------------------------------------------------|

  # REPORT_SECTION
   # Need to find a dynamic way of reporting it_hat and annual_mean_weight
  # NOTE: ADREPORT keeps track of standard error. REPORT does not and is quicker.
  # Can't ADREPORT lists or use loops
  # Hardwire for now
  it_hat_all <- c(it_hat[[1]],it_hat[[2]],it_hat[[3]],it_hat[[4]],it_hat[[5]])
  annual_mean_weight_all <- annual_mean_weight[[1]]

  # REPORTS (no std error - obj function components and the MCMC posteriors we want to see)
  # Includes some variables we need to pass to projection model
  # Objective function components
  REPORT(nlvec_dd_ct)
  REPORT(nlvec_dd_it)
  REPORT(nlvec_dd_rt)
  REPORT(nlvec_dd_wt)
  REPORT(priors)
  REPORT(qvec)
  REPORT(pvec)
  # Quantities needed for reporting and for the projection model
  # Biomass, recruits, q ... get R0, M, h, F and recdevs directly from tmbstan
  REPORT(biomass)
  REPORT(numbers)
  REPORT(Ft)
  REPORT(rt)
  REPORT(q)
  REPORT(surv)
  REPORT(tau)
  REPORT(alpha.sr)
  REPORT(beta.sr)

  # ADREPORT, includes standard error
  ADREPORT(biomass)
  ADREPORT(numbers)
  ADREPORT(rt)
  ADREPORT(delta)
  ADREPORT(q)
  # Predicted Catch, Indices and Annual Mean Weight
  ADREPORT(ct)
  # Figure out how to report these lists without creating errors
  ADREPORT(it_hat_all)
  ADREPORT(annual_mean_weight_all)

  # RETURN OBJECTIVE FUNCTION VALUE
  objfun # return joint neg log likelihood
} # end model
