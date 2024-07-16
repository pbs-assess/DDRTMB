# Likelihood functions that emulate what ADMB is doing for iscam
# Author: Robyn Forrest
# Date created: July 16 2024
# Date last modified: July 16, 2024

# Emulate what admb is doing when given a vector of residuals



# ~~~ dnorm (statslib/dnorm.cpp) ~~~
# dnorm with a constant estimate, constant mean and constant standard deviation
# See https://github.com/admb-project/admb/blob/dd6ccb3a46d44582455f76d9569b012918dc2338/contrib/statslib/dnorm.cpp#L46
admb_dnorm_const_const <- function(x,mu, std){
  negloglike <- 0.5*log(2.*pi)+log(std)+(0.5*(x-mu)^2)/(std*std)
  negloglike
}

# dnorm with a vector of residuals and constant standard deviation
# See https://github.com/admb-project/admb/blob/dd6ccb3a46d44582455f76d9569b012918dc2338/contrib/statslib/dnorm.cpp#L259
admb_dnorm_vector_const <- function(resid, std){
  n <- length(resid)
  SS <- sum(resid^2) # in ADMB: norm2(x-mu);
  negloglike <- n*(0.5*log(2.*pi)+log(std))+0.5*SS/(std*std)
  negloglike
}

# dnorm with a vector of residuals and vector of standard deviations (both size=n)
# See https://github.com/admb-project/admb/blob/dd6ccb3a46d44582455f76d9569b012918dc2338/contrib/statslib/dnorm.cpp#L311
admb_dnorm_vector_vector <- function(resid, std){
  n <- length(resid)
  if(length(std)!=n) stop("Residuals and st devs are different lengths in dnorm function \n")
  var <- std^2
  SS <- resid^2
  negloglike <- 0.5*n*log(2.*pi)+sum(log(std))+sum((SS/(2.*var)))
  negloglike
}

# ~~~ dlnorm (statslib/dlnorm.cpp) ~~~
# dlnorm with a constant estimate, constant mean and constant standard deviation
# See https://github.com/admb-project/admb/blob/dd6ccb3a46d44582455f76d9569b012918dc2338/contrib/statslib/dlnorm.cpp#L36
admb_dlnorm_const_const <- function(x,mu,std){
  negloglike <- 0.5*log(2.*pi)+log(std)+log(x)+(log(x)-mu)^2/(2.*std*std)
  negloglike
}

# ~~~ dbeta (statslib/dbeta.cpp) ~~~
# dbeta with a constant estimate, constant mean and constant standard deviation
# See https://github.com/admb-project/admb/blob/dd6ccb3a46d44582455f76d9569b012918dc2338/contrib/statslib/dbeta.cpp#L37
admb_dbeta_const_const <- function(x,a,b){
  if(x<=0 || x>=1.0)
  {
    stop("x is <=0 or >=1.0 in dbeta function\n")
  }

  if( a<=0 || b <=0 )
  {
    stop("a or b is <= 0 in dbeta function\n")
  }

  negloglike <-  -1.* lgamma(a+b)+(lgamma(a)+lgamma(b))-(a-1.)*log(x)-(b-1.)*log(1.-x)
  negloglike
}

# ~~~ dgamma (statslib/dbeta.cpp) ~~~
# dgamma with a constant estimate, constant mean and constant standard deviation
# See https://github.com/admb-project/admb/blob/dd6ccb3a46d44582455f76d9569b012918dc2338/contrib/statslib/dgamma.cpp#L37
admb_dgamma_const_const <- function(x,a,b){

  if(x<=0)
  {
    stop("x is <=0 in dgamma function\n")
  }

  if( a<=0 || b <=0 )
  {
    stop("a or b is <= 0 in dgamma function\n")
  }

  negloglike <- -a*log(b)+lgamma(a)-(a-1.)*log(x)+b*x
  negloglike

}

#########################################################################################
# TESTING
# gammln Log gamma function
# this returns the same as the lgamma function in R so not needed
# See L52 https://api.admb-project.org/combc_8cpp_source.html
# and https://api.admb-project.org/group__gammafunc.html#gae01374ddda83bc2216ba9df9c76d9995
# Used by dbeta and dgamma
gammln <- function(xx)
{
  cof <- c(76.18009173,-86.50532033,24.01409822,
           -1.231739516,0.120858003e-2,-0.536382e-5)

  x <- xx-1.0
  tmp <- x+5.5
  tmp <- tmp - (x+0.5)*log(tmp)
  ser <- 1.0
  for (j in 1:6)
  {
    x <- x + 1.0
    ser <- ser + cof[j]/x
  }

  return(-tmp+log(2.50662827465*ser))
}


