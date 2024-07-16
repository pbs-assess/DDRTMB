# Likelihood functions that emulate what ADMB is doing for iscam
# Author: Robyn Forrest
# Date created: July 16 2024
# Date last modified: July 16, 2024

# Emulate what admb is doing when given a vector of residuals
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
  if(length(std)!=n) stop("Residuals and st devs are different lengths.")
  var <- std^2
  SS <- resid^2
  negloglike <- 0.5*n*log(2.*pi)+sum(log(std))+sum((SS/(2.*var)))
  negloglike
}
