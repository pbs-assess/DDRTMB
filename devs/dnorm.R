# Comparing R probability distributions with ADMB implementations
# Author: RF
# Date: July 3 2024
# ADMB source for dnorm (by Steve Martell):
#   https://github.com/admb-project/admb/blob/main/contrib/statslib/dnorm.cpp
# dnorm is overloaded in ADMB to take different types of input (e.g., vector or variable)

# Notes from ADMB dnorm.cpp file:
# This file contains the negative loglikelihood
# functions for the normal distribution and is implemented to be consistent
# with the statistical program R with log=TRUE. The function
# dnorm is overloaded to accomodate single variables and vectors.
#
# There are also overloaded versions where the user can specify the likelihood (i.e., log=FALSE)
#
# The function is implemented as the negative log of the normal density function:

#   0.5ln(2pi) + ln(sigma) + 0.5frac{(x-mu)^2}{sigma^2}

#where mu is the mean and sigma is the standard deviation.

#The concentrated likelihood is implemented as:
#   0.5nln(sum_{i=1}^{n}epsilon^2)
# where epsilon is a vector of residuals with an assumed mean 0.
#

# Test
mu <- 0.5 # prior mean
std <- 0.1 # prior sd
x <- 0.4 # estimate

# first ADMB version
dnormADMB <- 0.5*log(2.*pi)+log(std)+(0.5*(x-mu)^2)/(std*std)
dnormADMB

dnormR <- dnorm(x,mu,std,log=TRUE)
dnormR

# OK, they return the same but it is negative in ADMB.
# So we do need to reverse the sign of jnll in DDRTMB after adding
# all the likelihood components together

# FOR RESIDUALS (i.e., rec devs)
# iscam has pvec(4) += dnorm(log_rec_devs,2.0) where 2.0 is the std

# dvariable dnorm( const dvar_vector& residual, const prevariable& std )
# {
#   RETURN_ARRAYS_INCREMENT();
#   long n=size_count(residual);
#   dvariable SS=norm2(residual);
#   dvariable tmp=n*(0.5*log(2.*M_PI)+log(std))+0.5*SS/(std*std);
#   RETURN_ARRAYS_DECREMENT();
#   return( tmp );
# }

# norm2 (from https://api.admb-project.org/group__matop.html#ga138dd90f1e3f89d4d4bc0b5d675c8eb4)
# Squared norm of a vector; constant objects.
# Computes the sum of squares of its vector argument.
# Parameters:
# t1	A vector, $v$.
# Returns
# $||v||^2 = v\times v = \sum_i v_i^2$
#   Definition at line 28 of file dvect3.cpp.
# Returns
# v^2 = v * v = sum_i(v_i^2)
#   Definition at line 28 of file dvect3.cpp.

# I think this is what we do in R
n <- 100
std <- 2
residual <- rnorm(n,0,1)
SS <- sum(residual^2)
tmp <- n*(0.5*log(2.*pi)+log(std))+0.5*SS/(std*std) # this is a large positive number - shouldn't it be negative?
