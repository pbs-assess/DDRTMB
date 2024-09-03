#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GET FT FOR THE DELAY DIFFERENCE MODEL FROM BARANOV EQUATION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get_ftdd (based on iscam get_ftdd (see devs\baranov.cpp)

# Uses Newton-Raphson algorithm to solve for F, given a Catch in a given year

# arguments:
# ct = catch
# m = natural mortality
# b = biomass

# Returns:
# a single value of fishing mortality

get_ftdd <- function(ct, m, b){

  if(ct<b){
    # Initial guess for F
    ft=ct/(b*exp(-m/2.))

    # now do iterative numerical search for ft. Should converge quickly.
    for(ii in 1:50){
      f <- ft
      z <- m+f
      s <- exp(-z)
      mort <- 1.-s

      # terms needed for derivative
      t1 <- f/z
      t2 <- t1*mort
      t3 <- mort*b
      # predicted catch
      pct <- t2*b

      # derivative of catch wrt ft (simple calculus using chain rule)
      dct <- (t3/z)  - (f*t3)/z^2 + (t1*s)*b

      # Newton step
      ft = ft - (pct-ct)/dct
    }# end for ii
  }else{
    # if ct > b, set F very high
    ft <- 20.
  } # end ifelse

  return(ft)
}


