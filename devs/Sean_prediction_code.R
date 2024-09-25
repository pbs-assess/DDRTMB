# Sean Anderson's code for predicting on a fitted model.
# Could be used for projections.
# https://gist.github.com/seananderson/fbd62dc226218a8b031a381854b599a4

library(RTMB)

nll <- function(p) {
  mu <- p$a + p$b * data$x
  nll <- -sum(dnorm(data$Y, mu, exp(p$logSigma), log = TRUE))
  REPORT(mu)
  REPORT(sigma)
  nll
}

set.seed(123)
x <- seq(0, 10, length.out = 10)
data <- list(Y = rnorm(length(x)) + x, x = x)
parameters <- list(a = 0, b = 0, logSigma = 0)
obj <- MakeADFun(nll, parameters)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr
lp <- obj$env$last.par.best
obj$report(lp)
set.seed(1)
obj$simulate()

obj$fn(lp)
obj$gr(lp)

# report() uses the updated data from the global environment:
set.seed(1234)
nd <- data.frame(Y = rep(NA, 11), x = rnorm(11, mean = 5), 11)
data <- nd
obj$report(lp)

# likewise, simulate() uses the new data:
set.seed(1)
obj$simulate()

# but the function and gradient evaluation use the original data:
obj$fn(lp)
obj$gr(lp)
