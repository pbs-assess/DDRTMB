# Testing stan. Code from Sean Anderson.
library(rstan)
library(tmbstan)
scode <- "
parameters {
  array[2] real y;
}
model {
  y[1] ~ normal(0, 1);
  y[2] ~ double_exponential(0, 2);
}
"
fit1 <- stan(model_code = scode, iter = 10, verbose = FALSE)
print(fit1)

# tmbstan
TMB::runExample("simple")
fit <- tmbstan(obj, chains=1)
