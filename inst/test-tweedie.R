setwd("inst")

## Load TMB TRUE
library(TMB)

## Make C++ file
# TMB::template("tweedie_test.cpp")

## Compile and load the model
compile("tweedie_test.cpp")
dyn.load(dynlib("tweedie_test"))

## Data and parameters
set.seed(1)
data <- list(y = tweedie:::rtweedie(1000, mu=2, phi=2, p=1.5))
parameters <- list(mu=0.8, phi=0.6, p=1.1)

## Make a function object
obj <- MakeADFun(data, parameters, DLL="tweedie_test")

## Call function minimizer
opt <- nlminb(obj$par, obj$fn, obj$gr)

## Get parameter uncertainties and convergence diagnostics
sdr <- sdreport(obj)
sdr

(p <- opt$par[["p"]])
(phi <- opt$par[["phi"]])
(mu <- opt$par[["mu"]])
-sum(fishMod::dTweedie(y = data$y, mu = mu, p = p, phi = phi, LOG = TRUE))

opt$objective
