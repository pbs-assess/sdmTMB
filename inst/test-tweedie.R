setwd("inst")

## Load TMB TRUE
library(TMB)

## Make C++ file
# TMB::template("tweedie_test.cpp")

## Compile and load the model
compile("tweedie_test.cpp")
dyn.load(dynlib("tweedie_test"))

## Data and parameters
data <- list(y = tweedie:::rtweedie(1000, mu=2, phi=2, p=1.5))
parameters <- list(ln_mu=0, ln_phi=0, thetaf=0)

## Make a function object
obj <- MakeADFun(data, parameters, DLL="tweedie_test")

## Call function minimizer
opt <- nlminb(obj$par, obj$fn, obj$gr)

## Get parameter uncertainties and convergence diagnostics
sdr <- sdreport(obj)
sdr

(phi <- exp(opt$par[["ln_phi"]]))
(p <- plogis(opt$par[["thetaf"]]) + 1)
(mu <- exp(opt$par[["ln_mu"]]))

nll <- -sum(fishMod::dTweedie(y = data$y, mu = mu, p = p, phi = phi, LOG = TRUE))
print(nll)
print(opt$objective)

setwd("..")
