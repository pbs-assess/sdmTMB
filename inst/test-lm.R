library(TMB)


TMB::compile("inst/test.cpp")
TMB::compile("inst/test_serial.cpp")
dyn.load(dynlib("inst/test"))
dyn.load(dynlib("inst/test_serial"))

## Simulate data
set.seed(123)
x <- seq(0, 10, length = 50001)
data <- list(Y = rnorm(length(x)) + x, x = x)
parameters <- list(a=0, b=0, logSigma=0)

## Fit model
obj <- MakeADFun(data, parameters, DLL="test")
obj_serial <- MakeADFun(data, parameters, DLL="test_serial")
opt <- nlminb(obj$par, obj$fn, obj$gr)

TMB::openmp(n = 2)

bench::mark(opt = nlminb(obj$par, obj$fn, obj$gr), opt_serial = nlminb(obj_serial$par, obj_serial$fn, obj_serial$gr))


d <- subset(pcod, year >= 2011) # subset for example speed
pcod_spde <- make_spde(d$X, d$Y, n_knots = 50) # only 50 knots for example speed
# Tweedie:


cores <- seq(1, 16)

fit <- function(cores) {
  TMB::openmp(n = cores)
  m <- sdmTMB(density ~ 0 + as.factor(year),
    data = d, time = "year", spde = pcod_spde, family = tweedie(link = "log"))
}


TMB::openmp(n = 1)

m <- sdmTMB(density ~ 0 + as.factor(year),
  data = d, time = "year", spde = pcod_spde, family = tweedie(link = "log"))

b <- bench::mark(fit(1), fit(2), fit(3), fit(4), fit(5), fit(6), fit(7), fit(8), fit(9), fit(10), fit(11), fit(12), fit(13), fit(14), fit(15), fit(16), iterations = 1, check = FALSE)

bb <- bench::mark(fit(1), fit(8), iterations = 2, check = FALSE)

bb <- b
b$expression[[1]] <- "01"
b$expression[[2]] <- "02"
b$expression[[3]] <- "03"
b$expression[[4]] <- "04"
b$expression[[5]] <- "05"
b$expression[[6]] <- "06"
b$expression[[7]] <- "07"
b$expression[[8]] <- "08"
b$expression[[9]] <- "09"
b$expression[[10]] <- "10"
b$expression[[11]] <- "11"
b$expression[[12]] <- "12"
b$expression[[13]] <- "13"
b$expression[[14]] <- "14"
b$expression[[15]] <- "15"
b$expression[[16]] <- "16"

plot(b)

TMB::openmp(n = 8)
tictoc::tic()
m <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
  data = d, time = "year", spde = pcod_spde, family = tweedie(link = "log"), cores = 1)
tictoc::toc()

tictoc::tic()
m <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
  data = d, time = "year", spde = pcod_spde, family = tweedie(link = "log"), cores = 8)
tictoc::toc()
