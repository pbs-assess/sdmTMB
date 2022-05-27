# Simulation proof of when you do/do not want the Jacobian adjustment
# when placing a prior on a transformed parameter
library(sdmTMB)
library(TMB)
library(tmbstan)
options(mc.cores = 1L)

# simualte spatial data from a single year
set.seed(123)

# make fake predictor(s) (a1) and sampling locations:
predictor_dat <- data.frame(
  X = runif(300), Y = runif(300),
  a1 = rnorm(300), year = 1
)
mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)
n_s <- nrow(mesh$mesh$loc)

sim_dat <- sdmTMB_simulate(
  formula = ~ 1,
  data = predictor_dat,
  time = "year",
  mesh = mesh,
  family = gaussian(),
  range = 0.1,
  sigma_E = 0.0,
  phi = 0.001,
  sigma_O = 0.2,
  seed = 42,
  B = c(2) # B0 = intercept, B1 = a1 slope
)
#ggplot(sim_dat, aes(X,Y, col = observed)) + geom_point()

# ln_tau_O = log(1/(0.2*0.2)) = 3.218876
# ln_kappa = log(0.1) = -2.302585
# data
dat <- list(
y = sim_dat$observed,
spde       = mesh$spde$param.inla[c("M0","M1","M2")],
A_st       = mesh$A_st,
A_spatial_index = mesh$sdm_spatial_id - 1L
)

# no prior -- this estimates B0 and ln_sigma well, but struggles with ln_tau_O and kappa
compile("inst/pc_matern.cpp")
dyn.load(dynlib("inst/pc_matern"))

data <- dat
parameters <- list(B0 = 0, ln_tau_O = 0, ln_kappa = 0, ln_sigma=-5,
                   omega_s = matrix(0,n_s, 1))
obj1 <- MakeADFun(data, parameters, DLL = "pc_matern", hessian = TRUE)
fit1 <- nlminb(obj$par, objective = obj$fn, gradient = obj$gr,
              control=list(eval.max=4000, iter.max=4000))
fit1

# with prior -- this estimates B0 well, but struggles with ln_tau_O and kappa
compile("inst/pc_matern_prior.cpp")
dyn.load(dynlib("inst/pc_matern_prior"))

data <- dat
data$matern_range <- 0.3 # greater than this is range_prob
data$range_prob <- 0.05
data$matern_SD <- 0.01 # less than this is SD_prob
data$SD_prob <- 0.05
obj2 <- MakeADFun(data, parameters, DLL = "pc_matern_prior", hessian = TRUE)
fit2 <- nlminb(obj$par, objective = obj$fn, gradient = obj$gr,
              control=list(eval.max=4000, iter.max=4000))
fit2

sdr1 = sdreport(obj1)
sdr2 = sdreport(obj2)
