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
  X = runif(3000), Y = runif(3000),
  year = 1
)
mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.05)
n_s <- nrow(mesh$mesh$loc)

sim_dat <- sdmTMB_simulate(
  formula = ~ 1,
  data = predictor_dat,
  time = "year",
  mesh = mesh,
  family = gaussian(),
  range = 0.1,
  sigma_E = 0.0,
  phi = 0.04,
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
parameters <- list(B0 = 0, ln_tau_O = 0, ln_kappa = 0, ln_sigma=0,
                   omega_s = matrix(0,n_s, 1))
obj1 <- MakeADFun(data, parameters, random = "omega_s", DLL = "pc_matern", hessian = TRUE)
fit1 <- nlminb(obj1$par, objective = obj1$fn, gradient = obj1$gr,
              control=list(eval.max=10000, iter.max=10000))
fit1

# with prior -- this estimates B0 well, but struggles with ln_tau_O and kappa
compile("inst/pc_matern_prior.cpp")
dyn.load(dynlib("inst/pc_matern_prior"))

data <- dat
data$matern_range <- 0.3 # greater than this is range_prob
data$range_prob <- 0.05
data$matern_SD <- 0.01 # less than this is SD_prob
data$SD_prob <- 0.05
obj2 <- MakeADFun(data, parameters, DLL = "pc_matern_prior", random = "omega_s", hessian = TRUE)
fit2 <- nlminb(obj2$par, objective = obj2$fn, gradient = obj2$gr,
              control=list(eval.max=4000, iter.max=4000))
fit2

sdr1 = sdreport(obj1)
sdr2 = sdreport(obj2)

# look at predictions
sim_dat$pred_1 = sdr1$value[which(names(sdr1$value)=="pred")]
sim_dat$pred_2 = sdr2$value[which(names(sdr2$value)=="pred")]

# predictions are perfectly correlated with / without prior
ggplot(sim_dat, aes(pred_1, pred_2)) +
  geom_point(alpha=0.1, col="blue")

# need to add 3rd case with Jacobian adj
