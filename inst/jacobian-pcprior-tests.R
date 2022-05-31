# Simulation proof of when you do/do not want the Jacobian adjustment
# when placing a prior on a transformed parameter
library(sdmTMB)
library(TMB)
library(tmbstan)

# simulate spatial data from a single year
set.seed(123)

# make fake predictor(s) (a1) and sampling locations:
predictor_dat <- data.frame(
  X = runif(800), Y = runif(800),
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

# Model 1: No PC prior
compile("inst/pc_matern.cpp")
dyn.load(dynlib("inst/pc_matern"))

data <- dat
parameters <- list(B0 = 0, ln_tau_O = 0, ln_kappa = 0, ln_sigma=0,
                   omega_s = matrix(0,n_s, 1))
obj1 <- MakeADFun(data, parameters, random = "omega_s", DLL = "pc_matern", hessian = TRUE)
fit1 <- nlminb(obj1$par, objective = obj1$fn, gradient = obj1$gr,
              control=list(eval.max=4000, iter.max=4000))
fit1

# Model 2: Including PC prior
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

# Model 3: Including PC prior and Jacobian adjustment
compile("inst/pc_matern_prior_jacobian.cpp")
dyn.load(dynlib("inst/pc_matern_prior_jacobian"))

obj3 <- MakeADFun(data, parameters, DLL = "pc_matern_prior_jacobian", random = "omega_s", hessian = TRUE)
fit3 <- nlminb(obj3$par, objective = obj3$fn, gradient = obj3$gr,
               control=list(eval.max=4000, iter.max=4000))
fit3

# Look at predictions from all models
sdr1 = sdreport(obj1)
sdr2 = sdreport(obj2)
sdr3 = sdreport(obj3)

# look at predictions
sim_dat$pred_1 = sdr1$value[which(names(sdr1$value)=="pred")]
sim_dat$pred_2 = sdr2$value[which(names(sdr2$value)=="pred")]
sim_dat$pred_3 = sdr3$value[which(names(sdr3$value)=="pred")]

# predictions are perfectly correlated with / without prior
ggplot(sim_dat, aes(pred_1, pred_2)) +
  geom_point(alpha=0.1, col="blue")

cor(sim_dat[,c("pred_1","pred_2","pred_3")])


# sample with Stan -- first model (no prior) throws ~ 200 divergent transitions
options(mc.cores = parallel::detectCores())
#m1 <- tmbstan(obj1, chains = 3, iter = 5000)
m2 <- tmbstan(obj2, chains = 3, iter = 5000)
m3 <- tmbstan(obj3, chains = 3, iter = 5000)

#save(m1,m2,m3,file="m1m2m3.rds")

out <- data.frame(
  type =
    c("Prior no adjustment", "Prior with adjustment")
)

out$mean_B0 <- c(mean((extract(m2)$B0)), mean((extract(m3)$B0)))
out$median_B0 <- c(median((extract(m2)$B0)), median((extract(m3)$B0)))
out$mean_kappa <- c(mean((exp(extract(m2)$ln_kappa))), mean((exp(extract(m3)$ln_kappa))))
out$median_kappa <- c(median((exp(extract(m2)$ln_kappa))),
                           median((exp(extract(m3)$ln_kappa))))
out$mean_range <- c(mean(sqrt(8)/(exp(extract(m2)$ln_kappa))),
                         mean(sqrt(8)/(exp(extract(m3)$ln_kappa))))
out$median_range <- c(median(sqrt(8)/(exp(extract(m2)$ln_kappa))),
                           median(sqrt(8)/(exp(extract(m3)$ln_kappa))))
out$mean_sigma <- c(mean((exp(extract(m2)$ln_sigma))),
                         mean((exp(extract(m3)$ln_sigma))))
out$median_sigma <- c(median((exp(extract(m2)$ln_sigma))),
                           median((exp(extract(m3)$ln_sigma))))
#sigma_O = 1/ [exp(ln_tau_O) * (sqrt(4 * pi) * kappa[1]]
out$mean_sigmaO <- c(mean(1/(exp(extract(m2)$ln_tau_O) * sqrt(4*pi) * exp(extract(m2)$ln_kappa))),
                      mean(1/(exp(extract(m3)$ln_tau_O) * sqrt(4*pi) * exp(extract(m3)$ln_kappa))))
out$median_sigmaO <- c(median(1/(exp(extract(m2)$ln_tau_O) * sqrt(4*pi) * exp(extract(m2)$ln_kappa))),
                          median(1/(exp(extract(m3)$ln_tau_O) * sqrt(4*pi) * exp(extract(m3)$ln_kappa))))

# This shows:
# B0 is off without the PC prior, estimated very well with prior, and perfectly with adjustment
# range parameter seems consistently off
# sigma obs is estimated well for all models
# sigma_O is biased high for all (0.27)
