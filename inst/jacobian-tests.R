# Simulation proof of when you do/do not want the Jacobian adjustment
# when placing a prior on a transformed parameter
# SA, 2021-06-11

library(TMB)
library(tmbstan)
options(mc.cores = 1L)

# no prior
compile("inst/lm.cpp")
dyn.load(dynlib("inst/lm"))
set.seed(123)
N <- 100L
data <- list(
  Y = rnorm(N, sd = 1) + seq(0, 1, length.out = N),
  x = seq(0, 1, length.out = N))
parameters <- list(a = 0, b = 1, logSigma = log(1))
obj <- MakeADFun(data, parameters, DLL = "lm", hessian = TRUE)
opt <- nlminb(obj$par, objective = obj$fn, gradient = obj$gr)
opt
exp(opt$par[[3]])

# prior on sigma; no jacobian adjustment
compile("inst/lm_prior_nojb.cpp")
dyn.load(dynlib("inst/lm_prior_nojb"))
obj2 <- MakeADFun(data, parameters, DLL = "lm_prior_nojb", hessian = TRUE)
opt2 <- nlminb(obj2$par, objective = obj2$fn, gradient = obj2$gr)
opt2
exp(opt2$par[[3]])

# prior on sigma; *with* jacobian adjustment
compile("inst/lm_prior_jb.cpp")
dyn.load(dynlib("inst/lm_prior_jb"))
obj3 <- MakeADFun(data, parameters, DLL = "lm_prior_jb", hessian = TRUE)
opt3 <- nlminb(obj3$par, objective = obj3$fn, gradient = obj3$gr)
opt3
exp(opt3$par[[3]])

# prior on sigma; no exp() transformation; use limits instead
# this version should always be correct
compile("inst/lm_limits.cpp")
dyn.load(dynlib("inst/lm_limits"))
parameters4 <- list(a = 0, b = 1, sigma = 1)
obj4 <- MakeADFun(data, parameters4, DLL = "lm_limits", hessian = TRUE)
opt4 <- nlminb(obj4$par,
  objective = obj4$fn, gradient = obj4$gr,
  lower = c(-Inf, -Inf, 0), upper = c(Inf, Inf, Inf)
)
opt4
opt4$par[[3]]

# sample with Stan
m1 <- tmbstan(obj, chains = 4, iter = 8000)
m2 <- tmbstan(obj2, chains = 4, iter = 8000)
m3 <- tmbstan(obj3, chains = 4, iter = 8000)
m4 <- tmbstan(obj4,
  chains = 4, iter = 8000,
  lower = c(-Inf, -Inf, 0), upper = c(Inf, Inf, Inf)
)

# summarize
out <- data.frame(
  type =
    c("No prior", "Prior no adjustment", "Prior with adjustment", "Prior limits")
)
out$stan_mean <- NA
out$stan_median <- NA
out$nlminb <- NA

out$stan_mean[1] <- mean(exp(extract(m1)$logSigma))
out$stan_median[1] <- mean(exp(extract(m1)$logSigma))
out$stan_mean[2] <- mean(exp(extract(m2)$logSigma))
out$stan_median[2] <- median(exp(extract(m2)$logSigma))
out$stan_mean[3] <- mean(exp(extract(m3)$logSigma))
out$stan_median[3] <- median(exp(extract(m3)$logSigma))
out$stan_mean[4] <- mean(extract(m4)$sigma)
out$stan_median[4] <- median(extract(m4)$sigma)

out$nlminb[1] <- exp(opt$par[[3]])
out$nlminb[2] <- exp(opt2$par[[3]])
out$nlminb[3] <- exp(opt3$par[[3]])
out$nlminb[4] <- opt4$par[[3]]

tibble::as_tibble(out)

# Conclusion:
# - you don't want the Jacobian adjustment for MAP/MLE

# What about simulating from the MVN covariance matrix?

sim_N <- 4e5
sd1 <- sdreport(obj)
sims1 <- mvtnorm::rmvnorm(sim_N, mean = sd1$par.fixed, sigma = sd1$cov.fixed)
sd2 <- sdreport(obj2)
sims2 <- mvtnorm::rmvnorm(sim_N, mean = sd2$par.fixed, sigma = sd2$cov.fixed)
sd3 <- sdreport(obj3)
sims3 <- mvtnorm::rmvnorm(sim_N, mean = sd3$par.fixed, sigma = sd3$cov.fixed)
sd4 <- sdreport(obj4)
sims4 <- mvtnorm::rmvnorm(sim_N, mean = sd4$par.fixed, sigma = sd4$cov.fixed)

out$sims_mean <- NA
out$sims_median <- NA

out$sims_mean[1] <- mean(exp(sims1[, 3]))
out$sims_median[1] <- median(exp(sims1[, 3]))

out$sims_mean[2] <- mean(exp(sims2[, 3]))
out$sims_median[2] <- median(exp(sims2[, 3]))

out$sims_mean[3] <- mean(exp(sims3[, 3]))
out$sims_median[3] <- median(exp(sims3[, 3]))

out$sims_mean[4] <- mean(I(sims4[, 3]))
out$sims_median[4] <- median(I(sims4[, 3]))

print(tibble::as_tibble(out))

# Conclusion:
# - sim from MVN is same as nlminb optimization
# - you only want Jacobian when invoking MCMC from Stan
