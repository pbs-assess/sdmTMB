# Simulation proof of when you do/do not want the Jacobian adjustment
# when placing a prior on a transformed parameter
library(sdmTMB)
library(TMB)
library(tmbstan)

# simulate spatial data from a single year
set.seed(123)

compile("inst/pc_matern_prior_jacobian2.cpp")
dyn.load(dynlib("inst/pc_matern_prior_jacobian2"))

data <- list(
  matern_range = 5,
  range_prob = 0.05,
  matern_SD = 1,
  SD_prob = 0.05,
  include_jacob = 0L
)

pars <- list(
  ln_tau_O = 0,
  ln_kappa = 0
)

obj <- MakeADFun(data, pars, DLL = "pc_matern_prior_jacobian2")

range_expected <- 10
kappa_expected = sqrt(8) / range_expected

sigma_O_expected <- 0.5
ln_tau_O_expected <- log(1 / (sqrt(4 * pi) * kappa_expected * sigma_O_expected))

log(kappa_expected)
ln_tau_O_expected

pars <- list(
  ln_tau_O = ln_tau_O_expected,
  ln_kappa = log(kappa_expected)
)
obj$fn(pars)

# inla_docs_log <- function(range = 1, sigma = 0.6, matern_range = 5,
#   range_prob = 0.05, matern_SD = 2, SD_prob = 0.05) {
#   d <- 2
#   dhalf <- d / 2
#   lam1 <- -log(range_prob) * matern_range^dhalf
#   lam2 <- -log(SD_prob) / matern_SD
#   log(dhalf) + log(lam1) + log(range^(-1 - dhalf)) - lam1 * range^-dhalf +
#     log(lam2) - lam2 * sigma
# }

ranges <- seq(0.1, 10, length.out = 100)
out <- sapply(ranges, function(.r)
  sdmTMB:::inla_docs_log(.r, sigma = 0.6, matern_range = 5, range_prob = 0.05, matern_SD = 2, SD_prob = 0.05))
plot(ranges, out)

# -------------------

convert_range_sigma <- function(range, sigma) {
  kappa = sqrt(8) / range
  ln_tau_O <- log(1 / (sqrt(4 * pi) * kappa * sigma))
  list(
    ln_kappa = log(kappa_expected),
    ln_tau_O = ln_tau_O
  )
}

data <- list(
  matern_range = 5,
  range_prob = 0.05,
  matern_SD = 2,
  SD_prob = 0.05,
  include_jacob = 0L
)

sigma_O <- 0.6
out2 <- sapply(ranges, function(.r) {
  pars <- list(
    ln_tau_O = convert_range_sigma(.r, sigma_O)$ln_tau_O,
    ln_kappa = convert_range_sigma(.r, sigma_O)$ln_kappa
  )
  obj <- MakeADFun(data, pars, DLL = "pc_matern_prior_jacobian2")
  -obj$fn(pars)
})
plot(ranges, out2)

plot(out, out2)

# -------------------





fit <- nlminb(obj$par, obj$fn, obj$gr)
fit$par

obj$fn(fit$par)

fit_stan <- tmbstan(obj, iter = 100, chains = 1, control = list(max_treedepth = 15))

# ...


# ---------------------




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


# sample with Stan -- first model (no prior) throws 200+ divergent transitions
options(mc.cores = parallel::detectCores())
m1 <- tmbstan(obj1, chains = 3, iter = 5000)
m2 <- tmbstan(obj2, chains = 3, iter = 5000)
m3 <- tmbstan(obj3, chains = 3, iter = 5000)

#save(m1,m2,m3,file="m1m2m3.rds")

out <- data.frame(
  type =
    c("No prior", "Prior no adjustment", "Prior with adjustment")
)

out$mean_B0 <- c(mean((extract(m1)$B0)), mean((extract(m2)$B0)),
                      mean((extract(m3)$B0)))
out$median_B0 <- c(median((extract(m1)$B0)), median((extract(m2)$B0)),
                      median((extract(m3)$B0)))
out$mean_kappa <- c(mean((exp(extract(m1)$ln_kappa))), mean((exp(extract(m2)$ln_kappa))),
                         mean((exp(extract(m3)$ln_kappa))))
out$median_kappa <- c(median((exp(extract(m1)$ln_kappa))), median((exp(extract(m2)$ln_kappa))),
                           median((exp(extract(m3)$ln_kappa))))
out$mean_range <- c(mean(sqrt(8)/(exp(extract(m1)$ln_kappa))), mean(sqrt(8)/(exp(extract(m2)$ln_kappa))),
                         mean(sqrt(8)/(exp(extract(m3)$ln_kappa))))
out$median_range <- c(median(sqrt(8)/(exp(extract(m1)$ln_kappa))), median(sqrt(8)/(exp(extract(m2)$ln_kappa))),
                           median(sqrt(8)/(exp(extract(m3)$ln_kappa))))
out$mean_sigma <- c(mean((exp(extract(m1)$ln_sigma))), mean((exp(extract(m2)$ln_sigma))),
                         mean((exp(extract(m3)$ln_sigma))))
out$median_sigma <- c(median((exp(extract(m1)$ln_sigma))), median((exp(extract(m2)$ln_sigma))),
                           median((exp(extract(m3)$ln_sigma))))
#sigma_O = 1/ [exp(ln_tau_O) * (sqrt(4 * pi) * kappa[1]]
out$mean_sigmaO <- c(mean(1/(exp(extract(m1)$ln_tau_O) * sqrt(4*pi) * exp(extract(m1)$ln_kappa))),
                          mean(1/(exp(extract(m2)$ln_tau_O) * sqrt(4*pi) * exp(extract(m2)$ln_kappa))),
                      mean(1/(exp(extract(m3)$ln_tau_O) * sqrt(4*pi) * exp(extract(m3)$ln_kappa))))
out$median_sigmaO <- c(median(1/(exp(extract(m1)$ln_tau_O) * sqrt(4*pi) * exp(extract(m1)$ln_kappa))),
                          median(1/(exp(extract(m2)$ln_tau_O) * sqrt(4*pi) * exp(extract(m2)$ln_kappa))),
                          median(1/(exp(extract(m3)$ln_tau_O) * sqrt(4*pi) * exp(extract(m3)$ln_kappa))))

# This shows:
# B0 is off without the PC prior, estimated very well with prior, and perfectly with adjustment
# range parameter seems consistently off
# sigma obs is estimated well for all models
# sigma_O is biased high for all (0.27)
