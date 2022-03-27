library(sdmTMB)
pcod_spde <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 25)
m <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
  data = pcod_2011, time = "year", mesh = pcod_spde, family = tweedie(link = "log"))


nd <- matrix(NA, ncol = 30L, nrow = nrow(pcod_2011))
for (i in seq_len(ncol(nd))) {
  nd[,i] <- sdmTMB_simulate(previous_fit = m)$observed
}

s <- DHARMa::createDHARMa(
  simulatedResponse = nd,
  observedResponse = pcod_2011$density,
  fittedPredictedResponse = exp(predict(m, newdata = NULL)$est),
  method = "PIT"
)

plot(s)

plot(s, quantreg = FALSE, rank = TRUE)

r <- residuals(m, type = "randomized-quantile")
qqnorm(r);qqline(r)

set.seed(10292)
predictor_dat <- data.frame(
  X = runif(600), Y = runif(600),
  a1 = rnorm(600)
)
mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)
sim_dat <- sdmTMB_simulate(
  formula = ~ 1 + a1,
  data = predictor_dat,
  time = NULL,
  mesh = mesh,
  # family = Gamma(link = "log"),
  range = 0.5,
  sigma_E = NULL,
  phi = 0.1,
  sigma_O = 0.2,
  seed = 423,
  B = c(0.2, -0.4) # B0 = intercept, B1 = a1 slope
)
# fit <- sdmTMB(observed ~ a1, data = sim_dat, mesh = mesh, family = Gamma(link = "log"),)
fit <- sdmTMB(observed ~ a1, data = sim_dat, mesh = mesh)
# fit <- sdmTMB(observed ~ a1, data = sim_dat, mesh = mesh, family = nbinom2())
fit
r <- residuals(fit)
qqnorm(r);qqline(r)

mini_sim <- function(fit) {
  params <- fit$tmb_obj$env$parList()
  tmb_data <- fit$tmb_data
  # tmb_data$sim_re <- rep(1, length(tmb_data$sim_re))
  newobj <- TMB::MakeADFun(data = tmb_data, map = fit$tmb_map,
    random = fit$tmb_random, parameters = params, DLL = "sdmTMB")
  s <- newobj$simulate()
  s$y_i
}
nd <- matrix(NA, ncol = 20L, nrow = nrow(sim_dat))
for (i in seq_len(ncol(nd))) {
  nd[,i] <- mini_sim(fit)
}
s <- DHARMa::createDHARMa(
  simulatedResponse = nd,
  observedResponse = sim_dat$observed,
  fittedPredictedResponse = fit$family$linkinv(predict(fit, newdata = NULL)$est),
  method = "PIT"
)
plot(s)

# library(tmbstan)
# m <- tmbstan(fit$tmb_obj, iter = 200, chains = 1, warmup = 180)
# p <- predict(fit, tmbstan_model = m)
# p <- p + rnorm(nrow(p), 0, sd = mean(exp(extract(m)[['ln_phi']])))
# p <- as.matrix(p)
# s <- DHARMa::createDHARMa(
#   simulatedResponse = p,
#   observedResponse = sim_dat$observed,
#   fittedPredictedResponse = as.numeric(p[,5]),
#   method = "PIT"
# )
# plot(s)
# or do any PIT residuals randomization here...

# e.g.

pcod_spde <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 10)
fit <- sdmTMB(
  density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
  data = pcod_2011, time = "year", mesh = pcod_spde,
  family = tweedie(link = "log")
)
r <- residuals(fit)
qqnorm(r, main = "Naive QQ plot")
qqline(r)

stan_fit <- tmbstan::tmbstan(fit$tmb_obj, iter = 100, chains = 1, warmup = 98)
eta <- predict(fit, tmbstan_model = stan_fit)
mu <- fit$family$linkinv(as.numeric(eta[,1]))
r_mcmc <- sdmTMB:::qres_tweedie(fit, y = pcod_2011$density, mu = mu)
qqnorm(r_mcmc, main = "MCMC-based QQ plot")
qqline(r_mcmc)
