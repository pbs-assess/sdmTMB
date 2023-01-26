library(sdmTMB)
mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
fit <- sdmTMB(
  present ~ as.factor(year),
  data = pcod_2011, mesh = mesh,
  family = binomial()
)
stan_args <- list(control = list(adapt_delta = 0.8, max_treedepth = 10))
r3 <- residuals(
  fit, type = "mle-mcmc", mcmc_iter = 11, mcmc_warmup = 10,
  stan_args = stan_args
)
