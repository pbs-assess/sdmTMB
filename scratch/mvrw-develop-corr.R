library(ggplot2)
library(sdmTMB)

# as mvrw-develp but with unstructured correlation matrix rather than ordinal

set.seed(12928)
stateDim <- 7
timeSteps <- 400
sdObs <- 0.3
A <- matrix(runif(stateDim^2) * 2 - 1, ncol = stateDim)
Sigma <- t(A) %*% A
par(mfrow = c(1, 1))
image(Sigma)
print(Sigma)
sqrt(diag(Sigma))
d <- matrix(NA, timeSteps, stateDim)
obs <- d
d[1, ] <- MASS::mvrnorm(1, rep(0, stateDim), Sigma = Sigma) # initial state
i <- 1
obs[i, ] <- d[i, ] + rnorm(stateDim, rep(0, stateDim), sdObs)
for (i in 2:timeSteps) {
  d[i, ] <- d[i - 1, ] + MASS::mvrnorm(1, rep(0, stateDim), Sigma = Sigma)
  obs[i, ] <- d[i, ] + rnorm(stateDim, rep(0, stateDim), sdObs)
}
matplot(d, type = "l")
truth <- d
matpoints(obs, pch = 21)

d <- data.frame(
  y = reshape2::melt(obs)[, 3],
  year = rep(1:nrow(obs), stateDim),
  group = rep(letters[1:stateDim], each = timeSteps)
)

fit <- sdmTMB(
  y ~ 1,
  data = d,
  spatial = "off",
  time = "year",
  spatiotemporal = "off",
  mvrw_category = "group",
  control = sdmTMBcontrol(
    map = list(b_j = factor(NA)),
    start = list(b_j = 0)
  )
)

fit_null <- sdmTMB(
  y ~ 1,
  data = d,
  spatial = "off",
  time = "year",
  spatiotemporal = "off"
)

# fit$tmb_params$mvrw_logsds
# fit$tmb_params$mvrw_u
# fit$tmb_params$mvrw_rho

fit$sd_report
r <- fit$tmb_obj$report(fit$tmb_obj$env$last.par.best)
r$mvrw_Sigma
expect_equal(nrow(r$mvrw_Sigma), stateDim)
expect_equal(ncol(r$mvrw_Sigma), stateDim)
expect_equal(diag(r$mvrw_Sigma), rep(1, stateDim))

pars <- fit$tmb_obj$env$parList()

pars$mvrw_rho
pars$mvrw_logsds
pars$mvrw_u[1:3, 1:3]

# https://stats.stackexchange.com/a/407954
# efficient but impossible to read!
cor2cov <- function(R, S) {
  sweep(sweep(R, 1, S, "*"), 2, S, "*")
}
Sigma_hat <- cor2cov(r$mvrw_Sigma, exp(pars$mvrw_logsds))
true_sds <- sqrt(diag(Sigma))
par(mfrow = c(1, 2))
plot(true_sds, exp(pars$mvrw_logsds), main = "SDs")
abline(0, 1)
plot(as.numeric(Sigma_hat), as.numeric(Sigma), main = "Covariances")
abline(0, 1)

expect_equal(exp(pars$ln_phi), sdObs, tolerance = 0.025)
expect_gt(cor(as.numeric(Sigma_hat), as.numeric(Sigma)), 0.98)
expect_gt(cor(true_sds, exp(pars$mvrw_logsds)), 0.98)


