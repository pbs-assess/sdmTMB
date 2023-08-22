library(ggplot2)
library(sdmTMB)

set.seed(12928)
rho <- 0.6
stateDim <- 10
timeSteps <- 200
sds <- rlnorm(stateDim, meanlog = log(0.4), sdlog = 0.4)
sdObs <- rep(0.2, stateDim)
corrMat <- matrix(0.0, stateDim, stateDim)
for (i in 1:stateDim) {
  for (j in 1:stateDim) {
    corrMat[i, j] <- rho^abs(i - j)
  }
}
Sigma <- corrMat * (sds %o% sds)
d <- matrix(NA, timeSteps, stateDim)
obs <- d
d[1, ] <- rnorm(stateDim, 0, 1) # initial state
i <- 1
obs[i, ] <- d[i, ] + rnorm(stateDim, rep(0, stateDim), sdObs)
for (i in 2:timeSteps) {
  d[i, ] <- d[i - 1, ] + MASS::mvrnorm(1, rep(0, stateDim), Sigma = Sigma)
  obs[i, ] <- d[i, ] + rnorm(stateDim, rep(0, stateDim), sdObs)
}
matplot(d, type = "l")
truth <- d
matpoints(obs)

d <- data.frame(
  y = reshape2::melt(obs)[,3],
  year = rep(1:nrow(obs), stateDim),
  group = rep(letters[1:stateDim], each = timeSteps)
)
head(d)

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
fit$sd_report

fit$tmb_data$year_i
fit$tmb_random

fit$tmb_params$mvrw_logsds
fit$tmb_params$mvrw_phi
fit$tmb_params$mvrw_u

p <- fit$tmb_obj$env$parList()

matplot(obs, type = "l", lwd = 1, lty = 2)
matpoints(t(p$mvrw_u), pch = 21)

nd <- expand.grid(year = unique(d$year), group = c("a", "b", "c"))
expect_error(predict(fit, newdata = nd), regexp = "missing")

nd <- expand.grid(year = unique(d$year), group = c(unique(d$group), "zz"))
expect_error(predict(fit, newdata = nd), regexp = "extra")

nd <- expand.grid(year = unique(d$year), group = unique(d$group))
pred <- predict(fit, newdata = nd)

head(pred)

ggplot(pred, aes(year, est, colour = group)) +
  geom_line() +
  geom_point(data = d, mapping = aes(y = y))

2 * plogis(p$mvrw_phi) - 1
head(t(p$mvrw_u))
s <- as.list(fit$sd_report, "Std. Error")
lwr <- 2 * plogis(p$mvrw_phi - 2 * s$mvrw_phi) - 1
upr <- 2 * plogis(p$mvrw_phi + 2 * s$mvrw_phi) - 1
expect_lt(lwr, rho)
expect_gt(upr, rho)

expect_gt(cor(pred$est, d$y), 0.99)
for (i in 1:stateDim) {
  expect_gt(cor(p$mvrw_u[i,], truth[,i]), 0.98)
}

plot(log(sds), p$mvrw_logsds)
abline(0, 1)
expect_gt(cor(log(sds), p$mvrw_logsds), 0.98)
