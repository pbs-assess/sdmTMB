set.seed(1)
library(MASS)
simdata <- function() {
  local(
    {
      rho <- 0.9
      sds <- seq(0.5, 2, length = stateDim)
      sdObs <- rep(1, stateDim)
      corrMat <- matrix(0.0, stateDim, stateDim)
      for (i in 1:stateDim) {
        for (j in 1:stateDim) {
          corrMat[i, j] <- rho^abs(i - j)
        }
      }
      Sigma <- corrMat * (sds %o% sds)
      d <- matrix(NA, timeSteps, stateDim)
      obs <- d
      ## init state
      d[1, ] <- rnorm(stateDim)
      i <- 1
      obs[i, ] <- d[i, ] + rnorm(stateDim, rep(0, stateDim), sdObs)
      for (i in 2:timeSteps) {
        d[i, ] <- d[i - 1, ] + mvrnorm(1, rep(0, stateDim), Sigma = Sigma)
        obs[i, ] <- d[i, ] + rnorm(stateDim, rep(0, stateDim), sdObs)
      }
      matplot(d, type = "l")
      matpoints(obs)
    },
    .GlobalEnv
  )
}
stateDim <- 3
timeSteps <- 100
simdata()
data <- list(obs = t(obs))
parameters <- list(
  u = data$obs * 0,
  transf_rho = 0.1,
  logsds = sds * 0,
  logsdObs = sdObs * 0
)

d <- data.frame(
  y = c(obs[, 1], obs[, 2], obs[, 3]),
  year = rep(1:nrow(obs), 3),
  group = c(
    rep("a", nrow(obs)),
    rep("b", nrow(obs)),
    rep("c", nrow(obs))
  )
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

nd <- expand.grid(year = 1:100, group = c("a", "b", "c"))
pred <- predict(fit, newdata = nd)

head(pred)

library(ggplot2)
ggplot(pred, aes(year, est, colour = group)) +
  geom_line() +
  geom_point(data = d, mapping = aes(y = y))
