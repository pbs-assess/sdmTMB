# test_that("Metabolic index threshold models fit" {
#   beta1 <- -0.4
#   beta3 <- 0.3
#   delta <- 0.5
#   Eo <- 0.1
#   x50 <- 1
#
#   N <- 200
#   set.seed(1)
#   invtemp <- rnorm(N)
#   po2 <- rlnorm(N)
#   mi <- po2 * exp(Eo * invtemp)
#   log_mu <- beta1 + beta3 * (1 / (1 + exp(-log(19) * (mi - x50) / delta)) - 1)
#   mu <- exp(log_mu)
#
#   # plot(mu)
#   # plot(mi, log_mu)
#   # plot(mi, mu)
#   obs <- rlnorm(N, log_mu, 0.1)
#   plot(mi, log(obs))
#
#   dat <- data.frame(y = obs, po2 = po2, invtemp = invtemp)
#   m <- sdmTMB(y ~ logistic(mi), data = dat, spatial = "off")
# })
