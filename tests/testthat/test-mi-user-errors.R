test_that("Missing input data when using metabolic index model", {
  beta1 <- -0.4
  beta3 <- 0.3
  delta <- 3
  Eo <- 0.1
  x50 <- 5
  
  N <- 3000
  set.seed(123)
  invtemp <- rnorm(N)
  po2 <- rlnorm(N)
  mi <- po2 * exp(Eo * invtemp)
  log_mu <- beta1 + beta3 * (1 / (1 + exp(-log(19) * (mi - x50) / delta)) - 1)
  mu <- exp(log_mu)
  
  # plot(mu)
  # plot(mi, log_mu)
  # plot(mi, mu)
  sigma <- 0.05
  obs <- rlnorm(N, log_mu - 0.5 * sigma^2, sigma)
  if (FALSE) {
    plot(mi, log(mu))
    plot(mi, log(obs))
  }
#Dataframe missing po2
dat <- data.frame(y = obs, invtemp = invtemp)

#Dataframe missing invtemp
dat2 <- data.frame(y = obs, po2 = po2)

#Dataframe includes mi
dat3 <- data.frame(y = obs, po2 = po2, invtemp=invtemp, mi=mi)

#Starting parameters
start <- matrix(0, ncol = 1, nrow = 4)
start[1,1] <- x50
start[2,1] <- delta
start[3,1] <- beta3
start[4,1] <- Eo

expect_error({
m2 <- sdmTMB(y ~ logistic(mi), data = dat, spatial = "off",
             family = lognormal(),
             control = sdmTMBcontrol(start = list(b_threshold = start)))
})

expect_error({
  m3 <- sdmTMB(y ~ logistic(mi), data = dat2, spatial = "off",
               family = lognormal(),
               control = sdmTMBcontrol(start = list(b_threshold = start)))
})

expect_error({
    m4 <- sdmTMB(y ~ logistic(mi), data = dat3, spatial = "off",
                 family = lognormal(),
                 control = sdmTMBcontrol(start = list(b_threshold = start)))
})
})
