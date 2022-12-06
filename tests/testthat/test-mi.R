test_that("Metabolic index threshold models fit", {
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
  dat <- data.frame(y = obs, po2 = po2, invtemp = invtemp)

  start <- matrix(0, ncol = 1, nrow = 4)
  start[1,1] <- x50
  start[2,1] <- delta
  start[3,1] <- beta3
  start[4,1] <- Eo

  m2 <- sdmTMB(y ~ logistic(mi), data = dat, spatial = "off",
    family = lognormal(),
    control = sdmTMBcontrol(start = list(b_threshold = start)))
  tidy(m2)
  print(m2)

  b <- tidy(m2)
  bv <- tidy(m2, "ran_pars")
  expect_equal(b$estimate[b$term == "(Intercept)"], beta1, tolerance = 0.05)
  expect_equal(b$estimate[b$term == "mi-delta"], delta, tolerance = 0.1)
  expect_equal(b$estimate[b$term == "mi-s50"], x50, tolerance = 0.05)
  expect_equal(b$estimate[b$term == "mi-smax"], beta3, tolerance = 0.1)
  expect_equal(bv$estimate[bv$term == "phi"], sigma, tolerance = 0.05)

  p <- predict(m2, newdata = dat)
  expect_true(sum(is.na(p$est)) == 0L)

  plot(log_mu, p$est)
  expect_gt(cor(log_mu, p$est), 0.9)
  r <- residuals(m2)
  qqnorm(r);qqline(r)

  # with priors
  m3 <- sdmTMB(y ~ logistic(mi), data = dat, spatial = "off",
    family = lognormal(),
    priors = sdmTMBpriors(threshold = normal(start[,1], c(0.5, 0.5, 0.5, 0.5))),
    control = sdmTMBcontrol(start = list(b_threshold = start)))

  expect_equal(m2$model$par, m3$model$par, tolerance = 0.1)
  expect_false(identical(m2$model$par, m3$model$par))

  expect_error({
    m4 <- sdmTMB(y ~ logistic(mi), data = dat, spatial = "off",
      family = lognormal(),
      priors = sdmTMBpriors(threshold = normal(start[1:3,1], c(0.5, 0.5, 0.5))),
      control = sdmTMBcontrol(start = list(b_threshold = start)))
  }, regexp = "number of")
})

test_that("MI simulation works", {
  set.seed(123)
  # make fake predictor(s) (a1) and sampling locations:
  predictor_dat <- data.frame(
    X = runif(300), Y = runif(300), a1 = rnorm(300)
  )
  mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)

  beta1 <- -0.4
  beta3 <- 0.3
  delta <- 3
  Eo <- 0.1
  x50 <- 5

  N <- 300
  set.seed(123)
  invtemp <- rnorm(N)
  po2 <- rlnorm(N)
  mi <- po2 * exp(Eo * invtemp)
  log_mu <- beta1 + beta3 * (1 / (1 + exp(-log(19) * (mi - x50) / delta)) - 1)
  mu <- exp(log_mu)

  predictor_dat$po2 <- po2
  predictor_dat$invtemp <- invtemp

  tb <- numeric(4L)
  tb[1] <- x50
  tb[2] <- delta
  tb[3] <- beta3
  tb[4] <- Eo

  sim_dat <- sdmTMB_simulate(
    formula = ~ 1 + logistic(mi),
    data = predictor_dat,
    mesh = mesh,
    family = lognormal(),
    range = 0.5,
    phi = 0.1,
    sigma_O = 0,
    seed = 42,
    threshold_coefs = tb,
    B = beta1
  )

  if (FALSE) {
    plot(sim_dat$mu, mu)
    abline(0, 1)
  }

  expect_equal(sim_dat$mu, mu, tolerance = 0.0001)
})
