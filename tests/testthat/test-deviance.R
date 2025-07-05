test_that("Deviance calculations are correct", {
  # Poisson
  set.seed(123)
  x <- rnorm(10)
  y <- rpois(10, exp(x))
  d <- data.frame(x = x, y = y)
  m <- sdmTMB(y ~ 1 + x, family = poisson(), spatial = "off", data = d)
  mglm <- glm(y ~ 1 + x, family = poisson(), data = d)
  r1 <- unname(residuals(mglm))
  r2 <- residuals(m, type = "deviance")
  expect_equal(r1, r2, tolerance = 0.01)
  expect_equal(deviance(mglm), deviance(m), tolerance = 0.001)

  # Gamma
  set.seed(123)
  x <- rnorm(10)
  y <- rgamma(10, shape = 0.4, scale = exp(x) / 0.4)
  d <- data.frame(x = x, y = y)
  m <- sdmTMB(y ~ x, family = Gamma(link = "log"), spatial = "off", data = d)
  mglm <- glm(y ~ x, family = Gamma(link = "log"), data = d)
  r1 <- unname(residuals(mglm))
  r2 <- residuals(m, type = "deviance")
  expect_equal(r1, r2, tolerance = 0.01)
  expect_equal(deviance(mglm), deviance(m), tolerance = 0.001)

  # Binomial
  set.seed(123)
  x <- rnorm(10)
  y <- rbinom(10, size = 1, prob = plogis(x))
  d <- data.frame(x = x, y = y)
  m <- sdmTMB(y ~ x, family = binomial(link = "logit"), spatial = "off", data = d)
  mglm <- glm(y ~ x, family = binomial(link = "logit"), data = d)
  r1 <- unname(residuals(mglm))
  r2 <- residuals(m, type = "deviance")
  expect_equal(r1, r2, tolerance = 0.01)
  expect_equal(deviance(mglm), deviance(m), tolerance = 0.001)

  # Binomial, size > 1
  set.seed(123)
  x <- rnorm(10)
  w <- sample(1:5, size = 10, replace = TRUE)
  y <- rbinom(10, size = w, prob = plogis(x))
  d <- data.frame(x = x, y = y, prop = y / w)
  m <- sdmTMB(prop ~ x, family = binomial(link = "logit"), weights = w, spatial = "off", data = d)
  expect_error(residuals(m, type = "deviance"), regexp = "size")
  expect_error(deviance(m), regexp = "size")

  # # NB2
  # set.seed(123)
  # x <- rnorm(20)
  # y <- rnbinom(20, size = 0.4, mu = exp(x))
  # d <- data.frame(x = x, y = y)
  # m <- sdmTMB(y ~ x, family = nbinom2(), spatial = "off", data = d)
  # mglm <- MASS::glm.nb(y ~ x, data = d)
  # r1 <- unname(residuals(mglm))
  # r2 <- residuals(m, type = "deviance")
  # expect_equal(r1, r2, tolerance = 0.01)
  # expect_equal(deviance(mglm), deviance(m), tolerance = 0.01)

  # NB1
  set.seed(1)
  x <- rnorm(20)
  y <- rnbinom1(20, phi = 2, mu = exp(x))
  d <- data.frame(x = x, y = y)
  m <- sdmTMB(y ~ x, family = nbinom1(), spatial = "off", data = d)
  mglm <- glmmTMB::glmmTMB(y ~ x, family = glmmTMB::nbinom1(), data = d)
  r1 <- unname(residuals(mglm, type = "deviance"))
  r2 <- residuals(m, type = "deviance")
  expect_equal(r1, r2, tolerance = 0.01)
  expect_equal(deviance(mglm), deviance(m), tolerance = 0.01)

  # Tweedie
  set.seed(1)
  x <- rnorm(20)
  y <- mgcv::rTweedie(exp(x), 1.5, 1.2)
  d <- data.frame(x = x, y = y)
  m <- sdmTMB(y ~ x, family = tweedie(), spatial = "off", data = d)
  library(mgcv)
  mglm <- mgcv::gam(y ~ x, family = mgcv::tw(), data = d)
  r1 <- unname(residuals(mglm, type = "deviance"))
  r2 <- residuals(m, type = "deviance")
  expect_equal(r1, r2, tolerance = 0.1)
  expect_equal(deviance(mglm), deviance(m), tolerance = 0.01)

  # Gaussian
  set.seed(1)
  x <- rnorm(20)
  y <- rnorm(20, x * 0.3, sd = 0.6)
  d <- data.frame(x = x, y = y)
  m <- sdmTMB(y ~ x, spatial = "off", data = d)
  mglm <- glm(y ~ x, data = d)
  r1 <- unname(residuals(mglm, type = "deviance"))
  r2 <- residuals(m, type = "deviance")
  expect_equal(r1, r2, tolerance = 0.01)
  expect_equal(deviance(mglm), deviance(m), tolerance = 0.01)

  # # lognormal
  # set.seed(1)
  # x <- rnorm(20)
  # y <- rlnorm(20, x * 0.3, sdlog = 0.2)
  # d <- data.frame(x = x, y = y)
  # m <- sdmTMB(y ~ x, family = lognormal(), spatial = "off", data = d)
  # mglm <- tinyVAST::tinyVAST(y ~ x, data = d, family = lognormal())
  # r1 <- unname(residuals(mglm, type = "deviance"))
  # r2 <- residuals(m, type = "deviance")
  # expect_equal(r1, r2, tolerance = 0.01)
  # expect_equal(sum(r1^2), deviance(m), tolerance = 0.01)

  # # delta gamma
  # m <- sdmTMB(density ~ depth_scaled, family = sdmTMB::delta_gamma(), spatial = "off", data = pcod)
  # dev1 <- deviance(m)
  # m2 <- tinyVAST::tinyVAST(density ~ depth_scaled,
  #   family = tinyVAST::delta_gamma(),
  #   data = as.data.frame(pcod), delta_options = list(formula = ~depth_scaled)
  # )
  # r2 <- m2$obj$report(m2$obj$env$last.par.best)
  # dev2 <- r2$deviance
  # expect_equal(dev1, dev2)

  # # delta lognormal
  # m <- sdmTMB(density ~ depth_scaled, family = sdmTMB::delta_lognormal(), spatial = "off", data = pcod)
  # dev1 <- deviance(m)
  # m2 <- tinyVAST::tinyVAST(density ~ depth_scaled,
  #   family = tinyVAST::delta_lognormal(),
  #   data = as.data.frame(pcod), delta_options = list(formula = ~depth_scaled)
  # )
  # r2 <- m2$obj$report(m2$obj$env$last.par.best)
  # dev2 <- r2$deviance
  # expect_equal(dev1, dev2)

  # # delta gamma poisson link
  # mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 8)
  # m0 <- sdmTMB(
  #   density ~ depth_scaled, mesh = mesh,
  #   family = delta_gamma(type = "poisson-link"), spatial = "on", data = pcod,
  # )
  # dev1 <- deviance(m0)
  # m2 <- tinyVAST::tinyVAST(
  #   density ~ depth_scaled,
  #   family = tinyVAST::delta_gamma(type = "poisson-link"), data = as.data.frame(pcod),
  #   delta_options = list(formula = ~depth_scaled)
  # )
  # r2 <- m2$obj$report()
  # dev2 <- r2$deviance
  # expect_equal(dev1, dev2)

  # # deviance explained
  # mnull <- sdmTMB(
  #   density ~ 1, spatial = "off", mesh = mesh,
  #   family = delta_gamma(type = "poisson-link"), data = pcod,
  # )
  # (devexplained <- 1 - deviance(m0) / deviance(mnull))
  #
  # expect_equal(devexplained, m2$deviance_explained, tolerance = 0.0001)
})
