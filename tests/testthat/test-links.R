test_that("cloglog works", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")

  m1 <- sdmTMB(present ~ 1, family = binomial(link = "cloglog"),
    spatial = "off", mesh = pcod_mesh_2011, data = pcod_2011)

  m2 <- stats::glm(present ~ 1, family = binomial(link = "cloglog"),
    data = pcod_2011)

  b <- tidy(m1)
  expect_equal(b$estimate, coef(m2)[[1]], tolerance = 1e-5)

  # wrong link:
  m3 <- stats::glm(present ~ 1, family = binomial(link = "logit"),
    data = pcod_2011)
  expect_true(abs(b$estimate - coef(m3)[[1]]) > 0.001) # should not match!

  set.seed(1)
  s <- simulate(m1, nsim = 100L)
  expect_equal(unique(as.numeric(s)), c(0, 1))

  r <- residuals(m1)
  qqnorm(r)

  # r <- residuals(m1, type = "mle-mcmc", mcmc_iter = 101, mcmc_warmup = 100)
  # r <- residuals(m1, type = "mle-mcmc")
  # r <- residuals(m3, type = "mle-mcmc")
})
