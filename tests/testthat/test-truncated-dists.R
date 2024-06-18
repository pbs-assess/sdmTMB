test_that("Delta-truncated NB2 works with various types of predictions", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")

  set.seed(1)
  d0 <- data.frame(y = rnbinom(n = 1e3, mu = 3, size = 5))

  fit <- sdmTMB(
    y ~ 1,
    data = d0,
    spatial = FALSE,
    family = delta_truncated_nbinom2()
  )
  fit$sd_report

  fit2 <- glmmTMB::glmmTMB(
    y ~ 1,
    ziformula = ~ 1,
    data = d0,
    family = glmmTMB::truncated_nbinom2()
  )

  p <- predict(fit)
  p2 <- predict(fit2)
  expect_equal(p2, p$est2, tolerance = 1e-4)

  p <- predict(fit, type = "response")
  p2 <- predict(fit2, type = "response")
  expect_equal(p2, p$est, tolerance = 1e-4)

  p <- predict(fit)
  p2 <- predict(fit2, type = "zlink")
  expect_equal(-p2, p$est1, tolerance = 1e-4)

  p <- predict(fit, type = "response")
  p2 <- predict(fit2, type = "conditional")
  expect_equal(p2, p$est2, tolerance = 1e-4)
})

test_that("Delta-truncated NB1 works with various types of predictions", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")

  set.seed(1)
  d0 <- data.frame(y = rnbinom(n = 1e3, mu = 3, size = 5))

  fit <- sdmTMB(
    y ~ 1,
    data = d0,
    spatial = FALSE,
    family = delta_truncated_nbinom1()
  )
  fit$sd_report

  fit2 <- glmmTMB::glmmTMB(
    y ~ 1,
    ziformula = ~ 1,
    data = d0,
    family = glmmTMB::truncated_nbinom1()
  )

  p <- predict(fit)
  p2 <- predict(fit2)
  expect_equal(p2, p$est2, tolerance = 1e-4)

  p <- predict(fit, type = "response")
  p2 <- predict(fit2, type = "response")
  expect_equal(p2, p$est, tolerance = 1e-4)

  p <- predict(fit)
  p2 <- predict(fit2, type = "zlink")
  expect_equal(-p2, p$est1, tolerance = 1e-4)

  p <- predict(fit, type = "response")
  p2 <- predict(fit2, type = "conditional")
  expect_equal(p2, p$est2, tolerance = 1e-4)
})

test_that("Truncated NB2 works with various types of predictions", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  set.seed(1)
  d0 <- data.frame(y = rnbinom(n = 1e3, mu = 3, size = 5))
  d0 <- subset(d0, y > 0)
  fit <- sdmTMB(
    y ~ 1,
    data = d0,
    spatial = FALSE,
    family = truncated_nbinom2()
  )
  fit2 <- glmmTMB::glmmTMB(
    y ~ 1,
    data = d0,
    family = glmmTMB::truncated_nbinom2()
  )
  p <- predict(fit)
  p2 <- predict(fit2)
  expect_equal(p2, p$est, tolerance = 1e-4)
  p <- predict(fit, type = "response")
  p2 <- predict(fit2, type = "response")
  expect_equal(p2, p$est, tolerance = 1e-4)
})

test_that("Truncated NB1 works with various types of predictions", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  set.seed(1)
  d0 <- data.frame(y = rnbinom(n = 1e3, mu = 3, size = 5))
  d0 <- subset(d0, y > 0)
  fit <- sdmTMB(
    y ~ 1,
    data = d0,
    spatial = FALSE,
    family = truncated_nbinom1()
  )
  fit2 <- glmmTMB::glmmTMB(
    y ~ 1,
    data = d0,
    family = glmmTMB::truncated_nbinom1()
  )
  p <- predict(fit)
  p2 <- predict(fit2)
  expect_equal(p2, p$est, tolerance = 1e-4)
  p <- predict(fit, type = "response")
  p2 <- predict(fit2, type = "response")
  expect_equal(p2, p$est, tolerance = 1e-4)
})
