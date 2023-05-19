test_that("Gamma, NB2, and lognormal mixtures fit and recover mixing proportion", {
  skip_on_ci()
  skip_on_cran()
  d <- pcod[pcod$year == 2017 & pcod$density > 0, ]
  spde <- make_mesh(d, c("X", "Y"), cutoff = 10)
  m <- sdmTMB(
    data = d, formula = density ~ 1,
    mesh = spde, family = gamma_mix(link = "log"),
    spatial = "off",
    control = sdmTMBcontrol(newton_loops = 1)
  )
  expect_true(all(!is.na(summary(m$sd_report)[, "Std. Error"])))
  # expect_length(residuals(m), nrow(d))

  # test non-spatial model
  set.seed(123)
  d <- pcod[pcod$density > 0, ]
  d$cluster <- sample(1:2, size = nrow(d), replace = TRUE, prob = c(0.9, 0.1))
  d$y <- rnorm(n = nrow(d), c(1, 4)[d$cluster], sd = 0.1)
  m <- sdmTMB(
    data = d, formula = y ~ 1,
    family = gamma_mix(),
    spatial = "off"
  )
  expect_equal(m$model$par[["logit_p_mix"]], stats::qlogis(0.1), tolerance = 0.1)
  expect_equal(exp(m$model$par[["log_ratio_mix"]]), 3.0, tolerance = 0.01)

  p <- predict(m, newdata = m$data)
  expect_equal(mean(p$est), log(mean(d$y)), tolerance = 0.001)

  # lognormal
  m <- sdmTMB(
    data = d, formula = y ~ 1,
    family = lognormal_mix(),
    spatial = "off"
  )
  expect_equal(m$model$par[["logit_p_mix"]], stats::qlogis(0.1), tolerance = 0.1)
  expect_equal(exp(m$model$par[["log_ratio_mix"]]), 3.0, tolerance = 0.01)

  # NB2
  set.seed(123)
  mix_ratio <- 10
  y_small <- rnbinom(5000, size = 2, mu = 4)
  y_large <- rnbinom(5000, size = 2, mu = 4 * mix_ratio)
  cluster <- sample(1:2, size = length(y_small), replace = TRUE, prob = c(0.9, 0.1))
  y <- ifelse(cluster == 1, y_small, y_large)
  d <- data.frame(y = y)
  m <- sdmTMB(
    data = d, formula = y ~ 1,
    family = nbinom2_mix(),
    spatial = "off"
  )
  expect_equal(m$model$par[["logit_p_mix"]], stats::qlogis(0.1), tolerance = 0.1)
  expect_equal(1 + exp(m$model$par[["log_ratio_mix"]]), mix_ratio, tolerance = 0.1)

  p <- predict(m, newdata = m$data)
  expect_equal(mean(p$est), log(mean(d$y)), tolerance = 0.001)
})

test_that("Test that residuals and prediction functions work with mixture models", {
  skip_on_ci()
  skip_on_cran()
  d <- pcod[pcod$year == 2017 & pcod$density > 0, ]
  m <- sdmTMB(
    data = d, formula = density ~ 1,
    family = lognormal_mix(link = "log"),
    spatial = "off"
  )
  expect_true(all(!is.na(summary(m$sd_report)[, "Std. Error"])))
  # expect_length(residuals(m), nrow(d))
  expect_length(predict(m)[["est"]], nrow(d))
  p <- predict(m, newdata = m$data)
  expect_equal(mean(p$est), log(mean(d$density)), tolerance = 0.01)
})

test_that("Test that delta Gamma mixture fits", {
  skip_on_ci()
  skip_on_cran()
  d <- pcod
  m <- sdmTMB(
    data = d, formula = density ~ 1,
    family = delta_gamma_mix(),
    spatial = "off"
  )
  p <- predict(m, newdata = m$data, type = "response")
  expect_equal(mean(p$est), mean(d$density), tolerance = 0.01)
  # expect_length(residuals(m), nrow(d))
  # set.seed(123)
  # d$test_gamma <- stats::rgamma(nrow(d), shape = 0.5, scale = 1 / 0.5)
  # m <- sdmTMB(data = d, formula = test_gamma ~ 1,
  #             mesh = spde, family = gamma_mix(link = "inverse"), spatiotemporal = "off",
  #             control = sdmTMBcontrol(newton_loops = 1))
})

test_that("Test that delta lognormal mixture fits", {
  skip_on_ci()
  skip_on_cran()

  set.seed(1)
  y1 <- stats::rlnorm(1000, 0, 0.5)
  y2 <- stats::rlnorm(1000, 1, 0.5)
  d <- data.frame(y = c(y1, y2))

  m <- sdmTMB(
    data = d,
    formula = y ~ 1,
    family = delta_lognormal_mix(),
    spatial = "off"
  )
  p <- predict(m, newdata = m$data, type = "response")
  expect_equal(mean(p$est), mean(d$y), tolerance = 0.01)
  # expect_length(residuals(m, model = 2), nrow(d))
})

test_that("Test that simulation functions work with mixture models", {
  skip_on_ci()
  skip_on_cran()

  set.seed(123)
  d <- pcod[pcod$density > 0, ]
  d$cluster <- sample(1:2, size = nrow(d), replace = T, prob = c(0.7, 0.3))
  phi <- 0.8
  d$y <- stats::rgamma(nrow(d), shape = phi, scale = c(2, 10)[d$cluster] / phi)

  m_gamma <- sdmTMB(
    data = d, formula = density ~ 1,
    family = gamma_mix(link = "log"),
    spatial = "off",
    control = sdmTMBcontrol(newton_loops = 1)
  )

  set.seed(123)
  d$cluster <- sample(1:2, size = nrow(d), replace = T, prob = c(0.7, 0.3))
  d$y <- stats::rlnorm(nrow(d), meanlog = log(c(2, 10))[d$cluster], sdlog = 0.5)

  m_logn <- sdmTMB(
    data = d, formula = y ~ 1,
    family = lognormal_mix(link = "log"),
    spatial = "off",
    control = sdmTMBcontrol(newton_loops = 1)
  )
  expect_length(simulate(m_gamma), nrow(d))
  expect_length(simulate(m_logn), nrow(d))
})
