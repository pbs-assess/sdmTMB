test_that("Gamma and lognormal mixtures fit and recover mixing proportion", {
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
  expect_length(residuals(m), nrow(d))

  # test non-spatial model
  set.seed(123)
  d <- pcod[pcod$density > 0, ]
  d$cluster <- sample(1:2, size = nrow(d), replace = T, prob = c(0.9, 0.1))
  d$y <- rnorm(n = nrow(d), c(1, 4)[d$cluster], sd = 0.1)
  m <- sdmTMB(
    data = d, formula = y ~ 1,
    family = gamma_mix(),
    spatial = "off",
    control = sdmTMBcontrol(newton_loops = 1)
  )
  expect_equal(m$model$par[["logit_p_mix"]], -2.2, tolerance = 0.1)
  expect_equal(exp(m$model$par[["log_ratio_mix"]]), 3.0, tolerance = 0.01)

  m <- sdmTMB(
    data = d, formula = y ~ 1,
    family = lognormal_mix(),
    spatial = "off",
    control = sdmTMBcontrol(newton_loops = 1)
  )
  expect_equal(m$model$par[["logit_p_mix"]], -2.2, tolerance = 0.1)
  # set.seed(123)
  # d$test_gamma <- stats::rgamma(nrow(d), shape = 0.5, scale = 1 / 0.5)
  # m <- sdmTMB(data = d, formula = test_gamma ~ 1,
  #             mesh = spde, family = gamma_mix(link = "inverse"), spatiotemporal = "off",
  #             control = sdmTMBcontrol(newton_loops = 1))
})

test_that("Test that residuals and prediction functions work with mixture models", {
  skip_on_ci()
  skip_on_cran()
  d <- pcod[pcod$year == 2017 & pcod$density > 0, ]
  spde <- make_mesh(d, c("X", "Y"), cutoff = 10)
  m <- sdmTMB(
    data = d, formula = density ~ 1,
    mesh = spde, family = lognormal_mix(link = "log"),
    spatial = "off",
    control = sdmTMBcontrol(newton_loops = 1)
  )
  expect_true(all(!is.na(summary(m$sd_report)[, "Std. Error"])))
  expect_length(residuals(m), nrow(d))
  expect_length(predict(m)[["est"]], nrow(d))
})

test_that("Test that delta Gamma mixture fits", {
  skip_on_ci()
  skip_on_cran()
  d <- pcod
  m <- sdmTMB(
    data = d, formula = density ~ 1,
    family = delta_gamma_mix(),
    spatial = "off",
    control = sdmTMBcontrol(newton_loops = 1)
  )
  expect_length(residuals(m), nrow(d))
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
    spatial = "off",
    control = sdmTMBcontrol(newton_loops = 1)
  )
  expect_length(residuals(m), nrow(d))
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
