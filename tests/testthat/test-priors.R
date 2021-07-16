test_that("Basic prior parsing works", {
  expect_equal(normal(0, 1), matrix(c(0, 1), ncol = 2L), ignore_attr = TRUE)
  expect_equal(halfnormal(0, 1), matrix(c(0, 1), ncol = 2L), ignore_attr = TRUE)
  expect_equal(pc_matern(5, 5, 0.05, 0.05), c(5, 5, 0.05, 0.05), ignore_attr = TRUE)

  expect_error(normal(NA, 1))
  expect_error(normal(0, -1))
  expect_error(normal(1, NA))
  expect_error(normal(c(1, 1), 1))

  expect_error(pc_matern(NA, 1))
  expect_error(pc_matern(1, NA))
  expect_error(pc_matern(-1, 1))
  expect_error(pc_matern(1, -1))
  expect_error(pc_matern(1, 1, 1.5, 0.05))
  expect_error(pc_matern(1, 1, -1.5, 0.05))
  expect_error(pc_matern(1, 1, 0.05, -0.05))
  expect_error(pc_matern(1, 1, 0.05, NA))
})

test_that("Prior fitting works", {
  skip_on_ci()
  skip_on_cran()
  skip_if_not_installed("INLA")
  d <- subset(pcod, year >= 2013)
  pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 30)

  # no priors
  m <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
    data = d, time = "year", spde = pcod_spde, family = tweedie(link = "log"),
    fields = "AR1",
    share_range = FALSE)

  # all the priors
  mp <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
    data = d, time = "year", spde = pcod_spde, family = tweedie(link = "log"),
    share_range = FALSE, fields = "AR1",
    priors = sdmTMBpriors(
      b = normal(c(0, 0, NA, NA, NA), c(2, 2, NA, NA, NA)),
      phi = halfnormal(0, 10),
      tweedie_p = normal(1.5, 2),
      ar1_rho = normal(0, 1),
      matern_s = pc_matern(range_gt = 5, sigma_lt = 1),
      matern_st = pc_matern(range_gt = 5, sigma_lt = 1))
  )

  expect_lt(abs(mp$model$par[1]), abs(m$model$par[1]))
  expect_gt(mp$model$par[["ln_tau_O"]], m$model$par[["ln_tau_O"]])
  expect_gt(mp$model$par[["ln_tau_E"]], m$model$par[["ln_tau_E"]])
})
