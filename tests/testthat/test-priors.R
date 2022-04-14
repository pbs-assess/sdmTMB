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
  d <- pcod_2011
  pcod_spde <- pcod_mesh_2011

  # no priors
  m <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
    data = d, time = "year", mesh = pcod_spde, family = tweedie(link = "log"),
    spatiotemporal = "AR1",
    share_range = FALSE
  )

  # population-effects missing a prior; should error
  expect_error(
    {
      mp <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
        data = d, time = "year", mesh = pcod_spde, family = tweedie(link = "log"),
        share_range = FALSE, spatiotemporal = "AR1",
        priors = sdmTMBpriors(
          b = normal(c(0, 0, NA, NA, NA), c(2, 2, NA, NA, NA))
        )
      )
    },
    regexp = "prior"
  )

  # all the priors
  mp <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
    data = d, time = "year", mesh = pcod_spde, family = tweedie(link = "log"),
    share_range = FALSE, spatiotemporal = "AR1",
    priors = sdmTMBpriors(
      b = normal(c(0, 0, NA, NA, NA, NA), c(2, 2, NA, NA, NA, NA)),
      phi = halfnormal(0, 10),
      tweedie_p = normal(1.5, 2),
      ar1_rho = normal(0, 1),
      matern_s = pc_matern(range_gt = 5, sigma_lt = 1),
      matern_st = pc_matern(range_gt = 5, sigma_lt = 1)
    )
  )

  expect_lt(abs(mp$model$par[1]), abs(m$model$par[1]))
  expect_gt(mp$model$par[["ln_tau_O"]], m$model$par[["ln_tau_O"]])
  expect_gt(mp$model$par[["ln_tau_E"]], m$model$par[["ln_tau_E"]])
})


test_that("Additional priors work", {
  skip_on_ci()
  skip_on_cran()
  skip_if_not_installed("INLA")
  d <- pcod_2011
  pcod_spde <- pcod_mesh_2011

  # priors with one covariate/intercept only and the normal scale not 1:
  m_norm <- sdmTMB(density ~ 1,
    data = d, mesh = pcod_spde, family = tweedie(link = "log"),
    priors = sdmTMBpriors(b = normal(0, 10)),
    spatial = "off", spatiotemporal = "off"
  )

  # univariate normal priors
  m_norm <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
    data = d, mesh = pcod_spde, family = tweedie(link = "log"),
    priors = sdmTMBpriors(b = normal(rep(0, 6), rep(1, 6))),
    spatial = "off", spatiotemporal = "off"
  )
  expect_identical(class(m_norm), "sdmTMB")
  m_mvn <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
    data = d, mesh = pcod_spde, family = tweedie(link = "log"),
    priors = sdmTMBpriors(b = mvnormal(rep(0, 6), diag(1, 6))),
    spatial = "off", spatiotemporal = "off"
  )
  expect_identical(class(m_mvn), "sdmTMB")
  m_mvn_na <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
    data = d, mesh = pcod_spde, family = tweedie(link = "log"),
    priors = sdmTMBpriors(b = mvnormal(c(NA, 0, 0, 0, 0, 0), diag(1, 6))),
    spatial = "off", spatiotemporal = "off"
  )
  expect_identical(class(m_mvn_na), "sdmTMB")

  m_norm_na <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
    data = d, mesh = pcod_spde, family = tweedie(link = "log"),
    priors = sdmTMBpriors(b = normal(c(NA, rep(0, 5)), c(NA, rep(1, 5)))),
    spatial = "off", spatiotemporal = "off"
  )
  expect_identical(class(m_norm_na), "sdmTMB")

  expect_equal(tidy(m_norm), tidy(m_mvn), tolerance = 0.0001)
  expect_equal(tidy(m_norm_na), tidy(m_mvn_na), tolerance = 0.0001)
  expect_true(
    abs(tidy(m_mvn)$estimate[tidy(m_mvn)$term == "depth_scaled"]) <
      abs(tidy(m_mvn_na)$estimate[tidy(m_mvn_na)$term == "depth_scaled"])
  )
})
