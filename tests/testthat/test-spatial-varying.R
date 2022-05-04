test_that("Spatially-varying coefficients are estimated correctly for binomial and delta models", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  local_edition(2)
  d <- pcod
  d$year_scaled <- as.numeric(scale(d$year))
  mesh10 <- make_mesh(d, c("X", "Y"), cutoff = 10)
  m1 <- sdmTMB(
    data = d,
    formula = present ~ 1 + year_scaled,
    spatial_varying = ~ 0 + year_scaled,
    mesh = mesh10,
    family = binomial()
  )
  b1 <- tidy(m1, effects = "ran_pars", conf.int = T)
  expect_equal(b1$estimate[3], 0.312, tolerance = 0.1)

  # also check that binomial portion of delta model matches the above
  m2 <- sdmTMB(
    data = d,
    formula = density ~ 1 + year_scaled,
    spatial_varying = ~ 0 + year_scaled,
    mesh = mesh10,
    family = delta_gamma()
  )
  b2 <- tidy(m2, effects = "ran_pars", conf.int = T)
  expect_equal(b2$estimate[3], b1$estimate[3], tolerance = 1e-3)
})
