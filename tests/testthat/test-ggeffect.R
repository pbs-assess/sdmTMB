test_that("ggeffects + sdmTMB", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  skip_if_not_installed("ggeffects")
  skip_if_not_installed("ggplot2")

  d <- pcod_2011
  d$fyear <- as.factor(d$year)
  fit <- sdmTMB(present ~ depth_scaled + I(depth_scaled^2) + fyear,
    data = d,
    mesh = pcod_mesh_2011,
    family = binomial()
  )
  g <- ggeffects::ggeffect(fit, "depth_scaled [-2.5:2.5, by=.1]")
  expect_true(inherits(g, "ggeffects"))
  plot(g)

  e <- effects::effect("depth_scaled", fit)
  e
  expect_true(inherits(e, "eff"))

  fit2 <- sdmTMB(present ~ depth_scaled + I(depth_scaled^2) + (1 | fyear),
    data = d,
    mesh = pcod_mesh_2011,
    family = binomial()
  )
  effects::effect("depth_scaled", fit2)
  g <- ggeffects::ggeffect(fit2, "depth_scaled [-2.5:2.5, by=.1]")
  plot(g)
  expect_true(inherits(g, "ggeffects"))

  fit3 <- sdmTMB(present ~ s(year, k = 3) + depth_scaled + I(depth_scaled^2),
    data = d,
    mesh = pcod_mesh_2011,
    family = binomial()
  )
  effects::effect("depth_scaled", fit3)
  g <- ggeffects::ggeffect(fit3, "depth_scaled [-2.5:2.5, by=.1]")
  plot(g)
  expect_true(inherits(g, "ggeffects"))

  fit4 <- sdmTMB(present ~ s(depth_scaled, k = 4),
    data = d,
    mesh = pcod_mesh_2011,
    family = binomial()
  )
  expect_error({
    effects::effect("depth_scaled", fit4)
  }, regexp = "missing")
})
