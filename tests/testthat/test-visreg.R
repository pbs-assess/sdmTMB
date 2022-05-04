test_that("visreg works", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  skip_if_not_installed("visreg")

  pcod_2011$fyear <- as.factor(pcod_2011$year)
  fit <- sdmTMB(
    density ~ s(depth_scaled) + fyear,
    data = pcod_2011, mesh = pcod_mesh_2011,
    spatial = "off",
    family = tweedie()
  )
  visreg::visreg(fit, xvar = "depth_scaled")
  visreg::visreg(fit, xvar = "fyear")
  visreg::visreg(fit, xvar = "depth_scaled", scale = "response")
  visreg::visreg2d(fit, xvar = "fyear", yvar = "depth_scaled")
  v <- visreg::visreg(fit, xvar = "depth_scaled")
  expect_identical(class(v), "visreg")

  # Delta model example:
  fit_dg <- sdmTMB(
    density ~ s(depth_scaled, year, k = 8),
    data = pcod_2011, mesh = pcod_mesh_2011,
    spatial = "off",
    family = delta_gamma()
  )
  visreg_delta(fit_dg, xvar = "depth_scaled", model = 1)
  visreg_delta(fit_dg, xvar = "depth_scaled", model = 2)
  visreg_delta(fit_dg,
    xvar = "depth_scaled", model = 1,
    scale = "response"
  )
  visreg_delta(fit_dg,
    xvar = "depth_scaled", model = 2,
    scale = "response"
  )
  visreg2d_delta(fit_dg,
    xvar = "depth_scaled", yvar = "year",
    model = 2, scale = "response"
  )
  visreg2d_delta(fit_dg,
    xvar = "depth_scaled", yvar = "year",
    model = 1, scale = "response"
  )
  visreg2d_delta(fit_dg, xvar = "depth_scaled", yvar = "year", model = 1)
  visreg2d_delta(fit_dg, xvar = "depth_scaled", yvar = "year", model = 2)

  v <- visreg_delta(fit_dg, xvar = "depth_scaled", model = 1)
  expect_identical(class(v), "visreg")

  v <- visreg2d_delta(fit_dg, xvar = "depth_scaled", yvar = "year", model = 2)
  expect_identical(class(v), "visreg2d")
})
