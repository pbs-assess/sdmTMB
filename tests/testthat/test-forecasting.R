test_that("Forecasting works with a time-varying parameter", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  spde <- make_mesh(pcod, c("X", "Y"), n_knots = 50, type = "kmeans")
  grid2019 <- qcs_grid[qcs_grid$year == max(qcs_grid$year), ]
  grid2019$year <- 2019L
  qcsgrid_forecast <- rbind(qcs_grid, grid2019)
  m <- sdmTMB(
    data = pcod, formula = density ~ 0,
    time_varying = ~ 1,
    fields = "AR1",
    extra_time = c(2016L, 2018L, 2019L),
    include_spatial = FALSE,
    time = "year",
    spde = spde,
    control = sdmTMBcontrol(mgcv = FALSE),
    family = tweedie(link = "log")
  )
  print(m)
  predictions <- predict(m, newdata = qcsgrid_forecast)
  expect_true("data.frame" %in% class(predictions))
})
