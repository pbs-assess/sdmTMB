test_that("Forecasting works with a time-varying parameter", {
  skip_on_cran()
  spde <- make_mesh(pcod, c("X", "Y"), n_knots = 50, type = "kmeans")
  nd <- replicate_df(qcs_grid, "year", unique(pcod$year))
  grid2019 <- nd[nd$year == max(nd$year), ]
  grid2019$year <- 2019L
  qcsgrid_forecast <- rbind(nd, grid2019)
  grid2019$year <- 2018L
  qcsgrid_forecast <- rbind(qcsgrid_forecast, grid2019)
  grid2019$year <- 2016L
  qcsgrid_forecast <- rbind(qcsgrid_forecast, grid2019)

  m <- sdmTMB(
    data = pcod, formula = density ~ 0,
    time_varying = ~ 1,
    spatiotemporal = "AR1",
    extra_time = c(2016L, 2018L, 2019L),
    spatial = FALSE,
    time = "year",
    mesh = spde,
    family = tweedie(link = "log")
  )
  print(m)
  predictions <- predict(m, newdata = qcsgrid_forecast)
  expect_true("data.frame" %in% class(predictions))
})
