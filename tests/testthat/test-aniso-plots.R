test_that("aniso plots work", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggforce")

  mesh <- make_mesh(pcod_2011, c("X", "Y"), n_knots = 70, type = "kmeans")
  fit <- sdmTMB(
    data = pcod_2011,
    formula = density ~ 1,
    mesh = mesh,
    family = tweedie(),
    share_range = FALSE,
    time = "year",
    anisotropy = TRUE #<
  )
  g <- plot_anisotropy(fit)
  expect_s3_class(g, "ggplot")

  d <- plot_anisotropy(fit, return_data = TRUE)
  expect_s3_class(d, "data.frame")

  plot_anisotropy2(fit)

  mesh <- make_mesh(pcod_2011, c("X", "Y"), n_knots = 60, type = "kmeans")
  fit <- sdmTMB(
    data = pcod_2011,
    formula = density ~ 1,
    mesh = mesh,
    family = delta_gamma(),
    anisotropy = TRUE #<
  )
  g <- plot_anisotropy(fit)
  expect_s3_class(g, "ggplot")

  plot_anisotropy2(fit)
  plot_anisotropy2(fit, model = 2)
})
