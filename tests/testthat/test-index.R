test_that("coef works", {
  skip_on_ci()
  skip_on_cran()

  pcod_spde <- make_mesh(pcod, c("X", "Y"), n_knots = 50, type = "kmeans")
  m <- sdmTMB(
    data = pcod,
    formula = density ~ 0 + as.factor(year),
    spatiotemporal = "off", # speed
    time = "year", mesh = pcod_spde,
    family = tweedie(link = "log")
  )
  nd <- replicate_df(qcs_grid, "year", unique(pcod$year))
  predictions <- predict(m, newdata = nd, return_tmb_object = TRUE)
  ind <- get_index(predictions)
  ind
  expect_s3_class(ind, "data.frame")

  ind <- get_index(predictions, bias_correct = TRUE)
  expect_s3_class(ind, "data.frame")

  cog <- get_cog(predictions)
  cog
  expect_s3_class(cog, "data.frame")

  cog <- get_cog(predictions, format = "wide")
  cog
  expect_s3_class(cog, "data.frame")

  expect_error(get_index(predictions, area = c(1, 2, 3)), regexp = "area")

  nd <- replicate_df(qcs_grid, "year", c(2003, 2004))

  expect_warning(
    predictions <- predict(m, newdata = nd, return_tmb_object = TRUE)
  )
  expect_error(ind <- get_index(predictions), regexp = "time")
})
