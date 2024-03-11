test_that("get_index works", {
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
})


test_that("get_index works with subsets of years", {
  skip_on_ci()
  skip_on_cran()

  m <- sdmTMB(
    density ~ 0 + as.factor(year),
    data = pcod_2011,
    time = "year",
    spatiotemporal = "off",
    spatial = 'off',
    family = delta_gamma()
  )
  nd <- replicate_df(qcs_grid, "year", unique(pcod_2011$year))
  p_full <- predict(m, newdata = nd, return_tmb_object = TRUE)

  nd_2011 <- replicate_df(qcs_grid, "year", 2011)
  p_2011 <- predict(m, newdata = nd_2011, return_tmb_object = TRUE)

  nd_2 <- replicate_df(qcs_grid, "year", c(2011, 2013))
  p_2 <- predict(m, newdata = nd_2, return_tmb_object = TRUE)

  nd_3 <- replicate_df(qcs_grid, "year", c(2015, 2011))
  p_3 <- predict(m, newdata = nd_3, return_tmb_object = TRUE)

  index_full <- get_index(p_full, bias_correct = TRUE)
  expect_equal(index_full$est, c(322529.7268, 293854.3696, 390942.2649, 184368.2756), tolerance = 0.01)
  index_2011 <- get_index(p_2011, bias_correct = TRUE)
  index_2 <- get_index(p_2, bias_correct = TRUE)
  index_3 <- get_index(p_3, bias_correct = TRUE)

  expect_equal(index_2011$est, subset(index_full, year == 2011)$est)
  expect_equal(index_2$est, subset(index_full, year %in% c(2011, 2013))$est)
  expect_equal(index_3$est, subset(index_full, year %in% c(2015, 2011))$est)

  index_apply <- lapply(unique(pcod_2011$year), \(y) {
    nd <- replicate_df(qcs_grid, "year", y)
    p <- predict(m, newdata = nd, return_tmb_object = TRUE)
    get_index(p, bias_correct = TRUE)
  })
  index_apply <- do.call(rbind, index_apply)
  expect_equal(index_apply, index_full)
})

# test_that("get_index faster epsilon bias correction", {
#   skip_on_ci()
#   skip_on_cran()
#
#   library(sdmTMB)
#   mesh <- make_mesh(pcod, c('X', 'Y'), cutoff = 5)
#
#   m <- sdmTMB(
#     density ~ factor(year),
#     data = pcod,
#     mesh = mesh,
#     time = "year",
#     family = delta_gamma()
#   )
#   nd <- replicate_df(qcs_grid, "year", unique(pcod$year))
#   p <- predict(m, newdata = nd, return_tmb_object = TRUE)
#
#   TRACE=TRUE
#
#   INTERN=TRUE
#   LOWRANK=TRUE
#   index <- get_index(p, bias_correct = TRUE)
#
#   INTERN=FALSE
#   LOWRANK=FALSE
#   index <- get_index(p, bias_correct = TRUE)
#
#   INTERN=FALSE
#   LOWRANK=TRUE
#   index <- get_index(p, bias_correct = TRUE)
#
#   INTERN=TRUE
#   LOWRANK=FALSE
#   index <- get_index(p, bias_correct = TRUE)
#
# })

