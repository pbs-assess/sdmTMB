context("Index and COG calculations work")

test_that("make_spde works but throws a warning", {
  expect_warning({pcod_spde <- make_spde(pcod$X, pcod$Y, n_knots = 30)})
})

pcod_spde <- make_mesh(pcod, c("X", "Y"), cutoff = 20)
m <- sdmTMB(
  data = pcod,
  formula = density ~ 0 + as.factor(year),
  time = "year", spde = pcod_spde, family = tweedie(link = "log")
)
predictions <- predict(m, newdata = qcs_grid, return_tmb_object = TRUE)

test_that("get_index() works", {
  ind <- get_index(predictions, bias_correct = FALSE)
  expect_gte(nrow(ind), 5L)
  expect_equal(class(ind), "data.frame")
})

test_that("get_cog() works", {
  cog <- get_cog(predictions, bias_correct = FALSE)
  expect_gte(nrow(cog), 5L)
  expect_equal(class(cog), "data.frame")
  p <- predict(m, newdata = qcs_grid, return_tmb_object = FALSE)
  expect_error(get_cog(p), regexp = "return_tmb_object")
})
