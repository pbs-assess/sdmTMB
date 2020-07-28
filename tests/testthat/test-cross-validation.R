context("Cross validation")

test_that("Basic cross validation works", {
  d <- subset(pcod, year >= 2015) # subset for example speed
  x <- sdmTMB_cv(
    formula = density ~ 1,
    d, family = tweedie(link = "log"), time = "year", x = "X", y = "Y",
    n_knots = 20, k_folds = 2, mgcv = FALSE
  )
  expect_equal(class(x$sum_loglik), "numeric")
  expect_true("data.frame" %in% class(x$data))
  expect_equal(class(x$models[[1]]), "sdmTMB")
})

test_that("Cross validation works with Gaussian and no time element", {
  set.seed(123)
  s <- sim(time_steps = 1, sigma_O = 0.2)
  x <- sdmTMB_cv(formula = observed ~ 0, data = s, x = "x", y = "y",
    family = gaussian(link = "identity"), n_knots = 30, k_folds = 2,
    mgcv = FALSE)
  expect_equal(class(x$sum_loglik), "numeric")
})
