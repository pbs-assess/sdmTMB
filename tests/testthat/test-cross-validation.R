context("Cross validation")

test_that("Basic cross validation works", {
  skip_on_ci()
  skip_on_cran()
  d <- subset(pcod, year >= 2011) # subset for example speed
  spde <- make_mesh(d, c("X", "Y"), cutoff = 20)

  set.seed(1)
  # library(future) # for parallel processing
  # plan(multisession) # for parallel processing
  x <- sdmTMB_cv(
    density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
    data = d, spde = spde,
    family = tweedie(link = "log"), time = "year", k_folds = 2
  )
  expect_equal(class(x$sum_loglik), "numeric")
  expect_true("data.frame" %in% class(x$data))
  expect_equal(class(x$models[[1]]), "sdmTMB")
})
