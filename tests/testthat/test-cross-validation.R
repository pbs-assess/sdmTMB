test_that("Basic cross validation works", {
  skip_on_ci()
  skip_on_cran()
  skip_if_not_installed("INLA")
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
  expect_equal(x$sum_loglik, sum(x$data$cv_loglik))
  expect_equal(x$sum_loglik, sum(x$fold_loglik))
  expect_true("data.frame" %in% class(x$data))
  expect_equal(class(x$models[[1]]), "sdmTMB")

  # Use fold_ids:
  x <- sdmTMB_cv(
    density ~ 0 + depth_scaled + depth_scaled2,
    data = pcod, spde = spde,
    family = tweedie(link = "log"),
    fold_ids = rep(seq(1, 2), nrow(pcod))[seq(1, nrow(pcod))])
  expect_equal(class(x$models[[1]]), "sdmTMB")
})
