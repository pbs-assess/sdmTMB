test_that("Basic cross validation works", {
  skip_on_ci()
  skip_on_cran()
  skip_if_not_installed("INLA")
  d <- pcod
  spde <- make_mesh(d, c("X", "Y"), cutoff = 15)

  set.seed(2)
  # library(future) # for parallel processing
  # plan(multisession) # for parallel processing
  x <- sdmTMB_cv(
    density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
    data = d, mesh = spde,
    family = tweedie(link = "log"), time = "year", k_folds = 2
  )
  expect_equal(class(x$sum_loglik), "numeric")
  expect_equal(x$sum_loglik, sum(x$data$cv_loglik))
  expect_equal(x$sum_loglik, sum(x$fold_loglik))
  expect_true("data.frame" %in% class(x$data))
  expect_equal(class(x$models[[1]]), "sdmTMB")

  # Use fold_ids:
  x <- sdmTMB_cv(
    density ~ 1,
    data = d, mesh = spde, spatial = "off",
    family = tweedie(link = "log"),
    fold_ids = rep(seq(1, 2), nrow(d))[seq(1, nrow(d))])
  expect_equal(class(x$models[[1]]), "sdmTMB")

  # student-t: was broken at one time, must deal with `df`
  d <- subset(d, density > 0)
  d$log_density <- log(d$density)
  spde <- make_mesh(d, c("X", "Y"), cutoff = 15)
  x <- sdmTMB_cv(
    log_density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
    data = d, mesh = spde, spatial = "off",
    family = sdmTMB::student(df = 9), time = "year", k_folds = 2
  )
})
