test_that("Basic cross validation works", {
  skip_on_ci()
  skip_on_cran()
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
  expect_equal(class(x$models[[1]]), "sdmTMB")

  # Try passing family as a variable -- this is per Issue #127
  fam <- gaussian(link = "identity")
  x <- sdmTMB_cv(
    log_density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
    data = d, mesh = spde, spatial = "off",
    family = fam, time = "year", k_folds = 2
  )
  expect_equal(class(x$models[[1]]), "sdmTMB")
})

test_that("Cross validation in parallel with globals", {
  skip_on_cran()
  # https://github.com/pbs-assess/sdmTMB/issues/127
  d <- pcod
  spde <- make_mesh(d, c("X", "Y"), cutoff = 15)
  set.seed(2)
  future::plan(future::multisession, workers = 2L)
  fam <- tweedie(link = "log")
  x <- sdmTMB_cv(
    density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
    data = d, mesh = spde,
    family = fam, time = "year", k_folds = 2L, future_globals = 'fam'
  )
  expect_s3_class(x$models[[1]], "sdmTMB")
  future::plan(future::sequential)
})

test_that("Leave future out cross validation works", {
  skip_on_cran()
  x <- sdmTMB_cv(
    present ~ 1,
    data = pcod_2011,
    mesh = pcod_mesh_2011,
    lfo = TRUE,
    lfo_forecast = 1,
    lfo_validations = 2,
    family = binomial(),
    time = "year"
  )
  expect_equal(class(x$sum_loglik), "numeric")
  expect_equal(x$sum_loglik, sum(x$data$cv_loglik))
  expect_equal(x$sum_loglik, sum(x$fold_loglik))
  expect_true("data.frame" %in% class(x$data))
  expect_true(inherits(x$models[[1]], "sdmTMB"))

  # Can see how the folds are assigned with
  table(x$models[[1]]$data$cv_fold, x$models[[1]]$data$year)

  expect_equal(length(x$models), 2)
  expect_equal(length(x$fold_loglik), 2)
  expect_equal(length(x$max_gradients), 2)
  expect_equal(cor(x$data$cv_fold[x$data$cv_fold!=1], x$data$year[x$data$cv_fold!=1]), 1.0)
})

test_that("Cross validation with offsets works", {
  skip_on_cran()
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")

  set.seed(1)
  d <- pcod_2011
  d$log_effort <- rnorm(nrow(d))

  # library(future)
  # future::plan(future::multisession)

  expect_error(
    fit_cv_off1 <- sdmTMB_cv(
      density ~ 1,
      data = d,
      mesh = pcod_mesh_2011,
      family = tweedie(),
      spatial = "off",
      offset = d$log_effort, #<
      k_folds = 2,
      parallel = TRUE #<
    ))

  set.seed(1)
  fit_cv_off1 <- sdmTMB_cv(
    density ~ 1,
    data = d,
    mesh = pcod_mesh_2011,
    family = tweedie(),
    spatial = "off",
    offset = "log_effort", #<
    k_folds = 2,
    parallel = TRUE #<
  )

  set.seed(1)
  fit_cv_off2 <- sdmTMB_cv(
    density ~ 1,
    data = d,
    mesh = pcod_mesh_2011,
    family = tweedie(),
    spatial = "off",
    offset = "log_effort", #<
    k_folds = 2,
    parallel = FALSE #<
  )

  set.seed(1)
  fit_cv_off3 <- sdmTMB_cv(
    density ~ 1,
    data = d,
    mesh = pcod_mesh_2011,
    family = tweedie(),
    spatial = "off",
    offset = "log_effort", #<
    k_folds = 2,
    use_initial_fit = TRUE, #<
    parallel = TRUE #<
  )

  set.seed(1)
  fit_cv_off4 <- sdmTMB_cv(
    density ~ 1,
    data = d,
    mesh = pcod_mesh_2011,
    family = tweedie(),
    spatial = "off",
    offset = "log_effort", #<
    k_folds = 2,
    use_initial_fit = TRUE, #<
    parallel = FALSE #<
  )

  # now without offset
  set.seed(1)
  fit_cv_off5 <- sdmTMB_cv(
    density ~ 1,
    data = d,
    mesh = pcod_mesh_2011,
    family = tweedie(),
    spatial = "off",
    k_folds = 2,
    use_initial_fit = TRUE, #<
    parallel = TRUE #<
  )

  set.seed(1)
  fit_cv_off6 <- sdmTMB_cv(
    density ~ 1,
    data = d,
    mesh = pcod_mesh_2011,
    family = tweedie(),
    spatial = "off",
    k_folds = 2,
    use_initial_fit = TRUE, #<
    parallel = FALSE #<
  )

  set.seed(1)
  fit_cv_off7 <- sdmTMB_cv(
    density ~ 1,
    data = d,
    mesh = pcod_mesh_2011,
    family = tweedie(),
    spatial = "off",
    k_folds = 2,
    use_initial_fit = FALSE, #<
    parallel = FALSE #<
  )

  set.seed(1)
  fit_cv_off8 <- sdmTMB_cv(
    density ~ 1,
    data = d,
    mesh = pcod_mesh_2011,
    family = tweedie(),
    spatial = "off",
    k_folds = 2,
    use_initial_fit = TRUE, #<
    parallel = FALSE #<
  )

  expect_equal(round(fit_cv_off1$models[[1]]$model$par, 4), round(fit_cv_off2$models[[1]]$model$par, 4))
  expect_equal(round(fit_cv_off1$models[[1]]$model$par, 4), round(fit_cv_off3$models[[1]]$model$par, 4))
  expect_equal(round(fit_cv_off1$models[[1]]$model$par, 4), round(fit_cv_off4$models[[1]]$model$par, 4))

  # with/without offset:
  expect_false(identical(fit_cv_off1$models[[1]]$model$par, fit_cv_off5$models[[1]]$model$par))

  expect_equal(round(fit_cv_off5$models[[1]]$model$par, 4), round(fit_cv_off6$models[[1]]$model$par, 4))
  expect_equal(round(fit_cv_off5$models[[1]]$model$par, 4), round(fit_cv_off7$models[[1]]$model$par, 4))
  expect_equal(round(fit_cv_off5$models[[1]]$model$par, 4), round(fit_cv_off8$models[[1]]$model$par, 4))

  # future::plan(future::sequential)
})

test_that("Delta model cross validation works", {
  skip_on_cran()
  set.seed(1)
  out_tw <- sdmTMB_cv(
    density ~ depth_scaled,
    data = pcod_2011, mesh = pcod_mesh_2011, spatial = "off",
    family = tweedie(), k_folds = 2
  )
  set.seed(1)
  out_dg <- sdmTMB_cv(
    density ~ depth_scaled,
    data = pcod_2011, mesh = pcod_mesh_2011, spatial = "off",
    family = delta_gamma(), k_folds = 2
  )
  diff_ll <- out_tw$sum_loglik - out_dg$sum_loglik
  expect_equal(round(diff_ll, 4), round(-22.80799, 4))

  set.seed(1)
  out_dpg <- sdmTMB_cv(
    density ~ depth_scaled,
    data = pcod_2011, mesh = pcod_mesh_2011, spatial = "off",
    family = delta_gamma(type = "poisson-link"), k_folds = 2
  )
  diff_ll <- out_dpg$sum_loglik - out_dg$sum_loglik
  expect_equal(round(diff_ll, 4), round(-4.629497, 4))
})
