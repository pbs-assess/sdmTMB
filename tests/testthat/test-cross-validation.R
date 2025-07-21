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
  print(x)
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
  skip_on_ci()
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
  skip_on_ci()
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

test_that("LFO fold assignments follow documented structure", {
  skip_on_cran()
  skip_on_ci()
  
  # Test the documented example from help: 9 time steps, lfo_forecast = 2, lfo_validations = 3
  # Expected:
  # - Fit data to time steps 1 to 5, predict and validate step 7.
  # - Fit data to time steps 1 to 6, predict and validate step 8. 
  # - Fit data to time steps 1 to 7, predict and validate step 9.
  
  # Use pcod data which has 9 time steps: 2003, 2004, 2005, 2007, 2009, 2011, 2013, 2015, 2017
  time_steps <- sort(unique(pcod$year))
  expect_equal(length(time_steps), 9)
  
  # This should work if we have enough validation periods
  # With lfo_forecast = 2 and lfo_validations = 2 instead of 3
  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 25)
  x <- sdmTMB_cv(
    present ~ 1,
    data = pcod,
    mesh = mesh,
    lfo = TRUE,
    lfo_forecast = 2,
    lfo_validations = 2,  # Reduced to 2 to fit in 9 time steps
    family = binomial(),
    spatiotemporal = "off",
    time = "year"
  )
  
  # Check fold structure
  fold_table <- table(x$data$cv_fold, x$data$year)
  
  # Validation 1: Should train on years 2003-2011 (time steps 1-5), validate on 2015 (step 7)
  # Validation 2: Should train on years 2003-2013 (time steps 1-6), validate on 2017 (step 8)
  
  # Check that validation years have the right fold assignments
  # For lfo_forecast = 2, validation data should be in folds k + lfo_forecast
  expect_true(all(x$data$cv_fold[x$data$year == 2015] == 3))  # fold 1 + 2
  expect_true(all(x$data$cv_fold[x$data$year == 2017] == 4))  # fold 2 + 2
  
  # Check that training years are properly assigned
  # For validation 1: years 2003-2011 should have fold <= 1
  train_years_1 <- c(2003, 2004, 2005, 2007, 2009, 2011)
  expect_true(all(x$data$cv_fold[x$data$year %in% train_years_1] <= 1))
  
  # For validation 2: years 2003-2013 should have fold <= 2  
  train_years_2 <- c(2003, 2004, 2005, 2007, 2009, 2011, 2013)
  expect_true(all(x$data$cv_fold[x$data$year %in% train_years_2] <= 2))
})

test_that("LFO parameter validation works", {
  skip_on_cran()
  skip_on_ci()
  
  mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 25)
  
  # Should error when lfo_validations + lfo_forecast > number of time steps
  expect_error(
    sdmTMB_cv(
      present ~ 1,
      data = pcod_2011,  # Only 2 time steps  
      mesh = mesh,
      lfo = TRUE,
      lfo_forecast = 2,
      lfo_validations = 3,  # 3 + 2 = 5 > 2 time steps
      family = binomial(),
      time = "year"
    ),
    "Not enough time steps"
  )
})

test_that("Cross validation with offsets works", {
  skip_on_cran()
  skip_on_ci()
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
  skip_on_ci()
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
