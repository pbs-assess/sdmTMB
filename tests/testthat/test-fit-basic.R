context("basic model fitting and prediction tests")

SEED <- 1
set.seed(SEED)
x <- stats::runif(100, -1, 1)
y <- stats::runif(100, -1, 1)
loc <- data.frame(x = x, y = y)
spde <- make_mesh(loc, c("x", "y"), cutoff = 0.02)

test_that("sdmTMB model fit with a covariate beta", {
  initial_betas <- 0.5
  range <- 0.1
  sigma_O <- 0.3 # SD of spatial process
  sigma_E <- 0.3 # SD of spatial process
  phi <- 0.1 # observation error
  s <- sdmTMB_sim(
    x = x, y = y, mesh = spde, time_steps = 6L, betas = initial_betas,
    phi = phi, range = range, sigma_O = sigma_O, sigma_E = sigma_E,
    seed = SEED
  )
  spde <- make_mesh(s, c("x", "y"), cutoff = 0.02)
  plot(spde)
  m <- sdmTMB(data = s, formula = observed ~ 0 + cov1, time = "time", spde = spde)
  expect_output(print(m), "fit by")
  expect_output(summary(m), "fit by")
  p <- as.list(m$model$par)
  r <- m$tmb_obj$report()
  est <- tidy(m, "ran_pars")
  expect_identical(round(est[,"estimate", drop = TRUE], 9),
    c(0.303614773, 0.308394419, 0.106559617, 0.064617947))
  expect_equal(m$model$convergence, 0L)
  expect_equal((p$b_j - initial_betas)^2, 0, tol = 0.001)
  expect_equal((exp(p$ln_phi) - phi)^2, 0, tol = 0.002)
  expect_equal((r$sigma_O - sigma_O)^2, 0, tol = 0.002)
  expect_equal((r$sigma_E[1] - sigma_E)^2, 0, tol = 0.001)
  expect_equal(est$estimate[est$term == "range"], range, tol = 0.01)
  p <- predict(m)
  r <- residuals(m)
  expect_equal(mean((p$est - s$observed)^2), 0, tol = 0.002)
})

test_that("Anisotropy fits and plots", {
  skip_on_ci()
  skip_on_cran()
  m <- sdmTMB(data = pcod,
    formula = density ~ 0 + as.factor(year),
    spde = make_mesh(pcod, c("X", "Y"), n_knots = 50, type = "kmeans"),
    family = tweedie(link = "log"), anisotropy = TRUE)
  expect_identical(class(m), "sdmTMB")
  plot_anisotropy(m)
  expect_identical(round(m$tmb_obj$report()$H, 8),
    structure(c(0.66544071, 0.07912087, 0.07912087, 1.51217096), .Dim = c(2L, 2L)))
})

test_that("Regularization works", {
  skip_on_cran()
  skip_on_ci()
  d <- subset(pcod, year >= 2015)
  d$depth_scaled <- as.numeric(scale(d$depth_scaled))
  m1 <- sdmTMB(data = d,
    formula = density ~ 0 + depth_scaled + as.factor(year),
    spde = make_mesh(d, c("X", "Y"), n_knots = 50, type = "kmeans"),
    family = tweedie(link = "log"))

  # Bypassing via NAs:
  m2 <- sdmTMB(data = d,
    formula = density ~ 0 + depth_scaled + as.factor(year),
    spde = make_mesh(d, c("X", "Y"), n_knots = 50, type = "kmeans"),
    family = tweedie(link = "log"),
    penalties = c(NA, NA, NA))
  expect_equal(m1$sd_report, m2$sd_report)

  # Ridge regression on depth term:
  m2 <- sdmTMB(data = d,
    formula = density ~ 0 + depth_scaled + as.factor(year),
    spde = make_mesh(d, c("X", "Y"), n_knots = 50, type = "kmeans"),
    family = tweedie(link = "log"),
    penalties = c(1, NA, NA))
  b1 <- tidy(m1)
  b2 <- tidy(m2)
  expect_lt(b2$estimate[1], b1$estimate[2])
  expect_lt(b2$std.error[1], b1$std.error[2])
})

test_that("A model with splines works", {
  skip_on_cran()
  skip_on_ci()
  d <- subset(pcod, year >= 2015)
  m <- sdmTMB(data = d,
    formula = density ~ 1 + s(depth_scaled),
    spde = make_mesh(d, c("X", "Y"), n_knots = 50, type = "kmeans"),
    family = tweedie(link = "log"))
  expect_identical(class(m), "sdmTMB")
})

test_that("A spatiotemporal version works with predictions on new data points", {
  d <- subset(pcod, year >= 2015)
  pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 30)
  m <- sdmTMB(
    data = d,
    formula = density ~ 0 + as.factor(year),
    time = "year", spde = pcod_spde, family = tweedie(link = "log"),
    include_spatial = FALSE
  )
  # Predictions at original data locations:
  predictions <- predict(m)
  predictions$resids <- residuals(m) # randomized quantile residuals
  # Predictions onto new data:
  predictions <- predict(m, newdata = subset(qcs_grid, year >= 2015))
  expect_identical(class(predictions), "data.frame")
})

test_that("Predictions on the original data set as `newdata`` return the same predictions", {
  skip_on_cran()
  skip_on_ci()
  set.seed(1)
  x <- stats::runif(70, -1, 1)
  y <- stats::runif(70, -1, 1)
  loc <- data.frame(x = x, y = y)
  spde <- make_mesh(loc, c("x", "y"), cutoff = 0.02)
  dat <- sdmTMB_sim(x = x, y = y, mesh = spde,
    time_steps = 9, rho = 0, range = 0.2,
    sigma_O = 0, sigma_E = 0.3, phi = 0.1,
    seed = 1
  )
  spde <- make_mesh(dat, c("x", "y"), cutoff = 0.02)
  m <- sdmTMB(
    ar1_fields = FALSE, include_spatial = FALSE,
    data = dat, formula = observed ~ 1, time = "time",
    family = gaussian(link = "identity"), spde = spde
  )
  p <- predict(m)
  p_nd <- predict(m, newdata = dat)
  tidy(m)
  tidy(m, conf.int = TRUE)
  tidy(m, effects = "ran_par")
  tidy(m, effects = "ran_par", conf.int = TRUE)

  cols <- c("est", "est_non_rf", "est_rf", "omega_s", "epsilon_st")
  expect_equal(p[,cols], p_nd[,cols], tolerance = 1e-3)

  m <- sdmTMB(
    ar1_fields = TRUE, include_spatial = FALSE,
    data = dat, formula = observed ~ 1, time = "time",
    family = gaussian(link = "identity"), spde = spde
  )

  p <- predict(m)
  p_nd <- predict(m, newdata = dat)
  expect_equal(p[,cols], p_nd[,cols], tolerance = 1e-3)
})

test_that("A time-varying model fits and predicts appropriately", {
  skip_on_cran()
  skip_on_ci()
  SEED <- 42
  set.seed(SEED)
  x <- stats::runif(60, -1, 1)
  y <- stats::runif(60, -1, 1)
  initial_betas <- 0.5
  range <- 0.5
  sigma_O <- 0
  sigma_E <- 0.1
  phi <- 0.1
  sigma_V <- 0.3
  loc <- data.frame(x = x, y = y)
  spde <- make_mesh(loc, c("x", "y"), cutoff = 0.02)

  s <- sdmTMB_sim(
    x = x, y = y, mesh = spde, range = range,
    betas = initial_betas, time_steps = 12L, sigma_V = sigma_V,
    phi = phi, sigma_O = sigma_O, sigma_E = sigma_E,
    seed = SEED
  )
  spde <- make_mesh(s, c("x", "y"), cutoff = 0.02)
  m <- sdmTMB(data = s, formula = observed ~ 0, include_spatial = FALSE,
    time_varying = ~ 0 + cov1, time = "time", spde = spde, mgcv = FALSE)
  expect_equal(exp(m$model$par["ln_tau_V"])[[1]], sigma_V, tolerance = 0.05)
  tidy(m, effects = "ran_par")
  b_t <- dplyr::group_by(s, time) %>%
    dplyr::summarize(b_t = unique(b), .groups = "drop") %>%
    dplyr::pull(b_t)
  r <- m$tmb_obj$report()
  b_t_fit <- r$b_rw_t
  plot(b_t, b_t_fit, asp = 1);abline(a = 0, b = 1)
  expect_equal(mean((b_t- b_t_fit)^2), 0, tolerance = 1e-4)
  p <- predict(m)
  plot(p$est, s$observed, asp = 1);abline(a = 0, b = 1)
  expect_equal(mean((p$est - s$observed)^2), 0, tolerance = 0.01)

  cols <- c("est", "est_non_rf", "est_rf", "omega_s", "epsilon_st")
  p_nd <- predict(m, newdata = s)
  expect_equal(p[,cols], p_nd[,cols], tolerance = 1e-4)
})

test_that("Year indexes get created correctly", {
  expect_identical(make_year_i(c(1, 2, 3)),    c(0L, 1L, 2L))
  expect_identical(make_year_i(c(1L, 2L, 3L)), c(0L, 1L, 2L))
  expect_identical(make_year_i(c(1L, 2L, 4L)), c(0L, 1L, 2L))
  expect_identical(make_year_i(c(1, 2, 4, 2)), c(0L, 1L, 2L, 1L))
  expect_identical(make_year_i(c(1, 4, 2)),    c(0L, 2L, 1L))
})

test_that("A spatial trend model fits", {
  skip_on_cran()
  skip_on_ci()
  d <- subset(pcod, year >= 2011) # subset for speed
  pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 30)
  m <- sdmTMB(density ~ depth_scaled, data = d,
    spde = pcod_spde, family = tweedie(link = "log"),
    spatial_trend = TRUE, time = "year")
  expect_true(all(!is.na(summary(m$sd_report)[,"Std. Error"])))
})

test_that("A logistic threshold model fits", {
  skip_on_cran()
  skip_on_ci()
  d <- subset(pcod, year >= 2011) # subset for speed
  pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 30)
  m <- sdmTMB(density ~ 0 + as.factor(year) + logistic(depth_scaled), data = d,
    spde = pcod_spde, family = tweedie(link = "log"),
    spatial_trend = FALSE, time = "year")
  expect_true(all(!is.na(summary(m$sd_report)[,"Std. Error"])))
})

test_that("A linear threshold model fits", {
  skip_on_cran()
  skip_on_ci()
  d <- subset(pcod, year >= 2011) # subset for speed
  pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 30)
  m <- sdmTMB(density ~ 0 + as.factor(year) + breakpt(depth_scaled), data = d,
    spde = pcod_spde, family = tweedie(link = "log"),
    spatial_trend = FALSE, time = "year")
  expect_true(all(!is.na(summary(m$sd_report)[,"Std. Error"])))
})
