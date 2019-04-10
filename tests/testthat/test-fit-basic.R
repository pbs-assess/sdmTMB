context("basic model fitting and prediction tests")

SEED <- 1
set.seed(SEED)
x <- stats::runif(150, -1, 1)
y <- stats::runif(150, -1, 1)

test_that("sdmTMB model fit with a covariate beta", {
  initial_betas <- 0.5
  kappa <- 4 # decay of spatial correlation (smaller = slower decay)
  sigma_O <- 0.3 # SD of spatial process
  sigma_E <- 0.3 # SD of spatial process
  phi <- 0.1 # observation error
  s <- sim(
    x = x, y = y,
    initial_betas = initial_betas, time_steps = 9L,
    phi = phi, kappa = kappa, sigma_O = sigma_O, sigma_E = sigma_E,
    seed = SEED, plot = TRUE
  )
  spde <- make_spde(x = s$x, y = s$y, n_knots = 100)
  plot_spde(spde)
  m <- sdmTMB(data = s, formula = observed ~ 0 + cov1, time = "time", spde = spde)
  p <- as.list(m$model$par)
  r <- m$tmb_obj$report()
  expect_equal(m$model$convergence, 0L)
  expect_equal((p$b_j - initial_betas)^2, 0, tol = 0.05)
  expect_equal((exp(p$ln_phi) - phi)^2, 0, tol = 0.05)
  expect_equal((r$sigma_O - sigma_O)^2, 0, tol = 0.05)
  expect_equal((r$sigma_E - sigma_E)^2, 0, tol = 0.05)
  expect_equal(exp(p$ln_kappa), kappa, tol = 1.1)
  p <- predict(m)
  r <- residuals(m)
  expect_equal(mean((p$data$est - s$observed)^2), 0, tol = 0.02)
})

test_that("Anisotropy fits and plots", {
  d <- subset(pcod, year >= 2015)
  m <- sdmTMB(data = d,
    formula = density ~ 0 + as.factor(year),
    spde = make_spde(d$X, d$Y, n_knots = 40),
    family = tweedie(link = "log"), anisotropy = TRUE)
  expect_identical(class(m), "sdmTMB")
  g <- plot_anisotropy(m)
  expect_identical(class(g), c("gg", "ggplot"))
})

test_that("A spatiotemporal version works with predictions on new data points", {
  d <- subset(pcod, year >= 2015)
  pcod_spde <- make_spde(d$X, d$Y, n_knots = 40)
  m <- sdmTMB(
    data = d,
    formula = density ~ 0 + as.factor(year),
    time = "year", spde = pcod_spde, family = tweedie(link = "log"),
    include_spatial = FALSE
  )
  # Predictions at original data locations:
  predictions <- predict(m)$data
  predictions$resids <- residuals(m) # randomized quantile residuals
  # Predictions onto new data:
  predictions <- predict(m, newdata = subset(qcs_grid, year >= 2015))
  expect_identical(class(predictions$data), "data.frame")
})

test_that("AR1 models fit with and without R normalization", {
  set.seed(1)
  x <- stats::runif(70, -1, 1)
  y <- stats::runif(70, -1, 1)
  dat <- sim(x = x, y = y,
    time_steps = 9, ar1_fields = TRUE, ar1_phi = 0.0,
    plot = TRUE, sigma_O = 1e-6, sigma_E = 0.3, phi = 0.1,
    seed = 1
  )
  spde <- make_spde(x = dat$x, y = dat$y, n_knots = 40)
  m <- sdmTMB(
    ar1_fields = TRUE, include_spatial = FALSE,
    data = dat, formula = observed ~ 1, time = "time",
    family = gaussian(link = "identity"), spde = spde, normalize = FALSE
  )
  m_normalize <- sdmTMB(
    ar1_fields = TRUE, include_spatial = FALSE,
    data = dat, formula = observed ~ 1, time = "time",
    family = gaussian(link = "identity"), spde = spde, normalize = TRUE
  )
  expect_equal(m$model$objective, m_normalize$model$objective, tolerance = 1e-5)
  expect_equal(m$model$par, m_normalize$model$par, tolerance = 1e-5)

  p <- predict(m)
  p_normalize <- predict(m_normalize)
  expect_equal(p$data, p_normalize$data)

  p_nd <- predict(m, newdata = dat, xy_cols = c("x", "y"))
  p_normalize <- predict(m_normalize, newdata = dat, xy_cols = c("x", "y"))
  expect_equal(p_nd$data, p_normalize$data, tolerance = 1e-4)
})

test_that("Predictions on the original data set as `newdata`` return the same predictions", {
  set.seed(1)
  x <- stats::runif(70, -1, 1)
  y <- stats::runif(70, -1, 1)
  dat <- sim(x = x, y = y,
    time_steps = 9, ar1_fields = TRUE, ar1_phi = 0.0,
    plot = TRUE, sigma_O = 1e-6, sigma_E = 0.3, phi = 0.1,
    seed = 1
  )
  spde <- make_spde(x = dat$x, y = dat$y, n_knots = 40)
  m <- sdmTMB(
    ar1_fields = FALSE, include_spatial = FALSE,
    data = dat, formula = observed ~ 1, time = "time",
    family = gaussian(link = "identity"), spde = spde, normalize = FALSE
  )
  p <- predict(m)
  p_nd <- predict(m, newdata = dat, xy_cols = c("x", "y"))

  cols <- c("est", "est_re_s", "est_re_st", "est_re_s_trend")
  expect_equal(p$data[,cols], p_nd$data[,cols], tolerance = 1e-4)

  m <- sdmTMB(
    ar1_fields = TRUE, include_spatial = FALSE,
    data = dat, formula = observed ~ 1, time = "time",
    family = gaussian(link = "identity"), spde = spde, normalize = FALSE
  )

  p <- predict(m)
  p_nd <- predict(m, newdata = dat, xy_cols = c("x", "y"))
  expect_equal(p$data[,cols], p_nd$data[,cols], tolerance = 1e-4)
})
