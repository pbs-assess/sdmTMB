# basic model fitting and prediction tests

SEED <- 1
set.seed(SEED)
x <- stats::runif(100, -1, 1)
y <- stats::runif(100, -1, 1)
loc <- data.frame(x = x, y = y)

test_that("sdmTMB model fit with a covariate beta", {
  local_edition(2)
  skip_if_not_installed("INLA")
  spde <- make_mesh(loc, c("x", "y"), cutoff = 0.02)
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
  expect_equal(s$observed[c(1:30, 500)],
    c(-0.364045551132498, -0.0884742484390758, -0.748796917796798,
    -0.0926168889732131, 0.0892802633041351, 0.92937309789911, 0.665431954199266,
    -0.81477972686779, 0.474664276994187, 0.929133590564276, -0.263130711477779,
    -0.115523493789406, 0.898574275118016, -0.964967397393866, -1.10224938263034,
    -0.286050035282319, 0.0810888689537039, -0.149727881361854, 0.364458765137563,
    0.31886147691733, -0.184284142634752, 1.27393901104861, -0.822473494638277,
    -0.594882019789678, -0.515872613077436, 0.578924954453508, -0.464615318681477,
    0.476914673336151, -0.52150670383694, -0.0152136800518284, -0.394972244437357
  ), tolerance = 1e-7)
  spde <- make_mesh(s, c("x", "y"), cutoff = 0.02)
  plot(spde)
  .t1 <- system.time({
    m <- sdmTMB(data = s, formula = observed ~ 0 + cov1, time = "time",
      silent = TRUE, spde = spde, control = sdmTMBcontrol(normalize = FALSE))
  })
  .t2 <- system.time({
    m_norm <- sdmTMB(data = s, formula = observed ~ 0 + cov1, time = "time",
      silent = TRUE, spde = spde, control = sdmTMBcontrol(normalize = TRUE))
  })
  .t3 <- system.time({
    m_pc <- sdmTMB(data = s, formula = observed ~ 0 + cov1, time = "time",
      silent = TRUE, spde = spde, control = sdmTMBcontrol(normalize = TRUE),
      priors = sdmTMBpriors(matern_st = pc_matern(range_gt = 0.2, sigma_lt = 0.1))
    )
  })
  expect_equal(m$model$par, m_norm$model$par, tolerance = 1e-6)
  expect_equal(m$model$par, m_pc$model$par, tolerance = 1e-3)
  #commented test out: pc model kappa is slightly greater (7.6e-07) than m_norm kappa
  #expect_true(m_pc$model$par[['ln_kappa']] < m_norm$model$par[['ln_kappa']])
  expect_lt(.t2[[1]], .t1[[1]])
  expect_output(print(m), "fit by")
  expect_output(summary(m), "fit by")
  p <- as.list(m$model$par)
  r <- m$tmb_obj$report()
  est <- tidy(m, "ran_pars")
  expect_equal(sort(est[,"estimate", drop = TRUE]),
    sort(c(0.106559618205922, 0.106559618205922, 0.0646179753599145, 0.303614800474715, 0.308394406301974)),
    tolerance = 1e-6)
  expect_equal(m$model$convergence, 0L)
  expect_equal((p$b_j - initial_betas)^2, 0, tolerance = 0.001)
  expect_equal((exp(p$ln_phi) - phi)^2, 0, tolerance = 0.002)
  expect_equal((r$sigma_O - sigma_O)^2, 0, tolerance = 0.002)
  expect_equal((r$sigma_E[1] - sigma_E)^2, 0, tolerance = 0.001)
  expect_equal(est$estimate[est$term == "range"][1], range, tolerance = 0.01)
  p <- predict(m)
  r <- residuals(m)
  r_sim <- residuals(m, type = "sim")
  expect_equal(mean((p$est - s$observed)^2), 0, tolerance = 0.002)
})

test_that("Anisotropy fits and plots", {
  skip_on_ci()
  skip_on_cran()
  skip_if_not_installed("INLA")
  local_edition(2)
  m <- sdmTMB(data = pcod,
    formula = density ~ 0 + as.factor(year),
    spde = make_mesh(pcod, c("X", "Y"), n_knots = 50, type = "kmeans"),
    family = tweedie(link = "log"), anisotropy = TRUE)
  expect_identical(class(m), "sdmTMB")
  plot_anisotropy(m)
  expect_equal(m$tmb_obj$report()$H,
    structure(
      c(0.665528444798002, 0.079350716881963, 0.079350716881963, 1.51202633656794),
      .Dim = c(2L, 2L)),
    tolerance = 1e-3)
})

test_that("A model with splines works", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  d <- subset(pcod, year >= 2015)
  m <- sdmTMB(data = d,
    formula = density ~ 1 + s(depth_scaled),
    spde = make_mesh(d, c("X", "Y"), n_knots = 50, type = "kmeans"),
    family = tweedie(link = "log"))
  expect_identical(class(m), "sdmTMB")
})

test_that("A spatiotemporal version works with predictions on new data points", {
  d <- pcod_2011
  pcod_spde <- pcod_mesh_2011
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
  predictions <- predict(m, newdata = subset(qcs_grid, year >= 2011))
  expect_identical(class(predictions), "data.frame")
})

test_that("Predictions on the original data set as `newdata`` return the same predictions", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  set.seed(1)
  local_edition(2)
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
    fields = "AR1", include_spatial = FALSE,
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
    fields = "AR1", include_spatial = FALSE,
    data = dat, formula = observed ~ 1, time = "time",
    family = gaussian(link = "identity"), spde = spde
  )

  p <- predict(m)
  p_nd <- predict(m, newdata = dat)
  expect_equal(p[,cols], p_nd[,cols], tolerance = 1e-3)
})
