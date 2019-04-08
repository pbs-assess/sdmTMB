context("basic model fitting and prediction tests")


SEED <- 97999

x <- stats::runif(201, 0, 10)
y <- stats::runif(201, 0, 10)

# ------------------------------------------------------
# a sdmTMB model with one covariate

test_that("sdmTMB model fit with a covariate beta", {
  set.seed(SEED)
  initial_betas <- c(-0.5)
  kappa <- 0.5 # decay of spatial correlation (smaller = slower decay)
  sigma_O <- 0.3 # SD of spatial process
  phi <- 0.1 # observation error
  s <- sim(
    x = x, y = y,
    initial_betas = initial_betas,
    phi = phi, kappa = kappa, sigma_O = sigma_O, seed = SEED #, plot = TRUE
  )

  spde <- make_spde(x = s$x, y = s$y, n_knots = 200)
  plot_spde(spde)

  m <- sdmTMB(
    data = s, formula = observed ~ 0 + cov1,
    spde = spde
  )
  r <- m$tmb_obj$report()
  expect_equal(m$model$convergence, 0L)
  expect_equal((r$b_j - initial_betas)^2, 0, tol = 0.05)
  expect_equal((exp(r$ln_phi) - phi)^2, 0, tol = 0.05)
  # expect_equal((r$sigma_O - sigma_O)^2, 0, tol = 0.05)
  # FIXME: does it make sense that kappa cannot be predicted well in model with one time slice?
  # expect_equal((exp(r$ln_kappa) - kappa)^2, 0, tol = 0.05)
  m$model$par
  m$tmb_obj$gr(m$model$par)
  p <- predict(m)
  r <- residuals(m)
  expect_equal(mean((p$data$est - s$observed)^2), 0, tol = 0.05)
})

test_that("NB2 fits", {
  d <- pcod[pcod$year == 2017, ]
  d$density <- round(d$density)
  spde <- make_spde(d$X, d$Y, n_knots = 45)
  m <- sdmTMB(data = d, formula = density ~ 1,
    spde = spde, family = nbinom2(link = "log"))
  sdmTMBphi <- exp(m$model$par[["ln_phi"]])
  m2 <- glmmTMB::glmmTMB(density ~ 1, data = d,
    family = glmmTMB::nbinom2(link = "log"))
  glmmTMBphi <- exp(m2$fit$par[["betad"]])
  expect_equal(glmmTMBphi, sdmTMBphi, tol = 0.01)
})

test_that("Anisotropy fits and plots", {
  m <- sdmTMB(data = subset(pcod, year >= 2015),
    formula = density ~ 0 + as.factor(year),
    time = "year", spde = make_spde(pcod$X, pcod$Y, n_knots = 40),
    family = tweedie(link = "log"), anisotropy = TRUE,
    include_spatial = FALSE)
  expect_identical(class(m), "sdmTMB")
  g <- plot_anisotropy(m)
  expect_identical(class(g), c("gg", "ggplot"))
})

test_that("Families work", {
  .names <- c("family", "link", "linkfun", "linkinv")
  expect_named(student(link = "identity"), .names)
  expect_named(lognormal(link = "log"), .names)
  expect_named(tweedie(link = "log"), .names)
  expect_named(nbinom2(link = "log"), .names)
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
