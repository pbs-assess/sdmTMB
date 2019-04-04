context("basic model fitting and prediction tests")


SEED <- 97999

x <- stats::runif(201, 0, 10)
y <- stats::runif(201, 0, 10)

# ------------------------------------------------------
# a sdmTMB model with one covariate

test_that("sdmTMB model fit with a covariate beta", {
  skip_on_cran()
  #skip_on_travis()

  set.seed(SEED)

  initial_betas <- c(-0.5)
  kappa <- 0.5 # decay of spatial correlation (smaller = slower decay)
  sigma_O <- 0.3 # SD of spatial process
  phi <- 0.1 # observation error

  s <- sim(
    x = x, y = y,
    initial_betas = initial_betas,
    #  sigma_V = sigma_V, time_steps = time_steps, ar1_fields = ar1_fields, ar1_phi = ar1_phi, sigma_E = sigma_E,
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

  expect_equal(mean((p$data$est - s$observed)^2), 0, tol = 0.05)

  })
