test_that("TMB IID simulation works", {
  skip_on_ci()
  skip_on_cran()

  set.seed(1)
  predictor_dat <- data.frame(
    X = runif(2000), Y = runif(2000),
    a1 = rnorm(2000), year = rep(1:10, each = 200)
  )
  mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)

  sim_dat <- sdmTMB_simulate(
    formula = ~ 1 + a1,
    data = predictor_dat,
    time = "year",
    mesh = mesh,
    family = gaussian(),
    range = 0.5,
    sigma_E = 0.1,
    phi = 0.1,
    sigma_O = 0.2,
    seed = 42,
    B = c(0.2, -0.4) # B0 = intercept, B1 = a1 slope
  )
  fit <- sdmTMB(observed ~ a1, sim_dat, mesh = mesh, time = "year")
  b <- tidy(fit)
  b
  expect_equal(b$estimate[b$term == "a1"], -0.4, tolerance = 0.1)
  expect_equal(b$estimate[b$term == "(Intercept)"], 0.2, tolerance = 0.2)
  b <- tidy(fit, "ran_pars")
  b
})

test_that("TMB AR1 simulation works", {
  skip_on_ci()
  skip_on_cran()

  set.seed(1)
  predictor_dat <- data.frame(
    X = runif(2000), Y = runif(2000),
    a1 = rnorm(2000), year = rep(1:10, each = 200)
  )
  mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)
  sim_dat <- sdmTMB_simulate(
    formula = ~ 1 + a1,
    data = predictor_dat,
    time = "year",
    mesh = mesh,
    family = gaussian(),
    range = 0.5,
    sigma_E = 0.1,
    phi = 0.1,
    sigma_O = 0,
    seed = 42,
    rho = 0.8, #<
    B = c(0.2, -0.4) # B0 = intercept, B1 = a1 slope
  )
  fit <- sdmTMB(observed ~ a1, sim_dat,
    mesh = mesh, time = "year",
    spatiotemporal = "ar1", spatial = "off"
  )
  b <- tidy(fit)
  b
  b <- tidy(fit, "ran_pars")
  b
  rho_hat <- b$estimate[b$term == "rho"]
  expect_true(rho_hat > 0.7 && rho_hat < 0.9)
  sigma_E_hat <- b$estimate[b$term == "sigma_E"]
  expect_true(sigma_E_hat > 0.07 && sigma_E_hat < 0.13)
})

test_that("TMB RW simulation works", {
  skip_on_ci()
  skip_on_cran()

  set.seed(1)
  predictor_dat <- data.frame(
    X = runif(2000), Y = runif(2000),
    a1 = rnorm(2000), year = rep(1:20, each = 200)
  )
  mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)

  sim_dat <- sdmTMB_simulate(
    formula = ~1,
    data = predictor_dat,
    time = "year",
    mesh = mesh,
    family = gaussian(),
    range = 0.3,
    sigma_E = 0.1,
    phi = 0.05,
    sigma_O = 0,
    seed = 42,
    rho = 1, #<
    B = 0
  )
  fit_ar1 <- sdmTMB(observed ~ 0,
    sim_dat,
    mesh = mesh, time = "year",
    spatiotemporal = "ar1", #<
    spatial = "off"
  )
  b_ar1 <- tidy(fit_ar1, "ran_pars")
  b_ar1
  rho_hat <- b_ar1$estimate[b_ar1$term == "rho"]
  expect_true(rho_hat > 0.9)
  sigma_E_hat_ar1 <- b_ar1$estimate[b_ar1$term == "sigma_E"]

  fit_rw <- sdmTMB(observed ~ 0,
    sim_dat,
    mesh = mesh, time = "year",
    spatiotemporal = "rw", #<
    spatial = "off"
  )
  b_rw <- tidy(fit_rw, "ran_pars")
  b_rw
  sigma_E_hat_rw <- b_rw$estimate[b_rw$term == "sigma_E"]
  sigma_E_hat_rw
  expect_true(sigma_E_hat_rw > 0.07 && sigma_E_hat_rw < 1.13)
  expect_true(sigma_E_hat_rw < sigma_E_hat_ar1)
})

test_that("TMB breakpt sims work", {
  skip_on_ci()
  skip_on_cran()

  set.seed(1)
  predictor_dat <- data.frame(
    X = runif(1000), Y = runif(1000),
    a1 = rnorm(1000)
  )
  mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.2)
  sim_dat <- sdmTMB_simulate(
    formula = ~ 1 + breakpt(a1),
    data = predictor_dat,
    mesh = mesh,
    family = gaussian(),
    range = 0.5,
    phi = 0.02,
    sigma_O = 0.001,
    seed = 42,
    B = 0,
    threshold_coefs = c(0.3, 0)
  )
  plot(predictor_dat$a1, sim_dat$observed)
  expect_lt(max(sim_dat$observed), 0.1)
  sim_dat$a1 <- predictor_dat$a1
  fit <- sdmTMB(observed ~ 1 + breakpt(a1), sim_dat, mesh = mesh, spatial = "off",
    control = sdmTMBcontrol(newton_loops = 1L))
  b <- tidy(fit)
  expect_gt(b$estimate[b$term == "a1-breakpt"], -0.02)
  expect_lt(b$estimate[b$term == "a1-breakpt"], 0.02)
  expect_gt(b$estimate[b$term == "a1-slope"], 0.28)
  expect_lt(b$estimate[b$term == "a1-slope"], 0.32)

  sim_dat <- sdmTMB_simulate(
    formula = ~ 1 + logistic(a1),
    data = predictor_dat,
    mesh = mesh,
    family = gaussian(),
    range = 0.5,
    phi = 0.001,
    sigma_O = 0.001,
    seed = 42,
    B = 0,
    threshold_coefs = c(0.2, 0.4, 0.5)
  )
  plot(predictor_dat$a1, sim_dat$observed)
  expect_lt(max(sim_dat$observed), 0.53)
  expect_gt(min(sim_dat$observed), -0.05)
})
