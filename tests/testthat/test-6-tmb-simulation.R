test_that("TMB IID simulation works", {
  skip_on_ci()
  skip_on_cran()
  skip_if_not_installed("INLA")

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
  skip_if_not_installed("INLA")

  set.seed(1)
  predictor_dat <- data.frame(
    X = runif(2000), Y = runif(2000),
    a1 = rnorm(2000), year = rep(1:10, each = 200)
  )
  mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)
  sim_dat <- sdmTMB_simulate(
    formula = ~1,
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
    B = 0 # B0 = intercept, B1 = a1 slope
  )
  fit <- sdmTMB(observed ~ 0, sim_dat,
    mesh = mesh, time = "year",
    spatiotemporal = "ar1", spatial = "off",
  )
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
  skip_if_not_installed("INLA")

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


test_that("TMB (custom) AR1 simulation is unbiased", {
  # run many times; check for bias
  skip_on_ci()
  skip_on_cran()
  skip_if_not_installed("INLA")

  do_sim_fit <- function(i) {
    cat("Iteration", i, "\n")
    set.seed(i)
    predictor_dat <- data.frame(
      X = runif(1000), Y = runif(1000),
      a1 = rnorm(1000), year = rep(1:10, each = 100)
    )
    mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)
    sim_dat <- sdmTMB_simulate(
      formula = ~1,
      data = predictor_dat,
      time = "year",
      mesh = mesh,
      family = gaussian(),
      range = 0.5,
      sigma_E = 0.1,
      phi = 0.1,
      sigma_O = 0,
      seed = i,
      rho = 0.8,
      B = 0
    )
    fit <- sdmTMB(observed ~ 0, sim_dat,
      mesh = mesh, time = "year",
      spatiotemporal = "ar1", spatial = "off",
    )
    b <- tidy(fit, "ran_pars")
    rho_hat <- b$estimate[b$term == "rho"]
    range_hat <- b$estimate[b$term == "range"]
    sigma_E_hat <- b$estimate[b$term == "sigma_E"]
    data.frame(rho = rho_hat, sigma_E = sigma_E_hat, range = range_hat)
  }
  out <- lapply(seq_len(12L), function(i) do_sim_fit(i))
  out <- do.call("rbind", out)
  expect_true(median(out$rho) > 0.79 && median(out$rho) < 0.81)
  expect_true(median(out$range) > 0.48 && median(out$range) < 0.52)
  expect_true(median(out$sigma_E) > 0.09 && median(out$sigma_E) < 0.11)
})


# test_that("TMB Delta simulation works", {
#   skip_on_ci()
#   skip_on_cran()
#   skip_if_not_installed("INLA")
#
#   set.seed(1)
#   predictor_dat <- data.frame(
#     X = runif(2000), Y = runif(2000),
#     a1 = rnorm(2000), year = 1
#   )
#   mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)
#   sim_dat <- sdmTMB_simulate(
#     formula = ~1,
#     data = predictor_dat,
#     time = "year",
#     mesh = mesh,
#     family = nbinom2(),
#     range = 0.3,
#     phi = 1,
#     sigma_O = 0.1,
#     seed = 1,
#     B = 1.5
#   )
#   sim_dat$observed[sample(1:nrow(sim_dat), size = 0.7*nrow(sim_dat), replace=T)] = 0
#   fit <- sdmTMB(observed ~ 1, sim_dat,
#                 mesh = mesh,
#                 family = delta_truncated_nbinom2(),
#                 share_range = TRUE
#   )
#
#   sim_1 <- simulate(fit, model = 1, nsim = 100)
#   expect_equal(ncol(sim_1), 100)
#   empirical_z <- 1 - length(which(sim_dat$observed==0)) / nrow(sim_dat)
#   expect_lt(max(abs(apply(sim_1,2,mean) - empirical_z)), 0.05)
#
#   sim_2 <- simulate(fit, model = 2, nsim = 100)
#   expect_equal(ncol(sim_2), 100)
#   expect_lt(max(abs(apply(sim_2,2,mean,na.rm=T) - 5.67)), 0.65)
# })
