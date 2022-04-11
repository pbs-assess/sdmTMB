test_that("Simulated residuals work", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  fit <- sdmTMB(density ~ as.factor(year) + poly(depth, 2),
    data = pcod_2011, time = "year", mesh = pcod_mesh_2011,
    family = tweedie(link = "log"), spatial = "off",
    spatiotemporal = "off"
  )
  s <- simulate(fit, nsim = 50)
  dharma_residuals(s, fit)
  r <- dharma_residuals(s, fit, plot = FALSE)
  expect_equal(class(r), "data.frame")
  expect_error(dharma_residuals(c(1, 2, 3), fit))
  expect_error(dharma_residuals(matrix(c(1, 2, 3)), fit))
  expect_error(dharma_residuals(s, fit, plot = "test"))
  expect_error(dharma_residuals(s, 99))
})

test_that("randomized quantile residuals work,", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")

  predictor_dat <- data.frame(X = runif(2000), Y = runif(2000))
  mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.2)

  sim_dat <- function(family, phi = 0.2, sigma_O = 0.2) {
    set.seed(1)
    sdmTMB_simulate(
      formula = ~1,
      data = predictor_dat,
      mesh = mesh,
      family = family,
      range = 0.5,
      phi = phi,
      sigma_O = sigma_O,
      seed = 1,
      B = 0
    )
  }

  check_resids <- function(fit) {
    set.seed(123)
    r <- residuals(fit)
    qqnorm(r)
    qqline(r)
    # p <- stats::shapiro.test(r)
    # expect_gt(p$p.value, 0.05)
  }
  check_resids_dharma <- function(fit) {
    set.seed(1)
    dharma_residuals(simulate(fit, nsim = 500), fit)
  }

  d <- sim_dat(gaussian())
  fit <- sdmTMB(
    observed ~ 1,
    family = gaussian(),
    data = d, mesh = mesh
  )
  check_resids(fit)
  check_resids_dharma(fit)

  d <- sim_dat(lognormal())
  fit <- sdmTMB(
    observed ~ 1,
    family = lognormal(),
    data = d, mesh = mesh
  )
  check_resids(fit)
  check_resids_dharma(fit)

  d <- sim_dat(Gamma(link = "log"), phi = 0.3, sigma_O = 0.001)
  fit <- sdmTMB(
    observed ~ 0,
    family = Gamma(link = "log"),
    data = d, mesh = mesh, spatial = "off", spatiotemporal = "off"
  )
  check_resids(fit)
  check_resids_dharma(fit)

  d <- sim_dat(binomial(), sigma_O = 0.001)
  fit <- sdmTMB(
    observed ~ 1,
    family = binomial(),
    data = d, mesh = mesh, spatial = "off", spatiotemporal = "off"
  )
  check_resids(fit)
  check_resids_dharma(fit)

  d <- sim_dat(nbinom2(), sigma_O = 0.001)
  fit <- sdmTMB(
    observed ~ 0,
    family = nbinom2(),
    data = d, mesh = mesh, spatial = "off", spatiotemporal = "off"
  )
  check_resids(fit)
  check_resids_dharma(fit)

  d <- sim_dat(nbinom1(), sigma_O = 0.001)
  fit <- sdmTMB(
    observed ~ 0,
    family = nbinom1(),
    data = d, mesh = mesh, spatial = "off", spatiotemporal = "off"
  )
  check_resids(fit)
  check_resids_dharma(fit)

  d <- sim_dat(Beta(), phi = 5, sigma_O = 0.001)
  fit <- sdmTMB(
    observed ~ 1,
    family = Beta(),
    data = d, mesh = mesh, spatial = 'off', spatiotemporal = 'off'
  )
  check_resids(fit)
  check_resids_dharma(fit)

  d <- sim_dat(poisson())
  fit <- sdmTMB(
    observed ~ 1,
    family = poisson(),
    data = d, mesh = mesh
  )
  check_resids(fit)
  check_resids_dharma(fit)

  d <- sim_dat(student(df = 2))
  fit <- sdmTMB(
    observed ~ 1,
    family = student(df = 2),
    data = d, mesh = mesh
  )
  check_resids(fit)
  check_resids_dharma(fit)

  # wrong df:
  fit <- sdmTMB(
    observed ~ 1,
    family = student(df = 10),
    data = d, mesh = mesh
  )
  r <- residuals(fit)
  qqnorm(r)
  qqline(r)
  p <- stats::shapiro.test(r)
  expect_lt(p$p.value, 0.05) # less than 0.05!
  check_resids_dharma(fit)

  d <- sim_dat(tweedie())
  fit <- sdmTMB(
    observed ~ 1,
    family = tweedie(),
    data = d, mesh = mesh
  )
  check_resids(fit)
  check_resids_dharma(fit)

  d <- sim_dat(truncated_nbinom2())
  fit <- sdmTMB(
    observed ~ 1,
    family = truncated_nbinom2(),
    data = d, mesh = mesh,
    spatial = "off", spatiotemporal = "off"
  )
  expect_error(residuals(fit), regexp = "truncated_nbinom2")
  check_resids_dharma(fit)

  d <- sim_dat(truncated_nbinom1())
  fit <- sdmTMB(
    observed ~ 1,
    family = truncated_nbinom1(),
    data = d, mesh = mesh,
    spatial = "off", spatiotemporal = "off"
  )
  expect_error(residuals(fit), regexp = "truncated_nbinom1")
  check_resids_dharma(fit)
})

test_that("residuals() works", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  pcod_spde <- make_mesh(pcod, c("X", "Y"), cutoff = 15)
  fit <- sdmTMB(density ~ 1, spatial = "off",
    data = pcod, mesh = pcod_spde,
    family = tweedie()
  )
  r <- residuals(fit)
  expect_true(length(r) == nrow(pcod))
  expect_true(sum(is.na(r)) == 0L)

  fit <- sdmTMB(present ~ 1, spatial = "off",
    data = pcod, mesh = pcod_spde,
    family = binomial()
  )
  r <- residuals(fit)
  expect_true(length(r) == nrow(pcod))
  expect_true(sum(is.na(r)) == 0L)

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde,
    family = delta_gamma(), spatial = "off",
    control = sdmTMBcontrol(newton_loops = 1)
  )
  expect_error(r <- residuals(fit), regexp = "delta")

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde,
    family = delta_poisson_link_gamma(), spatial = "off",
    control = sdmTMBcontrol(newton_loops = 1)
  )
  expect_error(r <- residuals(fit), regexp = "delta")

  })
