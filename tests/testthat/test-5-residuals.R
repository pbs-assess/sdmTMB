# test_that("Simulated residuals work", {
#   skip_on_cran()
#   skip_on_ci()
#   skip_if_not_installed("INLA")
#   fit <- sdmTMB(density ~ as.factor(year) + poly(depth, 2),
#     data = pcod_2011, time = "year", mesh = pcod_mesh_2011,
#     family = tweedie(link = "log"), spatial = "off",
#     spatiotemporal = "off"
#   )
#   s <- simulate(fit, nsim = 50)
  # dharma_residuals(s, fit)
  # r <- dharma_residuals(s, fit, plot = FALSE)
  # expect_equal(class(r), "data.frame")
  # expect_error(dharma_residuals(c(1, 2, 3), fit))
  # expect_error(dharma_residuals(matrix(c(1, 2, 3)), fit))
  # expect_error(dharma_residuals(s, fit, plot = "test"))
  # expect_error(dharma_residuals(s, 99))
# })

test_that("randomized quantile residuals work,", {
  skip_on_cran()
  skip_on_ci()

  set.seed(1)
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
    p <- stats::shapiro.test(r)
    expect_gt(p$p.value, 0.01)
    invisible(r)
  }
  # check_resids_dharma <- function(fit) {
  #   set.seed(1)
  #   dharma_residuals(simulate(fit, nsim = 500), fit)
  # }

  d <- sim_dat(gaussian())
  fit <- sdmTMB(
    observed ~ 1,
    family = gaussian(),
    data = d, mesh = mesh
  )
  check_resids(fit)
  check_resids(fit)
  # check_resids_dharma(fit)

  d <- sim_dat(lognormal())
  fit <- sdmTMB(
    observed ~ 1,
    family = lognormal(),
    data = d, mesh = mesh
  )
  check_resids(fit)
  # check_resids_dharma(fit)

  d <- sim_dat(Gamma(link = "log"), phi = 0.3, sigma_O = 0.001)
  fit <- sdmTMB(
    observed ~ 0,
    family = Gamma(link = "log"),
    data = d, mesh = mesh, spatial = "off", spatiotemporal = "off"
  )
  check_resids(fit)
  # check_resids_dharma(fit)

  d <- sim_dat(binomial(), sigma_O = 0.001)
  fit <- sdmTMB(
    observed ~ 1,
    family = binomial(),
    data = d, mesh = mesh, spatial = "off", spatiotemporal = "off"
  )
  check_resids(fit)
  # check_resids_dharma(fit)

  d <- sim_dat(nbinom2(), sigma_O = 0.001)
  fit <- sdmTMB(
    observed ~ 0,
    family = nbinom2(),
    data = d, mesh = mesh, spatial = "off", spatiotemporal = "off"
  )
  check_resids(fit)
  # check_resids_dharma(fit)

  d <- sim_dat(nbinom1(), sigma_O = 0.001)
  fit <- sdmTMB(
    observed ~ 0,
    family = nbinom1(),
    data = d, mesh = mesh, spatial = "off", spatiotemporal = "off"
  )
  check_resids(fit)
  # check_resids_dharma(fit)

  d <- sim_dat(Beta(), phi = 5, sigma_O = 0.001)
  fit <- sdmTMB(
    observed ~ 1,
    family = Beta(),
    data = d, mesh = mesh, spatial = 'off', spatiotemporal = 'off'
  )
  check_resids(fit)
  # check_resids_dharma(fit)

  set.seed(1)
  d <- sim_dat(poisson())
  fit <- sdmTMB(
    observed ~ 1,
    family = poisson(),
    data = d, mesh = mesh
  )
  check_resids(fit)
  # check_resids_dharma(fit)

  d <- sim_dat(student(df = 2))
  fit <- sdmTMB(
    observed ~ 1,
    family = student(df = 2),
    data = d, mesh = mesh
  )
  check_resids(fit)
  # check_resids_dharma(fit)

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
  # check_resids_dharma(fit)

  d <- sim_dat(tweedie())
  fit <- sdmTMB(
    observed ~ 1,
    family = tweedie(),
    data = d, mesh = mesh
  )
  check_resids(fit)
  # check_resids_dharma(fit)

  d <- sim_dat(truncated_nbinom2())
  fit <- sdmTMB(
    observed ~ 1,
    family = truncated_nbinom2(),
    data = d, mesh = mesh,
    spatial = "off", spatiotemporal = "off"
  )
  expect_error(residuals(fit), regexp = "truncated_nbinom2")
  # check_resids_dharma(fit)

  d <- sim_dat(truncated_nbinom1())
  fit <- sdmTMB(
    observed ~ 1,
    family = truncated_nbinom1(),
    data = d, mesh = mesh,
    spatial = "off", spatiotemporal = "off"
  )
  expect_error(residuals(fit), regexp = "truncated_nbinom1")
  # check_resids_dharma(fit)
})

test_that("residuals() works", {
  skip_on_cran()
  pcod_spde <- make_mesh(pcod, c("X", "Y"), cutoff = 15)
  # fit <- sdmTMB(density ~ 1, spatial = "off",
  #   data = pcod, mesh = pcod_spde,
  #   family = tweedie()
  # )
  # NaNs on some systems!?
  # r <- residuals(fit)
  # expect_true(length(r) == nrow(pcod))
  # expect_true(sum(is.na(r)) == 0L)

  fit <- sdmTMB(present ~ 1, spatial = "off",
    data = pcod, mesh = pcod_spde,
    family = binomial()
  )
  r <- residuals(fit)
  expect_true(length(r) == nrow(pcod))
  expect_true(sum(is.na(r)) == 0L)

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde,
    family = delta_gamma(), spatial = "on",
    control = sdmTMBcontrol(newton_loops = 1)
  )
  r <- residuals(fit)
  r <- residuals(fit, model = 1)
  qqnorm(r)
  set.seed(1)
  r <- residuals(fit, model = 2)
  qqnorm(r[!is.na(r)])

  # matches the Gamma positive-only model:
  pos <- subset(pcod, density > 0)
  mesh <- make_mesh(pos, c("X", "Y"), mesh = pcod_spde$mesh)
  fit2 <- sdmTMB(density ~ 1,
    data = pos, mesh = mesh,
    family = Gamma(link = "log"), spatial = "on",
    control = sdmTMBcontrol(newton_loops = 1)
  )
  set.seed(1)
  rpos <- residuals(fit2)
  expect_equal(as.double(r[!is.na(r)]), rpos)

  })

test_that("Pearson residuals work", {
  # binomial proportion
  set.seed(1)
  w <- sample(1:9, size = 300, replace = TRUE)
  dat <- data.frame(y = stats::rbinom(300, size = w, 0.5))
  dat$prop <- dat$y / w
  m <- sdmTMB(prop ~ 1, data = dat, weights = w, family = binomial(), spatial = "off")
  m1 <- glmmTMB::glmmTMB(prop ~ 1, data = dat, weights = w, family = binomial())
  r <- residuals(m, type = "pearson")
  r1 <- residuals(m1, type = "pearson")
  expect_equal(as.numeric(r), as.numeric(r1))

  # binomial cbind
  dat$y0 <- w - dat$y
  m <- sdmTMB(cbind(y, y0) ~ 1, data = dat, family = binomial(), spatial = "off")
  m1 <- glmmTMB::glmmTMB(cbind(y, y0) ~ 1, data = dat, family = binomial())
  r <- residuals(m, type = "pearson")
  r1 <- residuals(m1, type = "pearson")
  expect_equal(as.numeric(r), as.numeric(r1))

  # gaussian
  m <- sdmTMB(prop ~ 1, data = dat, family = gaussian(), spatial = "off")
  m1 <- glmmTMB::glmmTMB(prop ~ 1, data = dat, family = gaussian())
  r <- residuals(m, type = "pearson")
  r1 <- residuals(m1, type = "pearson")
  expect_equal(as.numeric(r), as.numeric(r1))

  # gamma
  set.seed(1)
  dat$y <- rlnorm(300, 0.4, 0.3)
  m <- sdmTMB(y ~ 1, data = dat, family = Gamma(link = "log"), spatial = "off")
  m1 <- glmmTMB::glmmTMB(y ~ 1, data = dat, family = Gamma(link = "log"))
  r <- residuals(m, type = "pearson")
  r1 <- residuals(m1, type = "pearson")
  expect_equal(as.numeric(r), as.numeric(r1))

  # poisson
  set.seed(1)
  dat$y <- rpois(300, 0.4)
  m <- sdmTMB(y ~ 1, data = dat, family = poisson(link = "log"), spatial = "off")
  m1 <- glmmTMB::glmmTMB(y ~ 1, data = dat, family = poisson(link = "log"))
  r <- residuals(m, type = "pearson")
  r1 <- residuals(m1, type = "pearson")
  expect_equal(as.numeric(r), as.numeric(r1))

  # FIXME: add variance functions; requires fancy environment Theta passing
  set.seed(1)
  dat$y <- rnbinom(300, 0.4, 0.3)
  m <- sdmTMB(y ~ 1, data = dat, family = nbinom2(), spatial = "off")
  # m1 <- glmmTMB::glmmTMB(y ~ 1, data = dat, family = glmmTMB::nbinom2())
  expect_error(r <- residuals(m, type = "pearson"), regexp = "Variance")
  # r1 <- residuals(m1, type = "pearson")
  # expect_equal(as.numeric(r), as.numeric(r1))
})

test_that("MCMC residuals throw error as needed", {
  skip_on_cran()
  fit_dg <- sdmTMB(
    density ~ 1,
    data = pcod_2011,
    mesh = pcod_mesh_2011,
    family = delta_gamma()
  )
  expect_error({
    r <- residuals(fit_dg, type = "mle-mcmc", mcmc_iter = 101, mcmc_warmup = 100)
  }, regexp = "mcmc")
})

test_that("MCMC residuals work with sdmTMBextra", {
  skip_on_cran()
  skip_if_not_installed("sdmTMBextra")
  skip_if_not_installed("rstan")
  skip_on_ci()
  d <- pcod_2011
  set.seed(1)
  d$offset <- rnorm(nrow(d))
  fit_dg <- sdmTMB(
    density ~ 1,
    data = d,
    mesh = pcod_mesh_2011,
    offset = d$offset,
    family = delta_gamma()
  )
  set.seed(1)
  p0 <- sdmTMBextra::predict_mle_mcmc(fit_dg, mcmc_iter = 11, mcmc_warmup = 10)
  set.seed(1)
  p1 <- sdmTMBextra::predict_mle_mcmc(fit_dg,
    mcmc_iter = 11, mcmc_warmup = 10,
    model = 1
  )
  set.seed(1)
  p2 <- sdmTMBextra::predict_mle_mcmc(fit_dg,
    mcmc_iter = 11, mcmc_warmup = 10,
    model = 2
  )
  expect_equal(p0, p1)
  set.seed(1)
  r1 <- residuals(fit_dg, type = "mle-mcmc", mcmc_samples = p1)
  set.seed(1)
  r2 <- residuals(fit_dg, type = "mle-mcmc", mcmc_samples = p2)
  qqnorm(r1)
  qqnorm(r2)
  expect_false(identical(r1, r2))
  expect_true(max(r1) < 5)
  expect_true(min(r1) > -5)
  expect_true(max(r2) < 5)
  expect_true(min(r2) > -5)
})

test_that("predict_mle_mcmc() works with extra_time #297", {
  skip_on_cran()
  skip_if_not_installed("sdmTMBextra")
  skip_if_not_installed("rstan")
  skip_on_ci()
  fit <- sdmTMB(
    density ~ 1,
    time = "year",
    spatiotemporal = "rw",
    spatial = "on",
    mesh = pcod_mesh_2011,
    data = pcod_2011,
    family = tweedie(link = "log"),
    extra_time = c(2012, 2014, 2016)
  )
  set.seed(123)
  samps <- sdmTMBextra::predict_mle_mcmc(fit, mcmc_iter = 11, mcmc_warmup = 10)
  expect_true(nrow(samps) == 972L)
})
