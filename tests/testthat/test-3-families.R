test_that("Families return a name to list with the correct names", {
  .names <- c("family", "link", "linkfun", "linkinv")
  expect_true(all(.names %in% names(student(link = "identity"))))
  expect_true(all(.names %in% names(lognormal(link = "log"))))
  expect_true(all(.names %in% names(tweedie(link = "log"))))
  expect_true(all(.names %in% names(nbinom2(link = "log"))))
})

test_that("The supplementary families work with appropriate links", {
  expect_identical(class(tweedie(link = "log")), "family")
  expect_identical(class(tweedie(link = log)), "family")
  expect_error(class(tweedie(link = "banana")))
  expect_error(class(tweedie(link = banana)))

  expect_identical(class(lognormal(link = "log")), "list")
  expect_identical(class(lognormal(link = log)), "list")
  expect_error(class(lognormal(link = "banana")))
  expect_error(class(lognormal(link = banana)))

  expect_identical(class(nbinom2(link = "log")), "list")
  expect_identical(class(nbinom2(link = log)), "list")
  expect_error(class(nbinom2(link = "banana")))
  expect_error(class(nbinom2(link = inverse)))

  expect_identical(class(student(link = "identity")), "list")
  expect_identical(class(student(link = identity)), "list")
  expect_error(class(student(link = "banana")))
  expect_error(class(student(link = banana)))
})

set.seed(1)
x <- stats::runif(100, -1, 1)
y <- stats::runif(100, -1, 1)
loc <- data.frame(x = x, y = y)

if (suppressWarnings(require("INLA", quietly = TRUE))) {
  spde <- make_mesh(loc, c("x", "y"), n_knots = 50, type = "kmeans")

  test_that("Student family fits", {
    skip_on_ci()
    skip_on_cran()
    set.seed(3)
    initial_betas <- 0.5
    range <- 0.5
    sigma_O <- 0.3
    phi <- 0.01
    s <- sdmTMB_sim(x = x, y = y,
      betas = initial_betas, time = 1L,
      phi = phi, range = range, sigma_O = sigma_O, sigma_E = 0,
      seed = 1, mesh = spde
    )
    spde <- make_mesh(s, c("x", "y"), n_knots = 50, type = "kmeans")
    m <- sdmTMB(data = s, formula = observed ~ 1, mesh = spde,
      family = student(link = "identity", df = 7),
      spatial = "off", spatiotemporal = "off"
    )
    expect_true(all(!is.na(summary(m$sd_report)[,"Std. Error"])))
    expect_length(residuals(m), nrow(s))
  })

  test_that("Lognormal fits", {
    skip_on_ci()
    skip_on_cran()
    range <- 1
    x <- stats::runif(500, -1, 1)
    y <- stats::runif(500, -1, 1)
    loc <- data.frame(x = x, y = y)
    spde <- make_mesh(loc, c("x", "y"), n_knots = 70, type = "kmeans")
    sigma_O <- 0.3
    sigma_E <- 0
    phi <- 0.2
    s <- sdmTMB_sim(x = x, y = y, time = 1L, mesh = spde, family = lognormal(),
      phi = phi, range = range, sigma_O = sigma_O, sigma_E = sigma_E, seed = 1
    )
    spde <- make_mesh(s, c("x", "y"), n_knots = 70, type = "kmeans")
    mlog <- sdmTMB(data = s, formula = observed ~ 1, mesh = spde,
      family = lognormal(link = "log"))
    expect_equal(exp(mlog$model$par[["ln_phi"]]), phi, tolerance = 0.1)
  })

  test_that("NB2 fits", {
    skip_on_ci()
    skip_on_cran()
    set.seed(1)
    x <- stats::runif(300, -1, 1)
    y <- stats::runif(300, -1, 1)
    loc <- data.frame(x = x, y = y)
    spde <- make_mesh(loc, c("x", "y"), n_knots = 80, type = "kmeans")
    s <- sdmTMB_sim(x = x, y = y, betas = 0.4, time = 1L, phi = 1.5, range = 0.8,
      sigma_O = 0.4, sigma_E = 0, seed = 1, mesh = spde, family = nbinom2())
    m <- sdmTMB(data = s, formula = observed ~ 1,
      mesh = spde, family = nbinom2(),
      control = sdmTMBcontrol(newton_loops = 1))
    expect_equal(round(tidy(m)[,"estimate"], 6), 0.274008)
  })

  test_that("Truncated NB2, truncated NB1, and regular NB1 fit", {
    skip_on_ci()
    skip_on_cran()
    set.seed(1)
    x <- stats::runif(300, -1, 1)
    y <- stats::runif(300, -1, 1)
    loc <- data.frame(x = x, y = y)
    spde <- make_mesh(loc, c("x", "y"), n_knots = 80, type = "kmeans")
    s <- sdmTMB_sim(x = x, y = y, betas = 0.4, time = 1L, phi = 1.5, range = 0.8,
      sigma_O = 0.4, sigma_E = 0, seed = 1, mesh = spde, family = nbinom2())

    m_sdmTMB <- sdmTMB(data = s, formula = observed ~ 1,
      mesh = spde, family = nbinom1(),
      spatial = "off")
    m_glmmTMB <- glmmTMB::glmmTMB(data = s, formula = observed ~ 1,
      family = glmmTMB::nbinom1())
    expect_equal(m_glmmTMB$fit$par[[1]], m_sdmTMB$model$par[[1]], tolerance = 0.00001)
    expect_equal(m_glmmTMB$fit$par[[2]], m_sdmTMB$model$par[[2]], tolerance = 0.00001)

    s_trunc <- subset(s, observed > 0)
    spde <- make_mesh(s_trunc, c("x", "y"), n_knots = 80, type = "kmeans")
    m_sdmTMB <- sdmTMB(data = s_trunc, formula = observed ~ 1,
      mesh = spde, family = truncated_nbinom2(), spatial = "off")
    m_glmmTMB <- glmmTMB::glmmTMB(data = s_trunc, formula = observed ~ 1,
      family = glmmTMB::truncated_nbinom2())
    expect_equal(m_glmmTMB$fit$par[[1]], m_sdmTMB$model$par[[1]], tolerance = 0.00001)
    expect_equal(m_glmmTMB$fit$par[[2]], m_sdmTMB$model$par[[2]], tolerance = 0.00001)

    m_sdmTMB <- sdmTMB(data = s_trunc, formula = observed ~ 1,
      mesh = spde, family = truncated_nbinom1(), spatial = "off")
    m_glmmTMB <- glmmTMB::glmmTMB(data = s_trunc, formula = observed ~ 1,
      family = glmmTMB::truncated_nbinom1())
    expect_equal(m_glmmTMB$fit$par[[1]], m_sdmTMB$model$par[[1]], tolerance = 0.00001)
    expect_equal(m_glmmTMB$fit$par[[2]], m_sdmTMB$model$par[[2]], tolerance = 0.00001)
  })

  test_that("Poisson fits", {
    skip_on_ci()
    skip_on_cran()
    d <- pcod
    spde <- make_mesh(pcod, c("X", "Y"), cutoff = 10)
    set.seed(3)
    d$density <- rpois(nrow(pcod), 3)
    m <- sdmTMB(data = d, formula = density ~ 1,
      mesh = spde, family = poisson(link = "log"),
      control = sdmTMBcontrol(newton_loops = 1)
    )
    expect_true(all(!is.na(summary(m$sd_report)[,"Std. Error"])))
    expect_length(residuals(m), nrow(pcod))
  })

  test_that("Binomial fits", {
    skip_on_ci()
    skip_on_cran()
    d <- pcod[pcod$year == 2017, ]
    d$density <- round(d$density)
    spde <- make_mesh(d, c("X", "Y"), cutoff = 10)
    d$present <- ifelse(d$density > 0, 1, 0)
    m <- sdmTMB(data = d, formula = present ~ 1,
      mesh = spde, family = binomial(link = "logit"),
      control = sdmTMBcontrol(newton_loops = 1))
    expect_true(all(!is.na(summary(m$sd_report)[,"Std. Error"])))
    expect_length(residuals(m), nrow(d))
  })

  test_that("Gamma fits", {
    skip_on_ci()
    skip_on_cran()
    d <- pcod[pcod$year == 2017 & pcod$density > 0, ]
    spde <- make_mesh(d, c("X", "Y"), cutoff = 10)
    m <- sdmTMB(data = d, formula = density ~ 1,
      mesh = spde, family = Gamma(link = "log"),
      spatial = "off",
      control = sdmTMBcontrol(newton_loops = 1))
    expect_true(all(!is.na(summary(m$sd_report)[,"Std. Error"])))
    expect_length(residuals(m), nrow(d))
    set.seed(123)
    d$test_gamma <- stats::rgamma(nrow(d), shape = 0.5, scale = 1 / 0.5)
    m <- sdmTMB(data = d, formula = test_gamma ~ 1,
      mesh = spde, family = Gamma(link = "inverse"), spatiotemporal = "off",
      control = sdmTMBcontrol(newton_loops = 1))
  })

  test_that("Beta fits", {
    skip_on_ci()
    skip_on_cran()
    set.seed(1)
    x <- stats::runif(400, -1, 1)
    y <- stats::runif(400, -1, 1)
    loc <- data.frame(x = x, y = y)
    spde <- make_mesh(loc, c("x", "y"), n_knots = 90, type = "kmeans")
    s <- sdmTMB_sim(x = x, y = y, sigma_E = 0, mesh = spde, sigma_O = 0.2,
      range = 0.8, family = Beta(), phi = 4)
    m <- sdmTMB(data = s, formula = observed ~ 1,
      mesh = spde, family = Beta(link = "logit"),
      control = sdmTMBcontrol(newton_loops = 1),
      spatial = "off")
    expect_true(all(!is.na(summary(m$sd_report)[,"Std. Error"])))
    expect_length(residuals(m), nrow(s))

    m_glmmTMB<- glmmTMB::glmmTMB(data = s, formula = observed ~ 1,
      family = glmmTMB::beta_family(link = "logit"))

    expect_equal(m$model$par[[2]], m_glmmTMB$fit$par[[2]], tolerance = 1e-4)
    expect_equal(m$model$par[[1]], m_glmmTMB$fit$par[[1]], tolerance = 1e-4)
  })
}

test_that("Censored Poisson fits", {
  skip_on_ci()
  skip_on_cran()
  set.seed(1)

  predictor_dat <- data.frame(X = runif(300), Y = runif(300))
  mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.2)
  sim_dat <- sdmTMB_simulate(
    formula = ~1,
    data = predictor_dat,
    mesh = mesh,
    family = poisson(),
    range = 0.5,
    sigma_O = 0.2,
    seed = 1,
    B = 2 # B0 = intercept
  )
  m_pois <- sdmTMB(
    data = sim_dat, formula = observed ~ 1,
    mesh = mesh, family = poisson(link = "log")
  )
  m_nocens_pois <- sdmTMB(
    data = sim_dat, formula = observed ~ 1,
    mesh = mesh, family = censored_poisson(link = "log"),
    experimental = list(upr = sim_dat$observed, lwr = sim_dat$observed)
  )
  expect_equal(m_nocens_pois$tmb_data$lwr, m_nocens_pois$tmb_data$upr)
  expect_equal(m_nocens_pois$tmb_data$lwr, as.numeric(m_nocens_pois$tmb_data$y_i))
  expect_equal(names(m_nocens_pois$tmb_data$family), "censored_poisson")
  expect_equal(m_pois$model, m_nocens_pois$model)

  # left-censored version
  L_1 <- 5 # zeros and ones cannot be observed directly - observed as <= L1
  y <- sim_dat$observed
  lwr <- ifelse(y <= L_1, 0, y)
  upr <- ifelse(y <= L_1, L_1, y)
  m_left_cens_pois <- sdmTMB(
    data = sim_dat, formula = observed ~ 1,
    mesh = mesh, family = censored_poisson(link = "log"),
    experimental = list(lwr = lwr, upr = upr),
    spatial = "off"
  )

  # right-censored version
  U_1 <- 8 # U_1 and above cannot be directly observed - instead we see >= U1
  y <- sim_dat$observed
  lwr <- ifelse(y >= U_1, U_1, y)
  upr <- ifelse(y >= U_1, NA, y)
  m_right_cens_pois <- sdmTMB(
    data = sim_dat, formula = observed ~ 1,
    mesh = mesh, family = censored_poisson(link = "log"),
    experimental = list(lwr = lwr, upr = upr),
    spatial = "off"
  )

  # interval-censored tough example
  # unique bounds per observation with upper limit 500 to test numerical underflow issues
  set.seed(123)
  U_2 <- sample(c(5:9), size = length(y), replace = TRUE)
  L_2 <- sample(c(1, 2, 3, 4), size = length(y), replace = TRUE)
  lwr <- ifelse(y >= U_2, U_2, ifelse(y <= L_2, 0, y))
  upr <- ifelse(y >= U_2, 500, ifelse(y <= L_2, L_2, y))
  m_interval_cens_pois <- sdmTMB(
    data = sim_dat, formula = observed ~ 1,
    mesh = mesh, family = censored_poisson(link = "log"),
    experimental = list(lwr = lwr, upr = upr),
    spatial = "off"
  )
  expect_true(all(!is.na(summary(m_interval_cens_pois$sd_report)[, "Std. Error"])))

  # reversed upr and lwr:
  expect_error(
    m <- sdmTMB(
      data = sim_dat, formula = observed ~ 1,
      mesh = mesh, family = censored_poisson(link = "log"),
      experimental = list(lwr = upr, upr = lwr)
    ), regexp = "lwr")

  # wrong length lwr and upr
  expect_error(
    m <- sdmTMB(
      data = sim_dat, formula = observed ~ 1,
      mesh = mesh, family = censored_poisson(link = "log"),
      experimental = list(lwr = c(1, 2), upr = c(4, 5, 6))
    ), regexp = "lwr")

  # missing lwr/upr
  expect_error(
    m <- sdmTMB(
      data = sim_dat, formula = observed ~ 1,
      mesh = mesh, family = censored_poisson(link = "log"),
    ), regexp = "lwr")

})
