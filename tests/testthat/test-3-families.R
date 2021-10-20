test_that("Families return a name to list with the correct names", {
  .names <- c("family", "link", "linkfun", "linkinv")
  expect_named(student(link = "identity"), c(.names, "df"))
  expect_named(lognormal(link = "log"), .names)
  expect_named(tweedie(link = "log"), .names)
  expect_named(nbinom2(link = "log"), .names)
})

test_that("The supplementary families work with appropriate links", {
  expect_identical(class(tweedie(link = "log")), "list")
  expect_identical(class(tweedie(link = log)), "list")
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
    m <- sdmTMB(data = s, formula = observed ~ 1, spde = spde,
      family = student(link = "identity", df = 7),
      control = sdmTMBcontrol(map_rf = TRUE)
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
    mlog <- sdmTMB(data = s, formula = observed ~ 1, spde = spde,
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
      spde = spde, family = nbinom2())
    expect_equal(round(tidy(m)[,"estimate"], 6), 0.274008)
  })

  test_that("Truncated NB2 fits", {
    skip_on_ci()
    skip_on_cran()
    set.seed(1)
    x <- stats::runif(300, -1, 1)
    y <- stats::runif(300, -1, 1)
    loc <- data.frame(x = x, y = y)
    spde <- make_mesh(loc, c("x", "y"), n_knots = 80, type = "kmeans")
    s <- sdmTMB_sim(x = x, y = y, betas = 0.4, time = 1L, phi = 1.5, range = 0.8,
      sigma_O = 0.4, sigma_E = 0, seed = 1, mesh = spde, family = nbinom2())
    s_trunc <- subset(s, observed > 0)
    spde <- make_mesh(s_trunc, c("x", "y"), n_knots = 80, type = "kmeans")
    m_sdmTMB <- sdmTMB(data = s_trunc, formula = observed ~ 1,
      spde = spde, family = truncated_nbinom2(), control = sdmTMBcontrol(map_rf = TRUE))
    m_glmmTMB <- glmmTMB::glmmTMB(data = s_trunc, formula = observed ~ 1,
      family = glmmTMB::truncated_nbinom2())
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
      spde = spde, family = poisson(link = "log"))
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
      spde = spde, family = binomial(link = "logit"))
    expect_true(all(!is.na(summary(m$sd_report)[,"Std. Error"])))
    expect_length(residuals(m), nrow(d))
  })

  test_that("Gamma fits", {
    skip_on_ci()
    skip_on_cran()
    d <- pcod[pcod$year == 2017 & pcod$density > 0, ]
    spde <- make_mesh(d, c("X", "Y"), cutoff = 10)
    m <- sdmTMB(data = d, formula = density ~ 1,
      spde = spde, family = Gamma(link = "log"))
    expect_true(all(!is.na(summary(m$sd_report)[,"Std. Error"])))
    expect_length(residuals(m), nrow(d))
    set.seed(123)
    d$test_gamma <- stats::rgamma(nrow(d), shape = 0.5, scale = 1 / 0.5)
    m <- sdmTMB(data = d, formula = test_gamma ~ 1,
      spde = spde, family = Gamma(link = "inverse"), spatiotemporal = "off")
    expect_true(all(!is.na(summary(m$sd_report)[,"Std. Error"])))
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
      spde = spde, family = Beta(link = "logit"))
    expect_true(all(!is.na(summary(m$sd_report)[,"Std. Error"])))
    expect_length(residuals(m), nrow(s))
  })

}
