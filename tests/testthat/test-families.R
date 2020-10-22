context("Families")

test_that("Families return a name to list with the correct names", {
  .names <- c("family", "link", "linkfun", "linkinv")
  expect_named(student(link = "identity"), .names)
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
spde <- make_mesh(loc, c("x", "y"), n_knots = 50, type = "kmeans")

test_that("Student family fits", {
  skip_on_travis()
  set.seed(1)
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
    family = student(link = "identity"))
  expect_true(all(!is.na(summary(m$sd_report)[,"Std. Error"])))
  expect_length(residuals(m), nrow(s))
})

test_that("Lognormal fits with a mean matching the Gamma roughly", {
  skip_on_travis()
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
  mgamma <- sdmTMB(data = s, formula = observed ~ 1, spde = spde,
    family = Gamma(link = "log"))
  expect_equal(mlog$model$par, mgamma$model$par, tolerance = 0.01)
})

test_that("NB2 fits", {
  d <- pcod[pcod$year == 2017, ]
  d$density <- round(d$density)
  spde <- make_mesh(d, c("X", "Y"), cutoff = 10)
  m <- sdmTMB(data = d, formula = density ~ 1,
    spde = spde, family = nbinom2(link = "log"))
  sdmTMBphi <- exp(m$model$par[["ln_phi"]])
  m2 <- glmmTMB::glmmTMB(density ~ 1, data = d,
    family = glmmTMB::nbinom2(link = "log"))
  glmmTMBphi <- exp(m2$fit$par[["betad"]])
  expect_equal(glmmTMBphi, sdmTMBphi, tol = 0.01)
})

test_that("Poisson fits", {
  d <- pcod
  spde <- make_mesh(pcod, c("X", "Y"), cutoff = 20)
  set.seed(3)
  d$density <- rpois(nrow(pcod), 3)
  m <- sdmTMB(data = d, formula = density ~ 1,
    spde = spde, family = poisson(link = "log"))
  expect_true(all(!is.na(summary(m$sd_report)[,"Std. Error"])))
  expect_length(residuals(m), nrow(pcod))
})

test_that("Binomial fits", {
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
  d <- pcod[pcod$year == 2017 & pcod$density > 0, ]
  # d$density <- d$density / 100
  spde <- make_mesh(d, c("X", "Y"), cutoff = 10)
  m <- sdmTMB(data = d, formula = density ~ 1,
    spde = spde, family = Gamma(link = "log"))
  expect_true(all(!is.na(summary(m$sd_report)[,"Std. Error"])))
  expect_length(residuals(m), nrow(d))
})

set.seed(1)
x <- stats::runif(400, -1, 1)
y <- stats::runif(400, -1, 1)
loc <- data.frame(x = x, y = y)
spde <- make_mesh(loc, c("x", "y"), n_knots = 90, type = "kmeans")

test_that("Beta fits", {
  skip_on_travis()
  s <- sdmTMB_sim(x = x, y = y, sigma_E = 0, mesh = spde, sigma_O = 0.2, range = 0.8, family = Beta(), phi = 4)
  m <- sdmTMB(data = s, formula = observed ~ 1,
    spde = spde, family = Beta(link = "logit"))
  expect_true(all(!is.na(summary(m$sd_report)[,"Std. Error"])))

  # m2 <- glmmTMB::glmmTMB(observed ~ 1, data = s,
  #   family = glmmTMB::beta_family(link = "logit"))
  # glmmTMBphi <- exp(m2$fit$par[["betad"]])
  #
  # expect_equal(m$model$par[["ln_phi"]], log(glmmTMBphi), tol = 0.1)
  expect_length(residuals(m), nrow(s))
})
