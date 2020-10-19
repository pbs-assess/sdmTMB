context("Families")

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

test_that("Student and family fits", {
  set.seed(1)
  initial_betas <- 0.5
  kappa <- 2
  sigma_O <- 0.3
  phi <- 0.01
  x <- stats::runif(100, -1, 1)
  y <- stats::runif(100, -1, 1)
  s <- sim(x = x, y = y,
    initial_betas = initial_betas, time = 1L,
    phi = phi, kappa = kappa, sigma_O = sigma_O, sigma_E = 0.0001,
    seed = 1
  )
  spde <- make_mesh(s, c("x", "y"), n_knots = 50, type = "kmeans")
  m <- sdmTMB(data = s, formula = observed ~ 1, spde = spde,
    family = student(link = "identity", df = 7))
  expect_true(all(!is.na(summary(m$sd_report)[,"Std. Error"])))
  expect_length(residuals(m), nrow(s))
})

test_that("Lognormal fits with a mean matching the Gamma roughly", {
  kappa <- .5
  x <- stats::runif(500, -1, 1)
  y <- stats::runif(500, -1, 1)
  sigma_O <- 0.3
  sigma_E <- 0.0001
  phi <- 0.2
  s <- sim(x = x, y = y, time = 1L,
    phi = phi, kappa = kappa, sigma_O = sigma_O, sigma_E = sigma_E, seed = 1
  )
  spde <- make_mesh(s, c("x", "y"), n_knots = 70, type = "kmeans")
  s$observed <- stats::rlnorm(nrow(s), mean = log(exp(s$real)), sd = 0.2)
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

test_that("Beta fits", {
  s <- sim(sigma_O = 0.02)
  s$observed <- stats::plogis(s$observed * 7)
  spde <- make_mesh(s, c("x", "y"), cutoff = 0.02)
  m <- sdmTMB(data = s, formula = observed ~ 1,
    spde = spde, family = Beta(link = "logit"))
  expect_true(all(!is.na(summary(m$sd_report)[,"Std. Error"])))

  m2 <- glmmTMB::glmmTMB(observed ~ 1, data = s,
    family = glmmTMB::beta_family(link = "logit"))
  glmmTMBphi <- exp(m2$fit$par[["betad"]])

  expect_equal(m$model$par[["ln_phi"]], log(glmmTMBphi), tol = 0.1)
  expect_length(residuals(m), nrow(s))
})
