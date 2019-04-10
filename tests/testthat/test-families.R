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

test_that("Student and lognormal families fit", {
  set.seed(1)
  initial_betas <- 0.5
  kappa <- 2
  sigma_O <- 0.3
  phi <- 0.01
  x <- stats::runif(100, -1, 1)
  y <- stats::runif(100, -1, 1)
  s <- sim(
    initial_betas = initial_betas, time = 1L,
    phi = phi, kappa = kappa, sigma_O = sigma_O, sigma_E = sigma_E,
    seed = 1, plot = TRUE
  )
  spde <- make_spde(s$x, s$y, n_knots = 50)
  m <- sdmTMB(data = s, formula = observed ~ 1, spde = spde,
    family = student(link = "identity"))
  expect_true(all(!is.na(summary(m$sd_report)[,"Std. Error"])))
})

test_that("NB2 fits", {
  d <- pcod[pcod$year == 2017, ]
  d$density <- round(d$density)
  spde <- make_spde(d$X, d$Y, n_knots = 30)
  m <- sdmTMB(data = d, formula = density ~ 1,
    spde = spde, family = nbinom2(link = "log"))
  sdmTMBphi <- exp(m$model$par[["ln_phi"]])
  m2 <- glmmTMB::glmmTMB(density ~ 1, data = d,
    family = glmmTMB::nbinom2(link = "log"))
  glmmTMBphi <- exp(m2$fit$par[["betad"]])
  expect_equal(glmmTMBphi, sdmTMBphi, tol = 0.01)
})

test_that("Poisson fits", {
  d <- pcod[pcod$year == 2017, ]
  d$density <- round(d$density)
  spde <- make_spde(d$X, d$Y, n_knots = 30)
  m <- sdmTMB(data = d, formula = density ~ 1,
    spde = spde, family = poisson(link = "log"))
  expect_true(all(!is.na(summary(m$sd_report)[,"Std. Error"])))
})


test_that("Binomial fits", {
  d <- pcod[pcod$year == 2017, ]
  d$density <- round(d$density)
  spde <- make_spde(d$X, d$Y, n_knots = 30)
  d$present <- ifelse(d$density > 0, 1, 0)
  m <- sdmTMB(data = d, formula = present ~ 1,
    spde = spde, family = binomial(link = "logit"))
  expect_true(all(!is.na(summary(m$sd_report)[,"Std. Error"])))
})

test_that("Gamma fits", {
  d <- pcod[pcod$year == 2017 & pcod$density > 0, ]
  # d$density <- d$density / 100
  spde <- make_spde(d$X, d$Y, n_knots = 30)
  m <- sdmTMB(data = d, formula = density ~ 1,
    spde = spde, family = Gamma(link = "log"))
  expect_true(all(!is.na(summary(m$sd_report)[,"Std. Error"])))
})
