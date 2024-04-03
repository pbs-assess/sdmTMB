test_that("simulate() and dharma_residuals() work", {
  skip_on_cran()

  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 15)
  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = mesh,
    family = tweedie()
  )
  s <- simulate(fit, nsim = 100)
  expect_equal(ncol(s), 100)
  expect_equal(nrow(s), nrow(pcod))
  expect_warning(dharma_residuals(s, fit), regexp = "mle")

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = mesh,
    family = delta_gamma()
  )
  s <- simulate(fit, nsim = 100, type = 'mle-mvn')
  expect_equal(ncol(s), 100)
  expect_equal(nrow(s), nrow(pcod))
  dharma_residuals(s, fit)

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = mesh,
    spatial = "off",
    family = delta_gamma(type = "poisson-link")
  )
  s <- simulate(fit, nsim = 100, type = "mle-mvn")
  expect_equal(ncol(s), 100)
  expect_equal(nrow(s), nrow(pcod))
  dharma_residuals(s, fit)
  dharma_residuals(s, fit, plot = FALSE)

  s <- simulate(fit, nsim = 100, params = "mle-mvn")
  expect_equal(ncol(s), 100)
  s <- simulate(fit, nsim = 100, params = "mle-mvn", re_form = NA)
  expect_equal(ncol(s), 100)

  s <- simulate(fit, nsim = 100, type = 'mle-mvn')
  dharma_residuals(s, fit, expected_distribution = "normal")
  x <- dharma_residuals(s, fit, return_DHARMa = TRUE)
  expect_s3_class(x, "DHARMa")
})
