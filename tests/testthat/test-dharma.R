test_that("simulate() and dharma_residuals() work", {
  skip_on_ci()
  skip_on_cran()

  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 15)
  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = mesh,
    family = tweedie()
  )
  s <- simulate(fit, nsim = 100)
  expect_equal(ncol(s), 100)
  expect_equal(nrow(s), nrow(pcod))
  # dharma_residuals(s, fit)

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = mesh,
    family = delta_gamma()
  )
  s <- simulate(fit, nsim = 100)
  expect_equal(ncol(s), 100)
  expect_equal(nrow(s), nrow(pcod))
  # dharma_residuals(s, fit)

  # fit <- sdmTMB(density ~ 1,
  #   data = pcod, mesh = mesh,
  #   spatial = "off",
  #   family = delta_poisson_link_gamma()
  # )
  # s <- simulate(fit, nsim = 100)
  # expect_equal(ncol(s), 100)
  # expect_equal(nrow(s), nrow(pcod))
  # dharma_residuals(s, fit)

  s <- simulate(fit, nsim = 100, params = "mvn")
  expect_equal(ncol(s), 100)
  s <- simulate(fit, nsim = 100, params = "mvn", re_form = NA)
  expect_equal(ncol(s), 100)
})
