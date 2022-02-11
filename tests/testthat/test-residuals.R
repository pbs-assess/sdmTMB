test_that("Simulated residuals work", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  fit <- sdmTMB(density ~ as.factor(year) + poly(depth, 2),
    data = pcod_2011, time = "year", mesh = pcod_mesh_2011,
    family = tweedie(link = "log"), spatial = "off",
    spatiotemporal = "off")
  s <- simulate(fit, nsim = 50)
  dharma_residuals(s, fit)
  r <- dharma_residuals(s, fit, plot = FALSE)
  expect_equal(class(r), "data.frame")
  expect_error(dharma_residuals(c(1, 2, 3), fit))
  expect_error(dharma_residuals(matrix(c(1, 2, 3)), fit))
  expect_error(dharma_residuals(s, fit, plot = "test"))
  expect_error(dharma_residuals(s, 99))
})
