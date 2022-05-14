test_that("sanity() runs", {
  skip_on_cran()
  skip_if_not_installed("INLA")

  fit <- sdmTMB(density ~ 1, time = "year",
    data = pcod_2011, mesh = pcod_mesh_2011,
    family = tweedie()
  )
  s <- sanity(fit)
  expect_true(all(unlist(s)))

  fit$gradients[1] <- 0.9
  expect_message(sanity(fit), regexp = "newton")

  fit$gradients[1] <- 0.00001
  fit$bad_eig <- TRUE
  expect_message(sanity(fit), regexp = "model may not have converged")

  fit$bad_eig <- FALSE
  fit$model$convergence <- 1L
  expect_message(sanity(fit), regexp = "did not converge")
  fit$model$convergence <- 0L

  fit$pos_def_hessian <- FALSE
  expect_message(sanity(fit), regexp = "Non-positive-definite")

  fit <- sdmTMB(density ~ 1, time = "year",
    data = pcod_2011, mesh = pcod_mesh_2011,
    family = delta_gamma()
  )
  s <- sanity(fit)
  expect_true(all(unlist(s)))

  fit <- sdmTMB(density ~ s(depth), time = "year",
    data = pcod_2011, mesh = pcod_mesh_2011,
    family = delta_gamma()
  )
  expect_message(s <- sanity(fit), regexp = "standard error")
  expect_message(s <- sanity(fit), regexp = "may be large")
  expect_message(s <- sanity(fit, se_ratio = 9), regexp = "9x")
  expect_false(s$se_magnitude_ok)
  expect_false(s$all_ok)
})
