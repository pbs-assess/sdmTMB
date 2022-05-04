test_that("sanity() runs", {
  skip_on_cran()
  skip_if_not_installed("INLA")

  fit <- sdmTMB(density ~ 1, time = "year",
    data = pcod_2011, mesh = pcod_mesh_2011,
    family = tweedie()
  )
  s <- sanity(fit)
  expect_true(all(unlist(s)))

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
  expect_false(s$se_magnitude_ok)
  expect_false(s$all_ok)
})
