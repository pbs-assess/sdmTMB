test_that("cAIC and EDF work", {
  skip_on_cran()
  skip_on_ci()
  mesh <- make_mesh(dogfish, c("X", "Y"), cutoff = 15)
  suppressMessages(
    fit <- sdmTMB(catch_weight ~ s(log(depth)),
      time_varying = ~1,
      time_varying_type = "ar1",
      time = "year",
      spatiotemporal = "off",
      mesh = mesh,
      family = tweedie(),
      data = dogfish,
      offset = log(dogfish$area_swept)
    )
  )
  expect_equal(AIC(fit), 12192.9613, tolerance = 1e-4)
  expect_equal(cAIC(fit), 12071.4289, tolerance = 1e-4)
  edf <- cAIC(fit, what = "EDF")
  expect_equal(sum(edf), 54.3870623, tolerance = 1e-4)
})
