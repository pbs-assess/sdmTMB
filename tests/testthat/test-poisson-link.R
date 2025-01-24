test_that("Poisson-link prediction issues from GitHub issues #389,", {
  skip_on_cran()
  dogfish$log_depth <- log(dogfish$depth)
  mesh <- make_mesh(dogfish, c("X", "Y"), cutoff = 30)
  fit_dpg <- sdmTMB(catch_weight ~ 0 + as.factor(year) + s(log_depth),
    family = delta_gamma(type = "poisson-link"),
    spatial = "on",
    mesh = mesh,
    data = dogfish,
    offset = log(dogfish$area_swept)
  )
  p <- predict(fit_dpg, re_form = NA)
  p$est_est1est2 <- p$est1 + p$est2
  expect_equal(p$est, p$est_est1est2)
  expect_equal(p$est_est1est2[1:5],
    c(7.01028301430728, 2.34143881755487, 6.96979232578834, 6.99973559970208,
      7.03187981132451), tolerance = 1e-3)

  pp <- predict(fit_dpg, type = "response")
  expect_equal(pp$est[1:5],
    c(693.502617933497, 107.803298324919, 8414.54507288536, 5770.52404422525,
    6545.06096568627), tolerance = 1e-3)
})
