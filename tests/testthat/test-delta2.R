test_that("share_range mapping works with delta models", {
  skip_on_cran()
  skip_if_not_installed("INLA")

  pcod_spde <- make_mesh(pcod, c("X", "Y"), cutoff = 15)

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde, spatial = "off",
    spatiotemporal = list("off", "rw"), time = "year",
    do_fit = FALSE, family = delta_gamma(),
    share_range = TRUE,
  )
  expect_identical(fit$tmb_map$ln_kappa, as.factor(c(NA, NA, 1, 1)))

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde, spatial = "off",
    spatiotemporal = list("iid", "iid"), time = "year",
    do_fit = FALSE, family = delta_gamma(),
    share_range = list(TRUE, TRUE),
  )
  expect_identical(fit$tmb_map$ln_kappa, as.factor(c(1, 1, 2, 2)))

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde, spatial = "off",
    spatiotemporal = list("iid", "iid"), time = "year",
    do_fit = FALSE, family = delta_gamma(),
    share_range = list(FALSE, TRUE),
  )
  expect_identical(fit$tmb_map$ln_kappa, as.factor(c(1, 2, 3, 3)))

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde, spatial = "off",
    spatiotemporal = list("iid", "iid"), time = "year",
    do_fit = FALSE, family = delta_gamma(),
    share_range = list(FALSE, FALSE),
  )
  expect_identical(fit$tmb_map$ln_kappa, as.factor(c(1, 2, 3, 4)))

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde, spatial = "off",
    spatiotemporal = list("iid", "iid"), time = "year",
    do_fit = FALSE, family = delta_gamma(),
    share_range = FALSE,
  )
  expect_identical(fit$tmb_map$ln_kappa, as.factor(c(1, 2, 3, 4)))

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde, spatial = "off",
    spatiotemporal = list("iid", "iid"), time = "year",
    do_fit = FALSE, family = delta_gamma(),
    share_range = TRUE,
  )
  expect_identical(fit$tmb_map$ln_kappa, as.factor(c(1, 1, 2, 2)))

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde, spatial = "off",
    spatiotemporal = list("off", "off"), time = "year",
    do_fit = FALSE, family = delta_gamma(),
    share_range = TRUE,
  )
  expect_identical(fit$tmb_map$ln_kappa, as.factor(c(NA, NA, NA, NA)))

  # spatial models:
  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde,
    do_fit = TRUE, family = tweedie()
  )
  expect_identical(fit$tmb_map$ln_kappa, as.factor(c(1, 1)))

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde,
    do_fit = TRUE, family = delta_gamma()
  )
  expect_identical(fit$tmb_map$ln_kappa, as.factor(c(1, 1, 2, 2)))

})

test_that("spatiotemporal field mapping/specification works with delta models", {
  skip_on_cran()
  skip_if_not_installed("INLA")

  pcod_spde <- make_mesh(pcod, c("X", "Y"), cutoff = 20)

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde, spatial = "off",
    time = "year", family = delta_gamma()
  )
  s1 <- as.list(fit$sd_report, "Estimate")
  s1$ln_tau_E

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde, spatial = "off",
    time = "year", family = delta_gamma(),
    spatiotemporal = list("iid", "iid")
  )
  s2 <- as.list(fit$sd_report, "Estimate")
  s2$ln_tau_E
  expect_equal(s1, s2)

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde, spatial = "off",
    time = "year", family = delta_gamma(),
    spatiotemporal = list("iid", "off")
  )
  s <- as.list(fit$sd_report, "Estimate")
  expect_gt(s$ln_tau_E[1], 0)
  expect_equal(s$ln_tau_E[2], 0)
  expect_output(print(fit), regexp = "Spatiotemporal model")
  expect_output(print(fit), regexp = "Spatiotemporal SD")

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde, spatial = "off",
    time = "year", family = delta_gamma(),
    spatiotemporal = list("off", "iid")
  )
  s <- as.list(fit$sd_report, "Estimate")
  expect_gt(s$ln_tau_E[2], 0)
  expect_equal(s$ln_tau_E[1], 0)
  print(fit)

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde, spatial = "off",
    time = "year", family = delta_gamma(),
    spatiotemporal = list("off", "off")
  )
  s <- as.list(fit$sd_report, "Estimate")
  expect_equal(s$ln_tau_E[1], 0)
  expect_equal(s$ln_tau_E[2], 0)
  print(fit)
  expect_output(print(fit), regexp = "Model fit")

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde, spatial = "off",
    time = "year", family = delta_gamma(),
    spatiotemporal = list("ar1", "off")
  )
  s <- as.list(fit$sd_report, "Estimate")
  expect_gt(s$ar1_phi[1], 0)
  expect_equal(s$ar1_phi[2], 0)
  expect_identical(fit$tmb_map$ar1_phi, as.factor(c(1, NA)))
  tidy(fit, "ran_pars", model = 1)
  tidy(fit, "ran_pars", model = 2)
  print(fit)
  expect_output(print(fit), regexp = "rho")

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde, spatial = "off",
    time = "year", family = delta_gamma(),
    spatiotemporal = list("off", "ar1")
  )
  s <- as.list(fit$sd_report, "Estimate")
  expect_gt(s$ar1_phi[2], 0)
  expect_equal(s$ar1_phi[1], 0)
  expect_identical(fit$tmb_map$ar1_phi, as.factor(c(NA, 1)))
  print(fit)
  expect_output(print(fit), regexp = "rho")

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde, spatial = "off",
    time = "year", family = delta_gamma(),
    spatiotemporal = list("rw", "off")
  )
  s <- as.list(fit$sd_report, "Estimate")
  expect_gt(s$ln_tau_E[1], 0)
  expect_equal(s$ln_tau_E[2], 0)

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde, spatial = "off",
    time = "year", family = delta_gamma(),
    spatiotemporal = list("ar1", "ar1")
  )
  s1 <- as.list(fit$sd_report, "Estimate")
  print(fit)
  expect_output(print(fit), regexp = "rho")

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde, spatial = "off",
    time = "year", family = delta_gamma(),
    spatiotemporal = "ar1"
  )
  s2 <- as.list(fit$sd_report, "Estimate")
  expect_equal(s1, s2)
})

test_that("delta models work with different main effects", {
  skip_on_cran()
  skip_if_not_installed("INLA")

  mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)

  fit <- sdmTMB(
    formula = list(
      density ~ depth_scaled,
      density ~ depth_scaled + depth_scaled2
    ),
    data = pcod_2011,
    mesh = mesh,
    spatial = "off",
    family = delta_gamma()
  )
  s <- as.list(fit$sd_report, "Estimate")
  se <- as.list(fit$sd_report, "Std. Error")
  expect_true("b_j2" %in% names(s))
  expect_true(all(!is.na(se$b_j2)))
  expect_true(all(!is.na(se$b_j)))
  tidy(fit)
  fit
  p <- predict(fit)

  # still works?
  fit <- sdmTMB(
    formula = density ~ depth_scaled + depth_scaled2,
    data = pcod_2011,
    mesh = mesh,
    spatial = "off",
    family = delta_gamma()
  )
  fit
  p <- predict(fit)

  # one smoother works
  fit <- sdmTMB(
    formula = density ~ s(depth_scaled),
    data = pcod_2011,
    mesh = mesh,
    spatial = "off",
    family = delta_gamma()
  )

  # should be same:
  fit2 <- sdmTMB(
    formula = list(density ~ s(depth_scaled), density ~ s(depth_scaled)),
    data = pcod_2011,
    mesh = mesh,
    spatial = "off",
    family = delta_gamma()
  )
  expect_equal(fit$sd_report, fit2$sd_report)

  # diff smoothers throw an error for now
  expect_error(
    {
      fit <- sdmTMB(
        formula = list(density ~ s(depth_scaled), density ~ s(year, k = 3)),
        data = pcod_2011,
        mesh = mesh,
        spatial = "off",
        family = delta_gamma()
      )
    },
    regexp = "smooth"
  )

  # diff random intercepts throw an error for now
  expect_error(
    {
      pcod_2011$fyear <- as.factor(pcod_2011$year)
      fit <- sdmTMB(
        formula = list(density ~ 1 + (1 | fyear), density ~ 1),
        data = pcod_2011,
        mesh = mesh,
        spatial = "off",
        family = delta_gamma()
      )
    },
    regexp = "random intercepts"
  )

  # OK:
  fit <- sdmTMB(
    formula = list(density ~ 1 + (1 | fyear), density ~ 1 + (1 | fyear)),
    data = pcod_2011,
    mesh = mesh,
    spatial = "off",
    family = delta_gamma()
  )
})
