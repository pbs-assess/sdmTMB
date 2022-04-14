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
