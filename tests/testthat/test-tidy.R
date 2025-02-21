test_that("tidy works", {
  fit <- sdmTMB(
    density ~ depth_scaled,
    data = pcod_2011, mesh = pcod_mesh_2011,
    family = tweedie(link = "log"), time = "year",
  )
  x <- tidy(fit, "ran_pars", conf.int = TRUE)
  x
  expect_true(sum(is.na(x$std.error)) == 0L)

  x <- tidy(fit, conf.int = TRUE)
  expect_true(sum(is.na(x$std.error)) == 0L)

  d <- pcod_2011
  d$year_scaled <- as.numeric(scale(d$year))
  fit <- sdmTMB(
    density ~ year_scaled,
    data = d, mesh = pcod_mesh_2011,
    family = tweedie(link = "log"), time = "year",
    spatiotemporal = "off", spatial_varying = ~ 0 + year_scaled
  )

  x <- tidy(fit, "ran_pars", conf.int = TRUE)
  x
  expect_true(sum(is.na(x$std.error)) == 0L)

  fit <- sdmTMB(
    density ~ depth_scaled,
    data = pcod_2011, mesh = pcod_mesh_2011,
    family = delta_gamma()
  )
  x <- tidy(fit, "ran_pars", conf.int = TRUE)
  x
  expect_true(sum(is.na(x$std.error)) == 0L)

  x <- tidy(fit, "ran_pars", conf.int = TRUE, model = 2)
  x
  expect_true(sum(is.na(x$std.error)) == 0L)

  fit <- sdmTMB(
    density ~ depth_scaled,
    data = pcod_2011, mesh = pcod_mesh_2011, share_range = FALSE,
    family = tweedie(link = "log"), time = "year",
  )
  x <- tidy(fit, "ran_pars", conf.int = TRUE)
  x
  expect_true(sum(is.na(x$std.error)) == 0L)

  # test that parsing of time varying random values works
  mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
  fit <- sdmTMB(
    density ~ 0 + as.factor(year),
    time_varying = ~ 0 + depth_scaled + depth_scaled2, #<
    data = pcod_2011, time = "year", mesh = mesh,
    family = tweedie()
  )
  pars <- tidy(fit, "ran_vals")
  expect_equal(pars$estimate, c(-0.87, -0.81, -0.75, -1.11,
                                -1.92, -0.92, -1.59, -2.20), tolerance = 0.01)
  expect_equal(pars$term, c("depth_scaled:2011","depth_scaled2:2013","depth_scaled:2015","depth_scaled2:2017",
                            "depth_scaled:2011","depth_scaled2:2013","depth_scaled:2015","depth_scaled2:2017"))


})
