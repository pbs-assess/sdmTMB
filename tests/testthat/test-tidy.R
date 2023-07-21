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
})
