test_that("project() works with delta models", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  library(ggplot2)
  mesh <- make_mesh(dogfish, c("X", "Y"), cutoff = 35)
  historical_years <- 2004:2022
  to_project <- 1
  future_years <- seq(max(historical_years) + 1, max(historical_years) + to_project)
  all_years <- c(historical_years, future_years)
  proj_grid <- replicate_df(wcvi_grid, "year", all_years)

  # we could fit our model like this, but for long projections, this becomes slow:
  fit <- sdmTMB(
    catch_weight ~ 1,
    time = "year",
    offset = log(dogfish$area_swept),
    extra_time = all_years, #< note that all years here
    spatial = "off", # speed
    spatiotemporal = "ar1",
    data = dogfish,
    mesh = mesh,
    family = delta_gamma()
  )
  p <- predict(fit, newdata = proj_grid)

  # instead, we could fit our model like this and then take simulation draws
  # from the projection time period:
  fit2 <- sdmTMB(
    catch_weight ~ 1,
    time = "year",
    offset = log(dogfish$area_swept),
    extra_time = historical_years, #< does *not* include projection years
    spatial = "off", # speed
    spatiotemporal = "ar1",
    data = dogfish,
    mesh = mesh,
    family = delta_gamma()
  )
  set.seed(1)
  out <- project(fit2, newdata = proj_grid, nsim = 100, uncertainty = "none")
  none <- out
  expect_identical(names(out), c("est1", "est2", "epsilon_st1", "epsilon_st2"))

  proj_grid$est_mean1 <- apply(out$est1, 1, mean)
  proj_grid$est_mean2 <- apply(out$est2, 1, mean)
  proj_grid$eps_mean1 <- apply(out$epsilon_st1, 1, mean)
  proj_grid$eps_mean2 <- apply(out$epsilon_st2, 1, mean)

  # visualize:
  ggplot(subset(proj_grid, year > 2021), aes(X, Y, fill = est_mean1)) +
    geom_raster() +
    facet_wrap(~year) +
    coord_fixed() +
    scale_fill_viridis_c(limits = c(-4, 4)) +
    ggtitle("Projection simulation (mean)")

  # or with predict() method:
  ggplot(subset(p, year > 2021), aes(X, Y, fill = est1)) +
    geom_raster() +
    facet_wrap(~year) +
    coord_fixed() +
    scale_fill_viridis_c() +
    labs(fill = "est_mean") +
    scale_fill_viridis_c(limits = c(-4, 4)) +
    ggtitle("Projection simulation (mean)")

  i <- p$year == 2023
  plot(p$est1[i], proj_grid$est_mean1[i])
  plot(p$est2[i], proj_grid$est_mean2[i])

  plot(p$epsilon_st1[i], proj_grid$eps_mean1[i])
  plot(p$epsilon_st2[i], proj_grid$eps_mean2[i])

  OK_COR <- 0.98
  expect_gt(cor(p$est1[i], proj_grid$est_mean1[i]), OK_COR)
  expect_gt(cor(p$est2[i], proj_grid$est_mean2[i]), OK_COR)
  expect_gt(cor(p$epsilon_st1[i], proj_grid$eps_mean1[i]), OK_COR)
  expect_gt(cor(p$epsilon_st2[i], proj_grid$eps_mean2[i]), OK_COR)

  # test return_tmb_report:

  # if instead we wanted to grab, say, the spatiotemporal random field values,
  # we can return the report and work with the raw output ourselves:
  set.seed(1)
  out <- project(
    fit2,
    newdata = proj_grid,
    nsim = 100,
    uncertainty = "none",
    return_tmb_report = TRUE #< difference from above example
  )

  eps <- lapply(out, \(x) x[["epsilon_st_A_vec"]][, 1])
  eps <- do.call(cbind, eps)
  eps_mean <- apply(eps, 1, mean)
  plot(eps_mean, proj_grid$eps_mean1)
  expect_gt(cor(eps_mean, proj_grid$eps_mean1), OK_COR)

  eps <- lapply(out, \(x) x[["epsilon_st_A_vec"]][, 2])
  eps <- do.call(cbind, eps)
  eps_mean <- apply(eps, 1, mean)
  plot(eps_mean, proj_grid$eps_mean2)
  expect_gt(cor(eps_mean, proj_grid$eps_mean2), OK_COR)

  # test the types of uncertainty
  set.seed(1)
  both <- project(fit2, newdata = proj_grid, nsim = 50, uncertainty = "both")
  set.seed(1)
  suppressWarnings({
    random <- project(fit2, newdata = proj_grid, nsim = 50, uncertainty = "random")
  })

  sd_both <- mean(apply(both$est1, 1, sd)[i])
  sd_none <- mean(apply(none$est1, 1, sd)[i])
  sd_random <- mean(apply(random$est1, 1, sd)[i], na.rm = TRUE)
  sd_both
  sd_none
  sd_random
  # expect_gt(sd_both, sd_random) # !!?
  expect_gt(sd_both, sd_none)
  expect_gt(sd_random, sd_none)
})

test_that("project() works with non-delta models", {
  skip_on_cran()

  mesh <- make_mesh(dogfish, c("X", "Y"), cutoff = 35)
  historical_years <- 2004:2022
  to_project <- 1
  future_years <- seq(max(historical_years) + 1, max(historical_years) + to_project)
  all_years <- c(historical_years, future_years)
  proj_grid <- replicate_df(wcvi_grid, "year", all_years)

  fit <- sdmTMB(
    catch_weight ~ 1,
    time = "year",
    offset = log(dogfish$area_swept),
    extra_time = all_years, #< note that all years here
    spatial = "off", # speed
    spatiotemporal = "ar1",
    data = dogfish,
    mesh = mesh,
    family = tweedie()
  )
  p <- predict(fit, newdata = proj_grid)
  fit2 <- sdmTMB(
    catch_weight ~ 1,
    time = "year",
    offset = log(dogfish$area_swept),
    extra_time = historical_years, #< does *not* include projection years
    spatial = "off", # speed
    spatiotemporal = "ar1",
    data = dogfish,
    mesh = mesh,
    family = tweedie()
  )
  set.seed(1)
  out <- project(fit2, newdata = proj_grid, nsim = 100, uncertainty = "none")
  expect_identical(names(out), c("est", "epsilon_st"))

  i <- p$year == 2023
  proj_grid$est_mean <- apply(out$est, 1, mean)
  proj_grid$eps_mean <- apply(out$epsilon_st, 1, mean)
  plot(p$est[i], proj_grid$est_mean[i])
  plot(p$epsilon_st[i], proj_grid$eps_mean[i])
  expect_gt(cor(p$est[i], proj_grid$est_mean[i]), 0.98)
  expect_gt(cor(p$epsilon_st[i], proj_grid$est_mean[i]), 0.98)
})

test_that("project() works with time-varying effects", {
  skip_on_cran()
  mesh <- make_mesh(dogfish, c("X", "Y"), cutoff = 35)
  historical_years <- 2004:2022
  to_project <- 1
  future_years <- seq(max(historical_years) + 1, max(historical_years) + to_project)
  all_years <- c(historical_years, future_years)
  proj_grid <- replicate_df(wcvi_grid, "year", all_years)

  fit <- sdmTMB(
    catch_weight ~ 1,
    time = "year",
    offset = log(dogfish$area_swept),
    time_varying = ~1,
    time_varying_type = "ar1",
    extra_time = all_years, #< note that all years here
    spatial = "off",
    spatiotemporal = "off",
    data = dogfish,
    # mesh = mesh,
    family = tweedie()
  )
  p <- predict(fit, newdata = proj_grid)
  fit2 <- sdmTMB(
    catch_weight ~ 1,
    time = "year",
    offset = log(dogfish$area_swept),
    time_varying = ~1,
    time_varying_type = "ar1",
    extra_time = historical_years, #< does *not* include projection years
    spatial = "off",
    spatiotemporal = "off",
    data = dogfish,
    # mesh = mesh,
    family = tweedie()
  )
  set.seed(1)
  out <- project(fit2, newdata = proj_grid, nsim = 100, uncertainty = "none")
  expect_identical(names(out), c("est"))
  i <- p$year == 2023
  hist(out$est[i, ])
  abline(v = mean(p$est[i]))
  expect_equal(mean(p$est[i]), 5.983172, tolerance = 1e-3)
})


test_that("project() works/fails as expected in some less obvious situations", {
  skip_on_cran()

  mesh <- make_mesh(dogfish, c("X", "Y"), cutoff = 35)
  historical_years <- 2004:2022
  to_project <- 1
  future_years <- seq(max(historical_years) + 1, max(historical_years) + to_project)
  all_years <- c(historical_years, future_years)
  proj_grid <- replicate_df(wcvi_grid, "year", all_years)

  # no time model:
  fit <- sdmTMB(
    catch_weight ~ 1,
    time = "year",
    offset = log(dogfish$area_swept),
    extra_time = historical_years,
    spatial = "off",
    spatiotemporal = "off",
    data = dogfish,
    mesh = mesh,
    family = delta_gamma()
  )

  set.seed(1)
  expect_message(out <- project(fit, newdata = proj_grid, nsim = 2, uncertainty = "none"), "structures")
  expect_equal(unique(as.numeric(out$est1)), 0.8032636, tolerance = 1e-3)
  expect_message(out <- project(fit, newdata = proj_grid, nsim = 1, uncertainty = "both"), regexp = "random")
  expect_error(out <- project(fit, newdata = proj_grid, nsim = 1, uncertainty = "random"), regexp = "random")

  # no time argument specified errors:
  fit <- sdmTMB(
    catch_weight ~ 1,
    # time = "year",
    offset = log(dogfish$area_swept),
    extra_time = historical_years,
    spatial = "off",
    spatiotemporal = "off",
    data = dogfish,
    # mesh = mesh,
    family = delta_gamma()
  )

  expect_error(out2 <- project(fit, newdata = proj_grid, nsim = 1), regexp = "time")

  # no future time errors:
  fit <- sdmTMB(
    catch_weight ~ 0,
    time = "year",
    offset = log(dogfish$area_swept),
    extra_time = historical_years,
    spatial = "off",
    spatiotemporal = "off",
    time_varying = ~ 1,
    data = dogfish,
    family = tweedie()
  )
  expect_error(out <- project(fit, newdata = subset(proj_grid, year %in% historical_years), nsim = 1), "new")

  # newdata is missing a time step; make sure that's fine and matches not missing the time step:
  all_years <- c(historical_years, 2023:2025)
  proj_grid <- replicate_df(wcvi_grid, "year", all_years)
  set.seed(1)
  out <- project(fit, newdata = proj_grid, nsim = 1)

  all_years <- c(historical_years, c(2023, 2025))
  proj_grid2 <- replicate_df(wcvi_grid, "year", all_years)
  set.seed(1)
  out2 <- project(fit, newdata = proj_grid2, nsim = 1)

  i <- proj_grid$year %in% c(2023, 2025)
  i2 <- proj_grid2$year %in% c(2023, 2025)
  expect_equal(out$est[i,], out2$est[i2,])
})
