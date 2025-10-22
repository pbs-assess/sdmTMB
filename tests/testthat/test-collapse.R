test_that("collapsing spatial and spatiotemporal fields works", {
  set.seed(123)

  predictor_dat <- data.frame(
    X = runif(1000), Y = runif(1000),
    a1 = rnorm(1000), year = sample(1:5, size = 1000, replace = TRUE)
  )
  mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)

  sim_dat <- sdmTMB_simulate(
    formula = ~ 1 + a1,
    data = predictor_dat,
    time = "year",
    mesh = mesh,
    family = tweedie(),
    range = 2,
    sigma_E = 0,
    phi = 0.2,
    sigma_O = 0,
    seed = 42,
    B = c(0.2, -0.4) # B0 = intercept, B1 = a1 slope
  )
  # create some fake 0s
  sim_dat$observed[sample(1:500, size = 50, replace = FALSE)] <- 0

  # test spatial collapse with gaussian family
  fit_nospatial <- sdmTMB(observed ~ a1, data = sim_dat, mesh = mesh, time = "year",
                          spatiotemporal = "off", spatial="off")
  fit <- sdmTMB(observed ~ a1, data = sim_dat, mesh = mesh, time = "year",
                spatiotemporal = "off",
                control = sdmTMBcontrol(collapse_spatial_variance = TRUE))
  expect_equal(tidy(fit_nospatial), tidy(fit))
  expect_equal(tidy(fit_nospatial,"ran_pars"), tidy(fit,"ran_pars"))
  expect_equal(tidy(fit_nospatial,"ran_vals"), tidy(fit,"ran_vals"))

  # test spatial / spatiotemporal collapse with spatiotemporal on
  fit_nospatial <- sdmTMB(observed ~ a1, data = sim_dat, mesh = mesh, time = "year",
                          spatiotemporal = "off", spatial="off")
  fit <- sdmTMB(observed ~ a1, data = sim_dat, mesh = mesh, time = "year",
                control = sdmTMBcontrol(collapse_spatial_variance = TRUE))
  expect_equal(tidy(fit_nospatial), tidy(fit))
  expect_equal(tidy(fit_nospatial,"ran_pars"), tidy(fit,"ran_pars"))
  expect_equal(tidy(fit_nospatial,"ran_vals"), tidy(fit,"ran_vals"))

  # test spatial collapse with delta family
  fit_nospatial <- sdmTMB(observed ~ a1, data = sim_dat, mesh = mesh, time = "year",
                          spatiotemporal = "off", spatial="off", family=delta_gamma())
  fit <- sdmTMB(observed ~ a1, data = sim_dat, mesh = mesh, time = "year",
                spatiotemporal = "off", family=delta_gamma(),
                control = sdmTMBcontrol(collapse_spatial_variance = TRUE))
  expect_equal(tidy(fit_nospatial), tidy(fit))
  expect_equal(tidy(fit_nospatial,"ran_pars"), tidy(fit,"ran_pars"))
  expect_equal(tidy(fit_nospatial,"ran_vals"), tidy(fit,"ran_vals"))

  # test spatial / spatiotemporal collapse with delta family
  # fit_nospatial <- sdmTMB(observed ~ a1, data = sim_dat, mesh = mesh, time = "year",
  #                         spatiotemporal = "off", spatial="off", family=delta_gamma())
  # fit <- sdmTMB(observed ~ a1, data = sim_dat, mesh = mesh, time = "year",
  #               family=delta_gamma())
  # expect_equal(tidy(fit_nospatial), tidy(fit))
  # expect_equal(tidy(fit_nospatial,"ran_pars"), tidy(fit,"ran_pars"))
  # expect_equal(tidy(fit_nospatial,"ran_vals"), tidy(fit,"ran_vals"))
})
