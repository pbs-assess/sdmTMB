# basic model fitting and prediction tests

#
# test_that("Test that non-stationary model works with random effects in epsilon works", {
#   local_edition(2)
#   skip_on_cran()
#   skip_on_ci()
#   skip_if_not_installed("INLA")
#
#   mesh <- make_mesh(predictor_dat, xy_cols = c("x", "y"), cutoff = 0.1)
#   epsilons <- exp(rnorm(time_steps))
#   s <- sdmTMB_simulate(
#     formula = ~ 1,
#     mesh = mesh, data = predictor_dat,
#     B = c(0.2), rho = 0.5,
#     phi = 0.2, range = 0.8, sigma_O = 0, sigma_E = epsilons[1],
#     seed = 123, family = gaussian()
#   )
#   s$time <- predictor_dat$year
#   s$year_centered <- s$time - mean(s$time)
#
#   # fit non-stationary model - iid
#   m <- sdmTMB(
#     data = s, formula = observed ~ 1,
#     time = "time", mesh = mesh,
#     spatiotemporal = "IID", spatial = "off",
#     experimental = list(
#       epsilon_predictor = "year_centered",
#       epsilon_model = "re"
#     ),
#     control = sdmTMBcontrol(
#       lower = list(b_epsilon = -1, ln_epsilon_re_sigma = -3),
#       upper = list(b_epsilon = 1, ln_epsilon_re_sigma = 1)
#     )
#   )
#   idx <- grep("ln_epsilon_re_sigma", names(m$sd_report$value))
#
#   expect_equal(as.numeric(m$sd_report$value[idx]), -1.054972, tolerance = 0.002)
#
#   m <- sdmTMB(
#     data = s, formula = observed ~ 1,
#     time = "time", mesh = mesh,
#     spatiotemporal = "AR1", spatial = "off",
#     experimental = list(
#       epsilon_predictor = "year_centered",
#       epsilon_model = "re"
#     ),
#     control = sdmTMBcontrol(
#       lower = list(b_epsilon = -1, ln_epsilon_re_sigma = -3),
#       upper = list(b_epsilon = 1, ln_epsilon_re_sigma = 1)
#     )
#   )
#   idx <- grep("ln_epsilon_re_sigma", names(m$sd_report$value))
#
#   expect_equal(as.numeric(m$sd_report$value[idx]), -2.130359, tolerance = 0.002)
# })


## # test_that("Test that non-stationary model works with random effects in epsilon with trend works", {
## #   local_edition(2)
## #   skip_on_cran()
## #   skip_on_ci()
## #   skip_if_not_installed("INLA")
# ## #
#   set.seed(42)
#   mesh <- make_mesh(loc, xy_cols = c("x", "y"), cutoff = 0.1)
#   epsilons <- exp(rnorm(time_steps))
#   s <- sdmTMB_sim(
#     x = x, y = y, mesh = mesh, X = X,
#     betas = c(0.2), time_steps = time_steps, rho = 0.5,
#     phi = 0.2, range = 0.8, sigma_O = 0, sigma_E = epsilons,
#     seed = 123, family = gaussian()
#   )
#   s$year_centered <- s$time - mean(s$time)
#   mesh <- make_mesh(s, xy_cols = c("x", "y"), cutoff = 0.1)

## #   # fit non-stationary model - iid
## #   m <- sdmTMB(
## #     data = s, formula = observed ~ 1,
## #     time = "time", mesh = mesh,
## #     spatiotemporal = "IID", spatial = "off",
## #     experimental = list(
## #       epsilon_predictor = "year_centered",
## #       epsilon_model = "trend-re"
## #     ),
## #     control = sdmTMBcontrol(
## #       lower = list(b_epsilon = -1, ln_epsilon_re_sigma = -3),
## #       upper = list(b_epsilon = 1, ln_epsilon_re_sigma = 1)
## #     )
## #   )
## #   idx <- grep("ln_epsilon_re_sigma", names(m$sd_report$value))
## #
## #   expect_equal(as.numeric(m$sd_report$value[idx]), -2.537232, tolerance = 0.002)
## #   expect_equal(as.numeric(m$sd_report$value[idx]), -2.537232, tolerance = 0.002)
## #
## #   m <- sdmTMB(
## #     data = s, formula = observed ~ 1,
## #     time = "time", mesh = mesh,
## #     spatiotemporal = "AR1", spatial = "off",
## #     experimental = list(
## #       epsilon_predictor = "year_centered",
## #       epsilon_model = "trend-re"
## #     ),
## #     control = sdmTMBcontrol(
## #       lower = list(b_epsilon = -1, ln_epsilon_re_sigma = -3),
## #       upper = list(b_epsilon = 1, ln_epsilon_re_sigma = 1)
## #     )
## #   )
## #   idx <- grep("ln_epsilon_re_sigma", names(m$sd_report$value))
## #
## #   expect_equal(as.numeric(m$sd_report$value[idx]), -2.303735, tolerance = 0.002)
## # })
## #
## #

test_that("Test that non-stationary model works without spatial field and epsilon trend works", {
  local_edition(2)
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")

  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 20)

  pcod$fyear <- as.factor(pcod$year)
  pcod$time <- pcod$year - min(pcod$year) + 1
  pcod$time = scale(pcod$year)
  fit <- sdmTMB(
    density ~ s(depth),
    data = pcod, mesh = mesh,
    spatial="off",
    time = "year",
    spatiotemporal = "ar1",
    family = tweedie(link = "log"),
    experimental = list(epsilon_model = "trend", epsilon_predictor = "time"),
    control = sdmTMBcontrol(lower = list(b_epsilon = -1),
                            upper = list(b_epsilon = 1))
  )
  par <- fit$sd_report$value[which(names(fit$sd_report$value)=="b_epsilon")]
  expect_equal(as.numeric(par), -0.05852822, tolerance = 0.002)

  par <- fit$sd_report$value[which(names(fit$sd_report$value)=="log_sigma_E")]
  expect_equal(as.numeric(par), c(1.0534572, 1.0409799, 1.0285026, 1.0035480, 0.9785934, 0.9536388, 0.9286842, 0.9037296, 0.8787750), tolerance = 0.002)

  # fit non-stationary model - iid
  fit <- sdmTMB(
    density ~ s(depth),
    data = pcod, mesh = mesh,
    spatial="off",
    time = "year",
    spatiotemporal = "iid",
    family = tweedie(link = "log"),
    experimental = list(epsilon_model = "trend", epsilon_predictor = "time"),
    control = sdmTMBcontrol(lower = list(b_epsilon = -1),
                            upper = list(b_epsilon = 1))
  )
  par <- fit$sd_report$value[which(names(fit$sd_report$value)=="b_epsilon")]
  expect_equal(as.numeric(par), -0.04915406, tolerance = 0.002)

  par <- fit$sd_report$value[which(names(fit$sd_report$value)=="log_sigma_E")]
  expect_equal(as.numeric(par), c(1.1262184, 1.1157395, 1.1052607, 1.0843029, 1.0633452, 1.0423874, 1.0214297, 1.0004719, 0.9795142), tolerance = 0.002)
})

test_that("Test that non-stationary model works with spatial field and epsilon trend works", {
  local_edition(2)
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")

  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 20)

  pcod$fyear <- as.factor(pcod$year)
  pcod$time <- pcod$year - min(pcod$year) + 1
  pcod$time = scale(pcod$year)
  fit <- sdmTMB(
    density ~ s(depth),
    data = pcod, mesh = mesh,
    spatial="on",
    time = "year",
    spatiotemporal = "ar1",
    family = tweedie(link = "log"),
    experimental = list(epsilon_model = "trend", epsilon_predictor = "time"),
    control = sdmTMBcontrol(lower = list(b_epsilon = -1),
                            upper = list(b_epsilon = 1))
  )
  par <- fit$sd_report$value[which(names(fit$sd_report$value)=="b_epsilon")]
  expect_equal(as.numeric(par), -0.04435818, tolerance = 0.002)

  par <- fit$sd_report$value[which(names(fit$sd_report$value)=="log_sigma_E")]
  expect_equal(as.numeric(par), c(0.8882275, 0.8787710, 0.8693146, 0.8504016, 0.8314887, 0.8125758, 0.7936628, 0.7747499, 0.7558370), tolerance = 0.002)

  # fit non-stationary model - iid
  fit <- sdmTMB(
    density ~ s(depth),
    data = pcod, mesh = mesh,
    spatial="on",
    time = "year",
    spatiotemporal = "iid",
    family = tweedie(link = "log"),
    experimental = list(epsilon_model = "trend", epsilon_predictor = "time"),
    control = sdmTMBcontrol(lower = list(b_epsilon = -1),
                            upper = list(b_epsilon = 1))
  )
  par <- fit$sd_report$value[which(names(fit$sd_report$value)=="b_epsilon")]
  expect_equal(as.numeric(par), -0.0457674, tolerance = 0.002)

  par <- fit$sd_report$value[which(names(fit$sd_report$value)=="log_sigma_E")]
  expect_equal(as.numeric(par), c(0.8916701, 0.8819133, 0.8721564, 0.8526426, 0.8331288, 0.8136150, 0.7941013, 0.7745875, 0.7550737), tolerance = 0.002)
})


test_that("Test that non-stationary model works with epsilon trend and delta model", {
  local_edition(2)
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")

  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 20)

  pcod$fyear <- as.factor(pcod$year)
  pcod$time <- pcod$year - min(pcod$year) + 1
  pcod$time = scale(pcod$year)
  fit <- sdmTMB(
    density ~ s(depth),
    data = pcod, mesh = mesh,
    spatial="off",
    time = "year",
    spatiotemporal = "ar1",
    family = delta_gamma(),
    experimental = list(epsilon_model = "trend", epsilon_predictor = "time"),
    control = sdmTMBcontrol(lower = list(b_epsilon = -1),
                            upper = list(b_epsilon = 1))
  )

  par <- fit$sd_report$value[which(names(fit$sd_report$value)=="b_epsilon")]
  expect_equal(as.numeric(par), c(-0.07908264, -0.09297464), tolerance = 0.002)

})


test_that("Test that non-stationary model works without spatial field and random effects in epsilon", {
  local_edition(2)
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")

  set.seed(42)
  time_steps <- 20

  epsilons <- exp(rnorm(time_steps, mean = 0, sd = exp(-3)))
  # make fake predictor(s) (a1) and sampling locations:
  predictor_dat <- data.frame(
    X = runif(length(epsilons)*50), Y = runif(length(epsilons)*50),
    a1 = rnorm(length(epsilons)*50), year = rep(1:length(epsilons), each = 50)
  )
  mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)

  sim_dat <- sdmTMB_simulate(
    formula = ~ 1 + a1,
    data = predictor_dat,
    time = "year",
    mesh = mesh,
    family = gaussian(),
    range = 0.5,
    sigma_E = epsilons,
    phi = 0.01,
    sigma_O = 0,
    seed = 42,
    B = c(0.2, -0.4) # B0 = intercept, B1 = a1 slope
  )

  sim_dat$time <- sim_dat$year
  sim_dat$year_centered <- sim_dat$time - mean(sim_dat$time)

  fit <- sdmTMB(
    observed ~ a1,
    data = sim_dat, mesh = mesh,
    spatial="off",
    time = "year",
    spatiotemporal = "iid",
    experimental = list(epsilon_model = "re"),
    control = sdmTMBcontrol(lower = list(ln_epsilon_re_sigma = -20),
                            upper = list(ln_epsilon_re_sigma = -1))
  )

  par <- fit$sd_report$value[which(names(fit$sd_report$value)=="ln_epsilon_re_sigma")]
  expect_equal(as.numeric(par), -13.5, tolerance = 0.002)
  par <- fit$sd_report$par.fixed[1:2]
  expect_equal(as.numeric(par), c(0.2579745,-0.40099), tolerance = 0.002)
})


test_that("Test that non-stationary model works without spatial field and trend and random effects in epsilon", {
  local_edition(2)
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")

  set.seed(42)
  time_steps <- 20

  epsilons <- exp(rnorm(time_steps, mean = 0, sd = exp(-3)))
  # make fake predictor(s) (a1) and sampling locations:
  predictor_dat <- data.frame(
    X = runif(length(epsilons)*50), Y = runif(length(epsilons)*50),
    a1 = rnorm(length(epsilons)*50), year = rep(1:length(epsilons), each = 50)
  )
  mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)

  sim_dat <- sdmTMB_simulate(
    formula = ~ 1 + a1,
    data = predictor_dat,
    time = "year",
    mesh = mesh,
    family = gaussian(),
    range = 0.8,
    rho = 0.5,
    sigma_E = epsilons,
    phi = 0.2,
    sigma_O = 0,
    seed = 42,
    B = c(0.2, -0.4) # B0 = intercept, B1 = a1 slope
  )

  sim_dat$time <- sim_dat$year
  sim_dat$year_centered <- sim_dat$time - min(sim_dat$time)

  fit <- sdmTMB(
    observed ~ a1,
    data = sim_dat, mesh = mesh,
    spatial="off",
    time = "year",
    spatiotemporal = "iid",
    experimental = list(epsilon_model = "trend-re",
                        epsilon_predictor = "year_centered"),
    control = sdmTMBcontrol(lower = list(ln_epsilon_re_sigma = -15, b_epsilon=-1),
                            upper = list(ln_epsilon_re_sigma = -1, b_epsilon=1))
  )

  par <- fit$sd_report$value[which(names(fit$sd_report$value)=="b_epsilon")]
  expect_equal(as.numeric(par), 0.01257052, tolerance = 0.002)

})
