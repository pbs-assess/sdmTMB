# basic model fitting and prediction tests

set.seed(42)
x <- runif(50, -1, 1)
y <- runif(50, -1, 1)
N <- length(x)
time_steps <- 20
X <- model.matrix(~1, data.frame(x1 = rnorm(N * time_steps)))
loc <- data.frame(x = x, y = y)
mesh <- make_mesh(loc, xy_cols = c("x", "y"), cutoff = 0.1)
epsilons <- exp(rnorm(time_steps))
s <- sdmTMB_sim(
  x = x, y = y, mesh = mesh, X = X,
  betas = c(0.2), time_steps = time_steps, rho = 0.5,
  phi = 0.2, range = 0.8, sigma_O = 0, sigma_E = epsilons,
  seed = 123, family = gaussian()
)
s$year_centered <- s$time - mean(s$time)
mesh <- make_mesh(s, xy_cols = c("x", "y"), cutoff = 0.1)

test_that("Test that non-stationary model works with random effects in epsilon works", {
  local_edition(2)
  skip_if_not_installed("INLA")

  # fit non-stationary model - iid
  m <- sdmTMB(
    data = s, formula = observed ~ 1,
    time = "time", mesh = mesh,
    epsilon_predictor = "year_centered",
    spatiotemporal = "IID", spatial = "off",
    experimental = list(
      epsilon_predictor = "time",
      epsilon_model = "re"
    ),
    control = sdmTMBcontrol(
      lower = list(b_epsilon = -1, ln_epsilon_re_sigma = -3),
      upper = list(b_epsilon = 1, ln_epsilon_re_sigma = 1)
    )
  )
  idx <- grep("ln_epsilon_re_sigma", names(m$sd_report$value))

  expect_equal(as.numeric(m$sd_report$value[idx]), -1.054972, tolerance = 0.002)

  m <- sdmTMB(
    data = s, formula = observed ~ 1,
    time = "time", mesh = mesh,
    epsilon_predictor = "year_centered",
    spatiotemporal = "AR1", spatial = "off",
    experimental = list(
      epsilon_predictor = "time",
      epsilon_model = "re"
    ),
    control = sdmTMBcontrol(
      lower = list(b_epsilon = -1, ln_epsilon_re_sigma = -3),
      upper = list(b_epsilon = 1, ln_epsilon_re_sigma = 1)
    )
  )
  idx <- grep("ln_epsilon_re_sigma", names(m$sd_report$value))

  expect_equal(as.numeric(m$sd_report$value[idx]), -2.130359, tolerance = 0.002)
})


test_that("Test that non-stationary model works with random effects in epsilon with trend works", {
  local_edition(2)
  skip_if_not_installed("INLA")

  # fit non-stationary model - iid
  m <- sdmTMB(
    data = s, formula = observed ~ 1,
    time = "time", mesh = mesh,
    spatiotemporal = "IID", spatial = "off",
    experimental = list(
      epsilon_predictor = "year_centered",
      epsilon_model = "trend-re"
    ),
    control = sdmTMBcontrol(
      lower = list(b_epsilon = -1, ln_epsilon_re_sigma = -3),
      upper = list(b_epsilon = 1, ln_epsilon_re_sigma = 1)
    )
  )
  idx <- grep("ln_epsilon_re_sigma", names(m$sd_report$value))

  expect_equal(as.numeric(m$sd_report$value[idx]), -2.537232, tolerance = 0.002)

  m <- sdmTMB(
    data = s, formula = observed ~ 1,
    time = "time", mesh = mesh,
    spatiotemporal = "AR1", spatial = "off",
    experimental = list(
      epsilon_predictor = "year_centered",
      epsilon_model = "trend-re"
    ),
    control = sdmTMBcontrol(
      lower = list(b_epsilon = -1, ln_epsilon_re_sigma = -3),
      upper = list(b_epsilon = 1, ln_epsilon_re_sigma = 1)
    )
  )
  idx <- grep("ln_epsilon_re_sigma", names(m$sd_report$value))

  expect_equal(as.numeric(m$sd_report$value[idx]), -2.303735, tolerance = 0.002)
})


test_that("Test that non-stationary model works with in epsilon with trend works", {
  local_edition(2)
  skip_if_not_installed("INLA")

  # fit non-stationary model - iid
  m <- sdmTMB(
    data = s, formula = observed ~ 1,
    time = "time", mesh = mesh,
    spatiotemporal = "IID", spatial = "off",
    experimental = list(
      epsilon_predictor = "year_centered",
      epsilon_model = "trend"
    ),
    control = sdmTMBcontrol(
      lower = list(b_epsilon = -1),
      upper = list(b_epsilon = 1)
    )
  )
  idx <- grep("b_epsilon", names(m$sd_report$value))

  expect_equal(as.numeric(m$sd_report$value[idx]), -0.04554321, tolerance = 0.002)

  m <- sdmTMB(
    data = s, formula = observed ~ 1,
    time = "time", mesh = mesh,
    spatiotemporal = "AR1", spatial = "off",
    experimental = list(
      epsilon_predictor = "year_centered",
      epsilon_model = "trend"
    ),
    control = sdmTMBcontrol(
      lower = list(b_epsilon = -1),
      upper = list(b_epsilon = 1)
    )
  )
  idx <- grep("b_epsilon", names(m$sd_report$value))

  expect_equal(as.numeric(m$sd_report$value[idx]), -0.05376236, tolerance = 0.002)
})
