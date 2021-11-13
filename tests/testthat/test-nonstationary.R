# basic model fitting and prediction tests

set.seed(42)
x <- runif(50, -1, 1)
y <- runif(50, -1, 1)
N <- length(x)
time_steps <- 20
X <- model.matrix(~ 1, data.frame(x1 = rnorm(N * time_steps)))
loc <- data.frame(x = x, y = y)
mesh <- make_mesh(loc, xy_cols = c("x", "y"), cutoff = 0.1)
epsilons = exp(rnorm(time_steps))
s <- sdmTMB_sim(
  x = x, y = y, mesh = mesh, X = X,
  betas = c(0.2), time_steps = time_steps, rho = 0.5,
  phi = 0.2, range = 0.8, sigma_O = 0, sigma_E = epsilons,
  seed = 123, family = gaussian()
)
s$year_centered <- s$time - mean(s$time)

test_that("Test that non-stationary model works with random effects in epsilon works", {
  local_edition(2)
  skip_if_not_installed("INLA")

  # fit non-stationary model - iid
  m <- sdmTMB(
    data = s, formula = observed ~ 1,
    time = "time", mesh = mesh,
    fields = "AR1", include_spatial = FALSE,
    experimental = list(epsilon_predictor = "year_centered",epsilon_model = "re"),
    control = sdmTMBcontrol(lower = list(b_epsilon = -1,ln_epsilon_re_sigma=-3),
                            upper = list(b_epsilon = 1,ln_epsilon_re_sigma=1))
  )
  idx = grep("ln_epsilon_re_sigma",names(m$sd_report$value))

  expect_equal(as.numeric(m$sd_report$value[idx]), -2.13, tolerance = 0.002)

  m <- sdmTMB(
    data = s, formula = observed ~ 1,
    time = "time", mesh = mesh,
    fields = "IID", include_spatial = FALSE,
    experimental = list(epsilon_predictor = "year_centered",epsilon_model = "re"),
    control = sdmTMBcontrol(lower = list(b_epsilon = -1,ln_epsilon_re_sigma=-3),
                            upper = list(b_epsilon = 1,ln_epsilon_re_sigma=1))
  )
  idx = grep("ln_epsilon_re_sigma",names(m$sd_report$value))

  expect_equal(as.numeric(m$sd_report$value[idx]), -1.055, tolerance = 0.002)

})
