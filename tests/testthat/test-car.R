test_that("Forecasting works with a time-varying parameter", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  # Good Stan refs, e.g. https://mc-stan.org/users/documentation/case-studies/mbjoseph-CARStan.html
  # demonstrate simple CAR model
  set.seed(123)
  N = 20

  W = matrix(0, N, N)
  for(i in 1:(N-1)) {
    W[i+1,i] = W[i,i+1]=1
  }
  CAR_nb = W + 2*diag(N)
  CAR_nb[cbind(c(1,N), c(1,N))] = 1
  D = diag(CAR_nb)
  alpha = 0.5
  tau = 0.1
  sigma_inv = tau * (D - alpha*W)
  Sigma = as.matrix(Matrix::forceSymmetric(solve(sigma_inv)))
  re = mvtnorm::rmvnorm(1, mean = rep(0, N), sigma = Sigma)

  # spatial data
  spat_data = data.frame(car_region = 1:N,
                         year = 1,
                         lon = runif(10),
                         lat = runif(10),
                         re = c(re))
  df = data.frame(car_region=sample(1:N,1000,replace=T),
                  x = runif(1000, 0, 5),
                  resid = rnorm(1000, 0, 0.03))
  df = dplyr::left_join(df, spat_data)
  df$y = 3 + 0.1 * df$x + df$re + df$resid

  # try with limits, no priors
  m <- sdmTMB(y ~ 1, data = df, time = "year",
              mesh = make_mesh(df,c("lon","lat"),n_knots=8),
              spatiotemporal = "off",
              spatial = "on",
              CAR_neighbours = CAR_nb,
              control = sdmTMBcontrol(lower = list(car_alpha_s = -1),
                                      upper = list(car_alpha_s = 1)))

  expect_equal(m$sd_report$par.fixed[which(names(m$sd_report$par.fixed)=="ln_car_tau_s")], log(tau), tolerance = 1e-3)
  expect_equal(m$sd_report$par.fixed[which(names(m$sd_report$par.fixed)=="car_alpha_s")], alpha, tolerance = 1e-3)
})
