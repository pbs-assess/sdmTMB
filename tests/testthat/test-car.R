test_that("Forecasting works with a 1D Car model", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  # Good Stan refs, e.g. https://mc-stan.org/users/documentation/case-studies/mbjoseph-CARStan.html
  # demonstrate simple CAR model
  set.seed(123)
  N = 20

  grid <- expand.grid(x = 1:N, y = 1:N)
  n <- nrow(grid)
  distance <- as.matrix(dist(grid))
  W <- array(0, c(n, n))
  W[distance == 1] <- 1

  CAR_nb = W
  diag(CAR_nb) = rowSums(W)
  D = rowSums(W)
  alpha = 0.8
  tau = 2
  B0 = 30

  Tau <- tau*(diag(rowSums(W)) - alpha*W)

  re = mvtnorm::rmvnorm(1, mean = rep(0, nrow(W)), sigma = solve(Tau))

  # spatial data
  spat_data = data.frame(car_region = 1:nrow(W),
                         year = 1,
                         lon = runif(nrow(W)), # dummy
                         lat = runif(nrow(W))) # dummy
  ln_phi = log(0.01)
  df = data.frame(car_region = 1:nrow(W),
                  x = rnorm(nrow(W),3,1),
                  resid = rnorm(nrow(W),0,exp(ln_phi)))
  df = dplyr::left_join(df,spat_data)

  df$y = B0 + re[df$car_region] + df$resid

  # try with limits, no priors
  m <- sdmTMB(y ~ 1, data = df, time = "year",
              mesh = make_mesh(df,c("lon","lat"),n_knots=8),
              spatiotemporal = "off",
              spatial = "on",
              CAR_neighbours = CAR_nb)

  names(B0) = "b_j"
  expect_equal(m$sd_report$par.fixed[which(names(m$sd_report$par.fixed)=="b_j")], B0, tolerance = 1e-2)

  est_ln_tau_inv = as.numeric(m$sd_report$par.fixed[which(names(m$sd_report$par.fixed)=="ln_car_tau_s")])
  est_tau = (1/exp(est_ln_tau_inv))^2
  expect_equal(est_tau, tau, tolerance = 0.1)

  alpha_est = as.numeric(m$sd_report$par.fixed[which(names(m$sd_report$par.fixed)=="logit_car_alpha_s")])
  expect_equal(plogis(alpha_est), alpha, tolerance = 0.07)

  names(ln_phi) = "ln_phi"
  expect_equal(m$sd_report$par.fixed[which(names(m$sd_report$par.fixed)=="ln_phi")],ln_phi, tolerance = 1e-2)

})
