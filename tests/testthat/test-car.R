test_that("Forecasting works with a 1D Car model", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  # Good Stan refs, e.g. https://mc-stan.org/users/documentation/case-studies/mbjoseph-CARStan.html
  # demonstrate simple CAR model
  set.seed(123)
  N = 10

  grid <- expand.grid(x = 1:N, y = 1:N)
  n <- nrow(grid)
  distance <- as.matrix(dist(grid))
  W <- array(0, c(n, n))
  W[distance == 1] <- 1

  #W = matrix(0, N, N)
  #for(i in 1:(N-1)) {
  #  W[i+1,i] = W[i,i+1]=1
  #}
  #CAR_nb = W + 2*diag(N)
  #CAR_nb[cbind(c(1,N), c(1,N))] = 1
  CAR_nb = W
  diag(CAR_nb) = rowSums(W)
  D = rowSums(W)
  alpha = 0.8
  tau = 2
  ln_phi = 0.03
  Tau <- tau*(diag(rowSums(W)) - alpha*W)

  re = mvtnorm::rmvnorm(1, mean = rep(0, nrow(W)), sigma = solve(Tau))

  # spatial data
  spat_data = data.frame(car_region = 1:nrow(W),
                         year = 1,
                         lon = runif(nrow(W)),
                         lat = runif(nrow(W)),
                         x = rnorm(nrow(W),3,1))
  spat_data$y = 30 + spat_data$x*0.1 + c(re) + rnorm(nrow(W), 0, exp(ln_phi))

  # try with limits, no priors
  m <- sdmTMB(y ~ x, data = spat_data, time = "year",
              mesh = make_mesh(spat_data,c("lon","lat"),n_knots=8),
              spatiotemporal = "off",
              spatial = "on",
              CAR_neighbours = CAR_nb,
              control = sdmTMBcontrol(lower = list(car_alpha_s = 0.1),
                                      upper = list(car_alpha_s = 0.9)))

  ln_car_tau_s = log(tau)
  names(ln_car_tau_s) = "ln_car_tau_s"
  expect_equal(m$sd_report$par.fixed[which(names(m$sd_report$par.fixed)=="ln_car_tau_s")], ln_car_tau_s, tolerance = 1e-3)

  names(alpha) = "car_alpha_s"
  expect_equal(m$sd_report$par.fixed[which(names(m$sd_report$par.fixed)=="car_alpha_s")], alpha, tolerance = 1e-2)

  names(ln_phi) = "ln_phi"
  expect_equal(m$sd_report$par.fixed[which(names(m$sd_report$par.fixed)=="ln_phi")],ln_phi, tolerance = 1e-2)

})
