test_that("Test CAR model works with 1 time slice", {
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

  CAR_nb = W
  diag(CAR_nb) = rowSums(W)
  D = rowSums(W)
  alpha = 0.8
  sigma = 0.5
  tau <- 1/(sigma*sigma)
  B0 = 0.5

  Tau <- tau*(diag(rowSums(W)) - alpha*W)

  re = mvtnorm::rmvnorm(1, mean = rep(0, nrow(W)), sigma = solve(Tau))

  # spatial data
  spat_data = data.frame(car_region = 1:nrow(W),
                         year = 1,
                         lon = runif(nrow(W)), # dummy
                         lat = runif(nrow(W))) # dummy

  # library(ggplot2)
  # ggplot(spat_data, aes(lon, lat)) + geom_point()

  ln_phi = log(0.1)
  df = data.frame(car_region = 1:nrow(W),
                  x = rnorm(nrow(W),3,1),
                  resid = rnorm(nrow(W),0,exp(ln_phi)))
  df = dplyr::left_join(df,spat_data)

  df$mu = B0 + re[df$car_region]
  df$y = B0 + re[df$car_region] + df$resid

  # ggplot(df, aes(lon, lat, colour = mu)) + geom_point() + scale_color_viridis_c()
  # ggplot(df, aes(lon, lat, colour = y)) + geom_point() + scale_color_viridis_c()

  # try with limits, no priors
  m <- sdmTMB(y ~ 1, data = df, time = "year",
              spatiotemporal = "off",
              spatial = "on",
              CAR_neighbours = CAR_nb)
  expect_equal(class(m), "sdmTMB")

  expect_equal(class(predict(m, df)), "data.frame")

  # test with several years of data
  df = expand.grid(car_region = 1:nrow(W),
                  year = 1:5)
  #df = dplyr::left_join(df,spat_data)
  df$resid = rnorm(nrow(W),0,exp(-1.203973))
  df$mu = B0 + re[df$car_region]
  df$y = B0 + re[df$car_region] + df$resid

  #spde = make_mesh(df, c("car_region","car_region"),n_knots=4)
  m <- sdmTMB(y ~1, data = df, time = "year",
              spatiotemporal = "off",
              spatial = "on",
              CAR_neighbours = CAR_nb,
              control = sdmTMBcontrol(lower = list(logit_car_alpha_s = -10, ln_car_tau_s = -5),
                                      upper = list(logit_car_alpha_s = 10, ln_car_tau_s = 4)))

})

test_that("Test CAR model works with several time slices", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")

  set.seed(123)
  N = 10

  grid <- expand.grid(x = 1:N, y = 1:N)
  n <- nrow(grid)
  distance <- as.matrix(dist(grid))
  W <- array(0, c(n, n))
  W[distance == 1] <- 1

  CAR_nb = W
  diag(CAR_nb) = rowSums(W)
  D = rowSums(W)
  alpha = 0.8
  sigma = 0.5
  tau <- 1/(sigma*sigma)
  B0 = 0.5

  Tau <- tau*(diag(rowSums(W)) - alpha*W)

  re = mvtnorm::rmvnorm(1, mean = rep(0, nrow(W)), sigma = solve(Tau))

  # spatial data
  spat_data = data.frame(car_region = 1:nrow(W),
                         year = 1,
                         lon = runif(nrow(W)), # dummy
                         lat = runif(nrow(W))) # dummy

  # library(ggplot2)
  # ggplot(spat_data, aes(lon, lat)) + geom_point()

  ln_phi = log(0.1)

  # test with several years of data
  df = expand.grid(car_region = 1:nrow(W),
                   year = 1:5)
  #df = dplyr::left_join(df,spat_data)
  df$resid = rnorm(nrow(W),0,ln_phi)
  df$mu = B0 + re[df$car_region]
  df$y = B0 + re[df$car_region] + df$resid

  #spde = make_mesh(df, c("car_region","car_region"),n_knots=4)
  m <- sdmTMB(y ~1, data = df, time = "year",
              spatiotemporal = "off",
              spatial = "on",
              CAR_neighbours = CAR_nb,
              control = sdmTMBcontrol(lower = list(logit_car_alpha_s = -10, ln_car_tau_s = -5),
                                      upper = list(logit_car_alpha_s = 10, ln_car_tau_s = 4)))

})


test_that("Test CAR model works with several time slices", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")

  set.seed(123)
  N = 10

  grid <- expand.grid(x = 1:N, y = 1:N)
  n <- nrow(grid)
  distance <- as.matrix(dist(grid))
  W <- array(0, c(n, n))
  W[distance == 1] <- 1

  CAR_nb = W
  diag(CAR_nb) = rowSums(W)
  D = rowSums(W)
  alpha = 0.8
  sigma = 0.5
  tau <- 1/(sigma*sigma)
  B0 = 0.5

  Tau <- tau*(diag(rowSums(W)) - alpha*W)

  re = mvtnorm::rmvnorm(1, mean = rep(0, nrow(W)), sigma = solve(Tau))

  # spatial data
  spat_data = data.frame(car_region = 1:nrow(W),
                         year = 1,
                         lon = runif(nrow(W)), # dummy
                         lat = runif(nrow(W))) # dummy

  # library(ggplot2)
  # ggplot(spat_data, aes(lon, lat)) + geom_point()

  ln_phi = log(0.1)

  # test with several years of data
  df = expand.grid(car_region = 1:nrow(W),
                   year = 1:5)
  #df = dplyr::left_join(df,spat_data)
  df$resid = rnorm(nrow(W),0,ln_phi)
  df$mu = B0 + re[df$car_region]
  df$y = B0 + re[df$car_region] + df$resid

  #spde = make_mesh(df, c("car_region","car_region"),n_knots=4)
  m <- sdmTMB(y ~1, data = df, time = "year",
              spatiotemporal = "iid",
              spatial = "off",
              CAR_neighbours = CAR_nb,
              control = sdmTMBcontrol(multiphase=FALSE,
                                      lower = list(logit_car_alpha_s = -10, ln_car_tau_s = -5),
                                      upper = list(logit_car_alpha_s = 10, ln_car_tau_s = 4)))

})

