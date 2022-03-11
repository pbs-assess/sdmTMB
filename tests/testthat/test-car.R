test_that("Test CAR model works with 1 time slice", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  # Good Stan refs, e.g. https://mc-stan.org/users/documentation/case-studies/mbjoseph-CARStan.html
  # demonstrate simple CAR model
  set.seed(123)
  N = 10 # number of cells in each direction

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
  #spat_data = data.frame(car_region = 1:nrow(W),
  #                       year = 1,
  #                       lon = runif(nrow(W)), # dummy
  #                       lat = runif(nrow(W))) # dummy

  # library(ggplot2)
  # ggplot(spat_data, aes(lon, lat)) + geom_point()

  ln_phi = log(0.1)
  df = data.frame(car_region = 1:nrow(W),
                  x = rnorm(nrow(W),3,1),
                  resid = rnorm(nrow(W),0,exp(ln_phi)))
  #df = dplyr::left_join(df,spat_data)

  df$mu = B0 + re[df$car_region] # prediction
  df$y = B0 + re[df$car_region] + df$resid # prediction + observation error
  df$year = 1
  # ggplot(df, aes(lon, lat, colour = mu)) + geom_point() + scale_color_viridis_c()
  # ggplot(df, aes(lon, lat, colour = y)) + geom_point() + scale_color_viridis_c()

  # try with no limits, no priors
  m <- sdmTMB(y ~ 1, data = df,
              spatiotemporal = "off",
              spatial = "on",
              CAR_neighbours = CAR_nb)
  expect_equal(class(m), "sdmTMB")

  # check that adding years doesn't affect the model
  df$year = sample(1:10,size=nrow(df), replace=T)
  m <- sdmTMB(y ~ 1, data = df,time="year",
              spatiotemporal = "off",
              spatial = "on",
              CAR_neighbours = CAR_nb)

  # add a beta prior to car_alpha_s
  m <- sdmTMB(y ~ 1, data = df,time="year",
              spatiotemporal = "off",
              spatial = "on",
              priors=sdmTMBpriors(car_alpha_s = beta(10,1)),
              CAR_neighbours = CAR_nb)
})

test_that("Test spatiotemporal CAR model works", {
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
  # spatial random component
  re_s = mvtnorm::rmvnorm(1, mean = rep(0, nrow(W)), sigma = solve(Tau))
  # spatiotemporal random component
  re_st = mvtnorm::rmvnorm(10, mean = rep(0, nrow(W)), sigma = solve(Tau))

  # spatial data
  spat_data = data.frame(car_region = 1:nrow(W),
                         year = 1,
                         lon = runif(nrow(W)), # dummy
                         lat = runif(nrow(W))) # dummy

  ln_phi = log(0.1)

  # test with several years of data
  df = expand.grid(car_region = 1:nrow(W),
                   year = 1:nrow(re_st))
  df$re_s = re_s[df$car_region]
  df$re_st = re_st[cbind(df$year, df$car_region)]

  df$resid = rnorm(nrow(W),0,exp(ln_phi))

  df$mu = B0 + df$re_s + df$re_st
  df$y = df$mu + df$resid

  #this is blowing up -- will look into why. it's specifically the stuff around 'nldens_st'
  m <- sdmTMB(y ~1, data = df, time="year",
              spatiotemporal = "iid",
              spatial = "off",
              CAR_neighbours = CAR_nb)

})
