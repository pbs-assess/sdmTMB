test_that("A logistic threshold model fits", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  d <- subset(pcod, year >= 2011) # subset for speed
  pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 30)
  m <- sdmTMB(density ~ 0 + as.factor(year) + logistic(depth_scaled), data = d,
    mesh = pcod_spde, family = tweedie(link = "log"),
    time = "year")
  expect_true(all(!is.na(summary(m$sd_report)[,"Std. Error"])))

  expect_true("depth_scaled-s50" %in% tidy(m)$term)
  expect_true("depth_scaled-s95" %in% tidy(m)$term)
  expect_true("depth_scaled-smax" %in% tidy(m)$term)
  expect_equal(tidy(m)[,"estimate"], c(1.555 , 1.655 , 1.718 , 1.138, -0.979, -3.173 , 1.760), tolerance = 1e-3)
})

test_that("A linear threshold model fits", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  d <- subset(pcod, year >= 2011) # subset for speed
  pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 30)
  m <- sdmTMB(density ~ 0 + as.factor(year) + breakpt(depth_scaled), data = d,
    mesh = pcod_spde, family = tweedie(link = "log"), spatial = "off")
  expect_true(all(!is.na(summary(m$sd_report)[,"Std. Error"])))

  expect_true("depth_scaled-slope" %in% tidy(m)$term)
  expect_true("depth_scaled-breakpt" %in% tidy(m)$term)
  expect_equal(tidy(m)[,"estimate"], c(4.798 , 4.779 , 4.768 , 4.112 , 1.085 ,-1.328), tolerance = 1e-3)
})

test_that("A linear threshold *delta* model fits", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")

  set.seed(1)
  predictor_dat <- data.frame(
    X = runif(10000), Y = runif(10000),
    a1 = rnorm(10000)
  )
  mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.2)
  s1 <- sdmTMB_simulate(
    formula = ~ 1 + breakpt(a1),
    data = predictor_dat,
    mesh = mesh,
    family = binomial(),
    range = 0.5,
    phi = 0.02,
    sigma_O = 0.01,
    seed = 42,
    B = 0,
    threshold_coefs = c(0.5, 0)
  )
  s2 <- sdmTMB_simulate(
    formula = ~ 1 + breakpt(a1),
    data = predictor_dat,
    mesh = mesh,
    family = Gamma(link = "log"),
    range = 0.5,
    phi = 1000,
    sigma_O = 0.01,
    seed = 42,
    B = 0,
    threshold_coefs = c(0.3, 0)
  )

  plot(predictor_dat$a1, s1$observed)
  plot(predictor_dat$a1, s2$observed)

  s <- s1
  s$oberved <- s1$observed * s2$observed
  s$a1 <- predictor_dat$a1

  suppressWarnings(
    fit <- sdmTMB(
      observed ~ 1 + breakpt(a1), s, mesh = mesh,
      family = delta_gamma(),
      spatial = "off", control = sdmTMBcontrol(newton_loops = 1L)
    )
  )
  # print(fit) # broken
})
