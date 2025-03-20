test_that("A logistic threshold model fits", {
  skip_on_cran()
  d <- subset(pcod, year >= 2011) # subset for speed
  pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 30)
  m <- sdmTMB(density ~ 0 + as.factor(year) + logistic(depth_scaled), data = d,
    mesh = pcod_spde, family = tweedie(link = "log"),
    time = "year")
  expect_true(all(!is.na(summary(m$sd_report)[,"Std. Error"])))

  expect_true("depth_scaled-s50" %in% tidy(m)$term)
  expect_true("depth_scaled-s95" %in% tidy(m)$term)
  expect_true("depth_scaled-smax" %in% tidy(m)$term)
  expect_equal(tidy(m)[,"estimate",drop=TRUE], c(1.555 , 1.655 , 1.718 , 1.138, -0.979, -0.937 , 1.760), tolerance = 1e-3)
})

test_that("A linear threshold model fits", {
  skip_on_cran()
  d <- subset(pcod, year >= 2011) # subset for speed
  pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 30)
  m <- sdmTMB(density ~ 0 + as.factor(year) + breakpt(depth_scaled), data = d,
    mesh = pcod_spde, family = tweedie(link = "log"), spatial = "off")
  expect_true(all(!is.na(summary(m$sd_report)[,"Std. Error"])))

  expect_true("depth_scaled-slope" %in% tidy(m)$term)
  expect_true("depth_scaled-breakpt" %in% tidy(m)$term)
  expect_equal(tidy(m)[,"estimate",drop=TRUE], c(4.798 , 4.779 , 4.768 , 4.112 , 1.085 ,-1.328), tolerance = 1e-3)
})

test_that("A linear threshold *delta* model fits", {
  skip_on_cran()

  set.seed(1)
  predictor_dat <- data.frame(
    X = runif(1000), Y = runif(1000),
    a1 = rnorm(1000)
  )
  mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.2)
  s1 <- sdmTMB_simulate(
    formula = ~ 1 + breakpt(a1),
    data = predictor_dat,
    mesh = mesh,
    family = binomial(),
    range = 0.5,
    phi = 0.001,
    sigma_O = 0.1,
    seed = 42,
    B = 0,
    threshold_coefs = c(0.5, 0.3)
  )
  s2 <- sdmTMB_simulate(
    formula = ~ 1 + breakpt(a1),
    data = predictor_dat,
    mesh = mesh,
    family = Gamma(link = "log"),
    range = 0.5,
    phi = 1000,
    sigma_O = 0.1,
    seed = 42,
    B = 0,
    threshold_coefs = c(0.3, 0.3)
  )

  plot(predictor_dat$a1, s1$observed)
  plot(predictor_dat$a1, s2$observed)

  s <- s1
  s$observed <- s1$observed * s2$observed
  s$a1 <- predictor_dat$a1
  s1$a1 <- predictor_dat$a1
  s2$a1 <- predictor_dat$a1

  ctrl <- sdmTMBcontrol(newton_loops = 1L)

  # binomial works:
  fit1 <- sdmTMB(observed ~ breakpt(a1),
    data = s1,
    family = binomial(),
    # mesh = mesh,
    spatial = "off",
    control = ctrl
  )
  print(fit1)

  s2_pos <- subset(s2, s1$observed > 0)
  # mesh2 <- make_mesh(s2_pos, xy_cols = c("X", "Y"), mesh = mesh$mesh)
  # Gamma works:
  fit2 <- sdmTMB(observed ~ breakpt(a1),
    data = s2_pos,
    family = Gamma(link = "log"),
    spatial = "off",
    # mesh = mesh2,
    control = ctrl
  )
  print(fit2)

  fit <- sdmTMB(
    observed ~ breakpt(a1),
    data = s,
    family = delta_gamma(),
    # mesh = mesh,
    spatial = "off",
    control = ctrl
  )
  print(fit)
  sanity(fit)

  t1 <- tidy(fit1)
  t2 <- tidy(fit2)
  td1 <- tidy(fit, model = 1)
  td2 <- tidy(fit, model = 2)

  expect_equal(t1$estimate, td1$estimate, tolerance = 1e-5)
  expect_equal(t2$estimate, td2$estimate, tolerance = 1e-5)
  expect_equal(t1$std.error, td1$std.error, tolerance = 1e-5)
  expect_equal(t2$std.error, td2$std.error, tolerance = 1e-5)
})

