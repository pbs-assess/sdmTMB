context("basic model fitting and prediction tests")

SEED <- 1
set.seed(SEED)
x <- stats::runif(100, -1, 1)
y <- stats::runif(100, -1, 1)
loc <- data.frame(x = x, y = y)
spde <- make_mesh(loc, c("x", "y"), cutoff = 0.02)

test_that("sdmTMB model fit with a covariate beta", {
  initial_betas <- 0.5
  range <- 0.1
  sigma_O <- 0.3 # SD of spatial process
  sigma_E <- 0.3 # SD of spatial process
  phi <- 0.1 # observation error
  s <- sdmTMB_sim(
    x = x, y = y, mesh = spde, time_steps = 6L, betas = initial_betas,
    phi = phi, range = range, sigma_O = sigma_O, sigma_E = sigma_E,
    seed = SEED
  )
  expect_equal(s$observed[c(1:30, 500)],
    c(-0.364045551132498, -0.0884742484390758, -0.748796917796798,
    -0.0926168889732131, 0.0892802633041351, 0.92937309789911, 0.665431954199266,
    -0.81477972686779, 0.474664276994187, 0.929133590564276, -0.263130711477779,
    -0.115523493789406, 0.898574275118016, -0.964967397393866, -1.10224938263034,
    -0.286050035282319, 0.0810888689537039, -0.149727881361854, 0.364458765137563,
    0.31886147691733, -0.184284142634752, 1.27393901104861, -0.822473494638277,
    -0.594882019789678, -0.515872613077436, 0.578924954453508, -0.464615318681477,
    0.476914673336151, -0.52150670383694, -0.0152136800518284, -0.394972244437357
  ), tolerance = 1e-7)
  spde <- make_mesh(s, c("x", "y"), cutoff = 0.02)
  plot(spde)
  m <- sdmTMB(data = s, formula = observed ~ 0 + cov1, time = "time", spde = spde)
  expect_output(print(m), "fit by")
  expect_output(summary(m), "fit by")
  p <- as.list(m$model$par)
  r <- m$tmb_obj$report()
  est <- tidy(m, "ran_pars")
  expect_equal(sort(est[,"estimate", drop = TRUE]),
    sort(c(0.106559618205922, 0.0646179753599145, 0.303614800474715, 0.308394406301974)),
    tolerance = 1e-7)
  expect_equal(m$model$convergence, 0L)
  expect_equal((p$b_j - initial_betas)^2, 0, tol = 0.001)
  expect_equal((exp(p$ln_phi) - phi)^2, 0, tol = 0.002)
  expect_equal((r$sigma_O - sigma_O)^2, 0, tol = 0.002)
  expect_equal((r$sigma_E[1] - sigma_E)^2, 0, tol = 0.001)
  expect_equal(est$estimate[est$term == "range"], range, tol = 0.01)
  p <- predict(m)
  r <- residuals(m)
  expect_equal(mean((p$est - s$observed)^2), 0, tol = 0.002)
})

test_that("Anisotropy fits and plots", {
  skip_on_ci()
  skip_on_cran()
  m <- sdmTMB(data = pcod,
    formula = density ~ 0 + as.factor(year),
    spde = make_mesh(pcod, c("X", "Y"), n_knots = 50, type = "kmeans"),
    family = tweedie(link = "log"), anisotropy = TRUE)
  expect_identical(class(m), "sdmTMB")
  plot_anisotropy(m)
  expect_equal(m$tmb_obj$report()$H,
    structure(
      c(0.665528444798002, 0.079350716881963, 0.079350716881963, 1.51202633656794),
      .Dim = c(2L, 2L)),
    tolerance = 1e-3)
})

test_that("Regularization works", {
  skip_on_cran()
  skip_on_ci()
  d <- subset(pcod, year >= 2015)
  d$depth_scaled <- as.numeric(scale(d$depth_scaled))
  m1 <- sdmTMB(data = d,
    formula = density ~ 0 + depth_scaled + as.factor(year),
    spde = make_mesh(d, c("X", "Y"), n_knots = 50, type = "kmeans"),
    family = tweedie(link = "log"))

  # Bypassing via NAs:
  m2 <- sdmTMB(data = d,
    formula = density ~ 0 + depth_scaled + as.factor(year),
    spde = make_mesh(d, c("X", "Y"), n_knots = 50, type = "kmeans"),
    family = tweedie(link = "log"),
    penalties = c(NA, NA, NA))
  expect_equal(m1$sd_report, m2$sd_report)

  # Ridge regression on depth term:
  m2 <- sdmTMB(data = d,
    formula = density ~ 0 + depth_scaled + as.factor(year),
    spde = make_mesh(d, c("X", "Y"), n_knots = 50, type = "kmeans"),
    family = tweedie(link = "log"),
    penalties = c(1, NA, NA))
  b1 <- tidy(m1)
  b2 <- tidy(m2)
  expect_lt(b2$estimate[1], b1$estimate[2])
  expect_lt(b2$std.error[1], b1$std.error[2])
})

test_that("A model with splines works", {
  skip_on_cran()
  skip_on_ci()
  d <- subset(pcod, year >= 2015)
  m <- sdmTMB(data = d,
    formula = density ~ 1 + s(depth_scaled),
    spde = make_mesh(d, c("X", "Y"), n_knots = 50, type = "kmeans"),
    family = tweedie(link = "log"))
  expect_identical(class(m), "sdmTMB")
})

test_that("A spatiotemporal version works with predictions on new data points", {
  d <- subset(pcod, year >= 2015)
  pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 30)
  m <- sdmTMB(
    data = d,
    formula = density ~ 0 + as.factor(year),
    time = "year", spde = pcod_spde, family = tweedie(link = "log"),
    include_spatial = FALSE
  )
  # Predictions at original data locations:
  predictions <- predict(m)
  predictions$resids <- residuals(m) # randomized quantile residuals
  # Predictions onto new data:
  predictions <- predict(m, newdata = subset(qcs_grid, year >= 2015))
  expect_identical(class(predictions), "data.frame")
})

test_that("Predictions on the original data set as `newdata`` return the same predictions", {
  skip_on_cran()
  skip_on_ci()
  set.seed(1)
  x <- stats::runif(70, -1, 1)
  y <- stats::runif(70, -1, 1)
  loc <- data.frame(x = x, y = y)
  spde <- make_mesh(loc, c("x", "y"), cutoff = 0.02)
  dat <- sdmTMB_sim(x = x, y = y, mesh = spde,
    time_steps = 9, rho = 0, range = 0.2,
    sigma_O = 0, sigma_E = 0.3, phi = 0.1,
    seed = 1
  )
  spde <- make_mesh(dat, c("x", "y"), cutoff = 0.02)
  m <- sdmTMB(
    ar1_fields = FALSE, include_spatial = FALSE,
    data = dat, formula = observed ~ 1, time = "time",
    family = gaussian(link = "identity"), spde = spde
  )
  p <- predict(m)
  p_nd <- predict(m, newdata = dat)
  tidy(m)
  tidy(m, conf.int = TRUE)
  tidy(m, effects = "ran_par")
  tidy(m, effects = "ran_par", conf.int = TRUE)

  cols <- c("est", "est_non_rf", "est_rf", "omega_s", "epsilon_st")
  expect_equal(p[,cols], p_nd[,cols], tolerance = 1e-3)

  m <- sdmTMB(
    ar1_fields = TRUE, include_spatial = FALSE,
    data = dat, formula = observed ~ 1, time = "time",
    family = gaussian(link = "identity"), spde = spde
  )

  p <- predict(m)
  p_nd <- predict(m, newdata = dat)
  expect_equal(p[,cols], p_nd[,cols], tolerance = 1e-3)
})

test_that("A time-varying model fits and predicts appropriately", {
  skip_on_cran()
  skip_on_ci()
  SEED <- 42
  set.seed(SEED)
  x <- stats::runif(60, -1, 1)
  y <- stats::runif(60, -1, 1)
  initial_betas <- 0.5
  range <- 0.5
  sigma_O <- 0
  sigma_E <- 0.1
  phi <- 0.1
  sigma_V <- 0.3
  loc <- data.frame(x = x, y = y)
  spde <- make_mesh(loc, c("x", "y"), cutoff = 0.02)

  s <- sdmTMB_sim(
    x = x, y = y, mesh = spde, range = range,
    betas = initial_betas, time_steps = 12L, sigma_V = sigma_V,
    phi = phi, sigma_O = sigma_O, sigma_E = sigma_E,
    seed = SEED
  )
  spde <- make_mesh(s, c("x", "y"), cutoff = 0.02)
  m <- sdmTMB(data = s, formula = observed ~ 0, include_spatial = FALSE,
    time_varying = ~ 0 + cov1, time = "time", spde = spde, mgcv = FALSE)
  expect_equal(exp(m$model$par["ln_tau_V"])[[1]], sigma_V, tolerance = 0.05)
  tidy(m, effects = "ran_par")
  b_t <- dplyr::group_by(s, time) %>%
    dplyr::summarize(b_t = unique(b), .groups = "drop") %>%
    dplyr::pull(b_t)
  r <- m$tmb_obj$report()
  b_t_fit <- r$b_rw_t
  plot(b_t, b_t_fit, asp = 1);abline(a = 0, b = 1)
  expect_equal(mean((b_t- b_t_fit)^2), 0, tolerance = 1e-4)
  p <- predict(m)
  plot(p$est, s$observed, asp = 1);abline(a = 0, b = 1)
  expect_equal(mean((p$est - s$observed)^2), 0, tolerance = 0.01)

  cols <- c("est", "est_non_rf", "est_rf", "omega_s", "epsilon_st")
  p_nd <- predict(m, newdata = s)
  expect_equal(p[,cols], p_nd[,cols], tolerance = 1e-4)
})

test_that("Year indexes get created correctly", {
  expect_identical(make_year_i(c(1, 2, 3)),    c(0L, 1L, 2L))
  expect_identical(make_year_i(c(1L, 2L, 3L)), c(0L, 1L, 2L))
  expect_identical(make_year_i(c(1L, 2L, 4L)), c(0L, 1L, 2L))
  expect_identical(make_year_i(c(1, 2, 4, 2)), c(0L, 1L, 2L, 1L))
  expect_identical(make_year_i(c(1, 4, 2)),    c(0L, 2L, 1L))
})

test_that("A spatial trend model fits", {
  skip_on_cran()
  skip_on_ci()
  d <- subset(pcod, year >= 2011) # subset for speed
  pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 30)
  m <- sdmTMB(density ~ depth_scaled, data = d,
    spde = pcod_spde, family = tweedie(link = "log"),
    spatial_trend = TRUE, time = "year")
  expect_true(all(!is.na(summary(m$sd_report)[,"Std. Error"])))
})

test_that("A logistic threshold model fits", {
  skip_on_cran()
  skip_on_ci()
  d <- subset(pcod, year >= 2011) # subset for speed
  pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 30)
  m <- sdmTMB(density ~ 0 + as.factor(year) + logistic(depth_scaled), data = d,
    spde = pcod_spde, family = tweedie(link = "log"),
    spatial_trend = FALSE, time = "year")
  expect_true(all(!is.na(summary(m$sd_report)[,"Std. Error"])))
})

test_that("A linear threshold model fits", {
  skip_on_cran()
  skip_on_ci()
  d <- subset(pcod, year >= 2011) # subset for speed
  pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 30)
  m <- sdmTMB(density ~ 0 + as.factor(year) + breakpt(depth_scaled), data = d,
    spde = pcod_spde, family = tweedie(link = "log"),
    spatial_trend = FALSE, time = "year")
  expect_true(all(!is.na(summary(m$sd_report)[,"Std. Error"])))
})

test_that("SPDE as generated by make_mesh is consistent", {
  set.seed(42)
  skip_on_cran()
  skip_on_ci()
  set.seed(1)
  x <- stats::runif(70, -1, 1)
  y <- stats::runif(70, -1, 1)
  loc <- data.frame(x = x, y = y)
  spde <- make_mesh(loc, c("x", "y"), cutoff = 0.02)

  # spot check first 10 locations
  loc_xy <- matrix(c(
    -0.468982674,-0.321854124,
    -0.255752201,0.678880700,
    0.145706727,-0.306633022,
    0.816415580,-0.332450138,
    -0.596636138,-0.047297510,
    0.796779370,0.784396672,
    0.889350537,0.728678941,
    0.321595585,-0.220020913,
    0.258228088,0.554641398,
    -0.876427459,0.921235994), ncol = 2, byrow = TRUE)

  expect_equal(c(spde$loc_xy[1:10,]), c(loc_xy[1:10,]), tolerance = 1e-8)

  # test the n
  expect_equal(spde$spde$param.inla$n, 141.0, tolerance = 1e-8)

  # test M0 values - first 10
  expect_equal(c(0.06889138, 0.04169878, 0.01547313, 0.05734409, 0.03271476, 0.03672571, 0.05200510, 0.01293629, 0.03631651, 0.04446844),
    as.numeric(spde$spde$param.inla$M0)[c(1,143,285,427,569,711,853,995,1137,1279)], tolerance=1e-4)
  # test M1 values - first 10
  expect_equal(c(3.602458, 5.105182, 4.571970, 4.353246, 4.496775, 4.802626, 4.505773, 4.935963, 4.419817, 5.169096),
    as.numeric(spde$spde$param.inla$M1)[c(1,143,285,427,569,711,853,995,1137,1279)], tolerance=1e-4)
  # test M2 values - first 10
  expect_equal(c(253.1738,814.3858,1754.2543,401.3752,865.0118,754.9644,543.8968,2425.8404,729.7511,1014.0049),
    as.numeric(spde$spde$param.inla$M2)[c(1,143,285,427,569,711,853,995,1137,1279)], tolerance=1e-4)

})

# test_that("The `map` argument works.", {
#   skip_on_cran()
#   skip_on_ci()
#   d <- subset(pcod, year >= 2013)
#   pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 50)
#
#   m <- sdmTMB(density ~ 0 + as.factor(year),
#     data = d, time = "year", spde = pcod_spde,
#     family = tweedie(link = "log")
#   )
#   p <- predict(m)
#
#   .map <- m$tmb_map
#   .map$ln_tau_O <- as.factor(NA)
#   .map$ln_tau_E <- as.factor(NA)
#   .map$omega_s <- factor(rep(NA, length(m$tmb_params$omega_s)))
#   .map$epsilon_st <- factor(rep(NA, length(m$tmb_params$epsilon_st)))
#   .map$ln_kappa <- as.factor(NA)
#
#   m.map <- sdmTMB(density ~ 0 + as.factor(year),
#     data = d, time = "year", spde = pcod_spde,
#     family = tweedie(link = "log"), map = .map
#   )
#   p.map <- predict(m.map)
#
#   expect_true(all(p.map$omega_s == 0))
#   expect_true(all(p.map$epsilon_st == 0))
#   expect_true(!identical(p$est, p.map$est))
#   })

test_that("The `map_rf` argument works.", {
  skip_on_cran()
  skip_on_ci()
  d <- subset(pcod, year >= 2013)
  pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 50)
  m <- sdmTMB(density ~ 0 + as.factor(year),
    data = d, time = "year", spde = pcod_spde,
    family = tweedie(link = "log")
  )
  p <- predict(m)
  m.map <- sdmTMB(density ~ 0 + as.factor(year),
    data = d, time = "year", spde = pcod_spde,
    family = tweedie(link = "log"), map_rf = TRUE
  )
  p.map <- predict(m.map)

  expect_true(all(p.map$omega_s == 0))
  expect_true(all(p.map$epsilon_st == 0))
  expect_true(!identical(p$est, p.map$est))
  expect_true(length(unique(p.map$est)) == length(unique(d$year)))

  # basic lm test:
  dpos <- d[d$density > 0, ]
  pcod_spde <- make_mesh(dpos, c("X", "Y"), cutoff = 50)
  m.sdmTMB.map <- sdmTMB(log(density) ~ depth, data = dpos,
    family = gaussian(), map_rf = TRUE, spde = pcod_spde)
  m.stats.glm <- stats::lm(log(density) ~ depth, data = dpos)
  m.glmmTMB <- glmmTMB::glmmTMB(log(density) ~ depth, data = dpos)
  .t <- tidy(m.sdmTMB.map)
  expect_equal(.t$estimate, as.numeric(coef(m.stats.glm)), tolerance = 1e-5)
  expect_equal(exp(m.sdmTMB.map$model$par[["ln_phi"]]),
    glmmTMB::sigma(m.glmmTMB), tolerance = 1e-5)

  # Bernoulli:
  pcod_binom <- pcod
  pcod_binom$present <- ifelse(pcod_binom$density > 0, 1L, 0L)
  pcod_spde <- make_mesh(pcod_binom, c("X", "Y"), cutoff = 50)
  m.sdmTMB.map <- sdmTMB(present ~ depth, data = pcod_binom,
    family = binomial(link = "logit"), map_rf = TRUE, spde = pcod_spde)
  m.stats.glm <- stats::glm(present ~ depth, data = pcod_binom,
    family = binomial(link = "logit"))
  m.glmmTMB <- glmmTMB::glmmTMB(present ~ depth, data = pcod_binom,
    family = binomial(link = "logit"))
  b1 <- as.numeric(unclass(glmmTMB::fixef(m.glmmTMB))$cond)
  b2 <- tidy(m.sdmTMB.map)$estimate
  b3 <- as.numeric(coef(m.stats.glm))
  expect_equal(b1, b2, tolerance = 1e-7)
  expect_equal(b3, b2, tolerance = 1e-6)
})

test_that("Large coordinates cause a warning.", {
  d <- subset(pcod, year == 2017)
  d$X <- d$X * 1000
  d$Y <- d$Y * 1000
  expect_warning(pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 70))
})
