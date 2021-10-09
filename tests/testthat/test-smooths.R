test_that("A model with 2 s() splines works", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  d <- subset(pcod, year >= 2000 & density > 0)
  pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 30)
  m <- sdmTMB(
    data = d,
    formula = log(density) ~ s(depth_scaled) + s(year, k = 5),
    spde = pcod_spde,
    control = sdmTMBcontrol(map_rf = TRUE)
  )
  expect_equal(ncol(m$tmb_data$X_ij), 1L)
  expect_equal(length(m$tmb_data$Zs), 2L)
  # head(m$tmb_data$Zs[[1]])
  # head(m$tmb_data$Zs[[2]])
  expect_equal(ncol(m$tmb_data$Xs), 2L)
  expect_equal(m$tmb_data$b_smooth_start, c(0L, 8L))
  # exp(as.list(m$sd_report, "Estimate")$ln_smooth_sigma)
  # as.list(m$sd_report, "Estimate")$b_smooth
  expect_equal(sum(is.na(as.list(m$sd_report, "Std. Error")$b_smooth)), 0L)
  p <- predict(m, newdata = NULL)
  m_mgcv <- mgcv::gam(log(density) ~ s(depth_scaled) + s(year, k = 5), data = d)
  p_mgcv <- predict(m_mgcv, newdata = d)
  plot(p$est, p_mgcv)
  abline(a = 0, b = 1)
  expect_gt(cor(p$est, p_mgcv), 0.999)
})

test_that("A model with t2() spline works", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  set.seed(2938)
  dat <- mgcv::gamSim(1, n = 400, dist = "normal", scale = 2)
  dat$X <- runif(nrow(dat))
  dat$Y <- runif(nrow(dat))
  spde <- make_mesh(dat, c("X", "Y"), cutoff = 0.1)
  m_mgcv <- mgcv::gam(y ~ t2(x0, x1, k = 7) + s(x2),
    data = dat,
    method = "REML"
  )
  p_mgcv <- predict(m_mgcv)
  m <- sdmTMB(y ~ t2(x0, x1, k = 7) + s(x2),
    data = dat,
    spde = spde, control = sdmTMBcontrol(map_rf = TRUE)
  )
  p <- predict(m, newdata = NULL)
  plot(p$est, p_mgcv)
  abline(a = 0, b = 1)
  expect_gt(cor(p$est, p_mgcv), 0.9999)
  expect_equal(as.numeric(p$est), as.numeric(p_mgcv), tolerance = 0.001)
})

test_that("A model with by in spline works", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  set.seed(19203)
  # examples from ?mgcv::gam.models
  # continuous by example:
  dat <- mgcv::gamSim(3, n = 800)
  m_mgcv <- mgcv::gam(y ~ s(x2, by = x1), data = dat)
  p_mgcv <- predict(m_mgcv)
  dat$X <- runif(nrow(dat))
  dat$Y <- runif(nrow(dat))
  spde <- make_mesh(dat, c("X", "Y"), cutoff = 0.1)
  m <- sdmTMB(y ~ s(x2, by = x1),
    data = dat,
    spde = spde, control = sdmTMBcontrol(map_rf = TRUE)
  )
  p <- predict(m, newdata = NULL)
  plot(p$est, p_mgcv)
  abline(a = 0, b = 1)
  expect_gt(cor(p$est, p_mgcv), 0.9999)

  # Factor `by' variable example (with a spurious covariate x0)
  set.seed(1)
  dat <- mgcv::gamSim(4)
  m_mgcv <- mgcv::gam(y ~ fac + s(x2, by = fac) + s(x0), data = dat)
  p_mgcv <- predict(m_mgcv)
  dat$X <- runif(nrow(dat))
  dat$Y <- runif(nrow(dat))
  spde <- make_mesh(dat, c("X", "Y"), cutoff = 0.1)
  m <- sdmTMB(y ~ fac + s(x2, by = fac) + s(x0),
    data = dat,
    spde = spde, control = sdmTMBcontrol(map_rf = TRUE)
  )
  p <- predict(m, newdata = NULL)
  plot(p$est, p_mgcv)
  abline(a = 0, b = 1)
  expect_gt(cor(p$est, p_mgcv), 0.9999)
})
