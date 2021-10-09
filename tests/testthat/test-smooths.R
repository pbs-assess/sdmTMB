# basic model fitting and prediction tests

test_that("A model without splines works", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")

  # test regular model w/o spline works
  d <- subset(pcod, year >= 2015)
  pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 25) # a coarse mesh for example speed
  m <- sdmTMB(data = d,
    formula = density ~ depth_scaled,
    spde = pcod_spde,
    family = tweedie(link = "log"),
    control = sdmTMBcontrol(map_rf = TRUE))
  m
  expect_identical(class(m), "sdmTMB")

})

test_that("A model with s() spline works", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")

  # test regular model with s() works
  d <- subset(pcod, year >= 2015 & density > 0)

  pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 25)
  m <- sdmTMB(data = d,
    formula = log(density) ~ s(depth_scaled),
    spde = pcod_spde,
    control = sdmTMBcontrol(map_rf = TRUE))

  head(m$tmb_data$X_ij)
  head(m$tmb_data$Zs[[1]])
  head(m$tmb_data$Xs)
  m$tmb_data$has_smooths
  m$tmb_data$smooth_matrix_dims

  m$sd_report
  exp(as.list(m$sd_report, "Estimate")$ln_smooth_sigma)
  as.list(m$sd_report, "Estimate")$b_smooth

  p1 <- predict(m, newdata = NULL)

  m_mgcv <- mgcv::gam(log(density) ~ s(depth_scaled), data = d)
  p2 <- predict(m_mgcv, newdata = d)

  plot(p1$est, p2);abline(a = 0, b = 1)

  m_gamm4 <- gamm4::gamm4(log(density) ~ s(depth_scaled), data = d)
  p2 <- predict(m_gamm4$gam)
  plot(p1$est, p2);abline(a = 0, b = 1)

  # str(m_gamm4$mer)
  m_gamm4$mer@beta

  # m_gamm4$mer
  # m_gamm4$gam$smooth

  # m_brms2 <- brms::brm(log(density) ~ year + s(depth_scaled), data = d,
  #   chains = 1, iter = 100, control = list(adapt_delta = 0.8))
  # m_brms2$model

  m_brms <- brms::brm(log(density) ~ s(depth_scaled), data = d,
    chains = 4, control = list(adapt_delta = 0.98), cores = 4, iter = 1000)

  p <- rstan::extract(m_brms$fit)
  mean(p$bs)
  z <- apply(p$s_1_1, 2, median)
  median(p$sds_1_1)

  dd <- brms::standata(m_brms)
  head(dd$Xs)
  # plot(dd$Xs)
  head(dd$Zs_1_1)

  head(m$tmb_data$Xsm[[1]])
  head(dd$X)

  sort(apply(dd$Zs_1_1, 2, sd))
  sort(apply(m$tmb_data$Zs[[1]], 2, sd))

  a <- predict(m_brms)
  plot(p1$est, a[,1]);abline(a = 0, b = 1)
  plot(p2, a[,1]);abline(a = 0, b = 1)
  plot(p2,p1$est);abline(a = 0, b = 1)

  plot(sort(z), sort(as.list(m$sd_report, "Estimate")$b_smooth));abline(a = 0, b = 1)

  expect_identical(class(m), "sdmTMB")
})

test_that("A model with te() spline works", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")

  # test regular model w/o spline works
  d <- subset(pcod, year >= 2015)
  pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 25) # a coarse mesh for example speed
  m <- sdmTMB(data = d,
              formula = density ~ te(depth_scaled),
              spde = pcod_spde,
              family = tweedie(link = "log"))
  expect_identical(class(m), "sdmTMB")

})


test_that("A model with by in spline works", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")

  # test regular model w/o spline works
  d <- subset(pcod, year >= 2015)
  pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 25) # a coarse mesh for example speed
  m <- sdmTMB(data = d,
              formula = density ~ s(depth_scaled,by=year),
              spde = pcod_spde,
              family = tweedie(link = "log"))
  expect_identical(class(m), "sdmTMB")

})

test_that("A model with splines, and no space terms work", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")

  # test regular model w/o spline works
  d <- subset(pcod, year >= 2015)
  m <- sdmTMB(data = d,
              formula = density ~ s(depth_scaled,k=4),
              spde = pcod_spde,
              family = tweedie(link = "log"))
  expect_identical(class(m), "sdmTMB")

})
