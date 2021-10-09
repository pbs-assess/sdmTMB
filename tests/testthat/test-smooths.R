test_that("A model with 2 s() splines works", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  d <- subset(pcod, year >= 2000 & density > 0)
  pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 30)
  m <- sdmTMB(data = d,
    formula = log(density) ~ s(depth_scaled) + s(year, k = 5),
    spde = pcod_spde,
    control = sdmTMBcontrol(map_rf = TRUE))
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
  plot(p$est, p_mgcv);abline(a = 0, b = 1)
  expect_gt(cor(p$est, p_mgcv), 0.999)
})

test_that("A model with t2() spline works", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  d <- subset(pcod, year >= 2000 & density > 0)
  pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 30)
  m <- sdmTMB(data = d,
    formula = log(density) ~ t2(depth_scaled, year),
    spde = pcod_spde, control = sdmTMBcontrol(map_rf = TRUE))
  p <- predict(m, newdata = NULL)
  m_mgcv <- mgcv::gam(log(density) ~ t2(depth_scaled, year), data = d,
    method = "REML")
  p_mgcv <- predict(m_mgcv, newdata = d)
  plot(p$est, p_mgcv);abline(a = 0, b = 1)
  expect_gt(cor(p$est, p_mgcv), 0.99)
})

test_that("A model with by in spline works", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  d <- subset(pcod, year >= 2006 & density > 0)
  pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 30)
  d$year_factor <- as.factor(d$year)
  expect_warning(m <- sdmTMB(data = d,
    formula = log(density) ~ s(depth_scaled, by = year_factor),
    spde = pcod_spde, control = sdmTMBcontrol(map_rf = TRUE)))
  expect_equal(ncol(m$tmb_data$X_ij), 1L)
  expect_equal(length(m$tmb_data$Zs), length(unique(d$year)))
  p <- predict(m, newdata = NULL)
  m_mgcv <- mgcv::gam(log(density) ~ s(depth_scaled, by = year_factor),
    data = d, method = "REML")
  p_mgcv <- predict(m_mgcv, newdata = d)
  plot(p$est, p_mgcv);abline(a = 0, b = 1)
  expect_gt(cor(p$est, p_mgcv), 0.85) # FIXME why slightly off? but matches brms better
  # m$tmb_data$b_smooth_start
  # m_brms <- brms::brm(log(density) ~ s(depth_scaled, by = year_factor),
  #   data = d, iter = 500, chains = 1, control = list(adapt_delta = 0.99))
  # p2 <- predict(m_brms)
  # plot(p_mgcv, p2[,1]);abline(a = 0, b = 1)
  # plot(p$est, p2[,1]);abline(a = 0, b = 1)
  # m_brms$model
  # dd <- brms::standata(m_brms)
  # head(dd$Xs)
  # head(m$tmb_data$Xs)
  # tail(dd$Xs)
  # tail(m$tmb_data$Xs)
  # head(dd$Zs_1_1)
  # head(m$tmb_data$Zs[[1]])
  # tail(dd$Zs_6_1)
  # tail(m$tmb_data$Zs[[6]])
  # e <- rstan::extract(m_brms$fit)
  # apply(e$bs, 2, median)
  # m$sd_report
  # median(e$sds_1_1)
  })
