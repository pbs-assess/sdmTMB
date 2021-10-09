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

#
# test_that("A model with t2() spline works", {
#   skip_on_cran()
#   skip_on_ci()
#   skip_if_not_installed("INLA")
#   d <- subset(pcod, year >= 2000 & density > 0)
#   pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 30)
#   m <- sdmTMB(data = d,
#     formula = log(density) ~ t2(depth_scaled, year),
#     spde = pcod_spde, control = sdmTMBcontrol(map_rf = TRUE))
#   expect_equal(ncol(m$tmb_data$X_ij), 1L)
#   expect_equal(length(m$tmb_data$Zs), 1L)
#   p <- predict(m, newdata = NULL)
#   m_mgcv <- mgcv::gam(log(density) ~ t2(depth_scaled, year), data = d)
#   p_mgcv <- predict(m_mgcv, newdata = d)
#   plot(p$est, p_mgcv);abline(a = 0, b = 1)
#   expect_gt(cor(p$est, p_mgcv), 0.99)
# })

# test_that("A model with by in spline works", {
#   skip_on_cran()
#   skip_on_ci()
#   skip_if_not_installed("INLA")
#
#   # test regular model w/o spline works
#   d <- subset(pcod, year >= 2015)
#   pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 25) # a coarse mesh for example speed
#   m <- sdmTMB(data = d,
#               formula = density ~ s(depth_scaled,by=year),
#               spde = pcod_spde,
#               family = tweedie(link = "log"))
#   expect_identical(class(m), "sdmTMB")
#
# })
