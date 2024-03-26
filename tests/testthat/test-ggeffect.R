## test_that("ggeffects + sdmTMB", {
##   skip_on_cran()
##   skip_if_not_installed("INLA")
##   skip_if_not_installed("ggeffects")
##   skip_if_not_installed("ggplot2")
##
##   pcod_2011$fyear <- as.factor(pcod_2011$year)
##   fit <- sdmTMB(present ~ depth_scaled + I(depth_scaled^2) + fyear,
##     data = pcod_2011,
##     mesh = pcod_mesh_2011,
##     family = binomial()
##   )
##   g <- ggeffects::ggeffect(fit, "depth_scaled [-2.5:2.5, by=.1]")
##   expect_s3_class(g, "data.frame")
##   plot(g)
##
##   e <- effects::effect("depth_scaled", fit)
##   e
##   expect_true(inherits(e, "eff"))
##
##   fit2 <- sdmTMB(present ~ depth_scaled + I(depth_scaled^2) + (1 | fyear),
##     data = pcod_2011,
##     mesh = pcod_mesh_2011,
##     family = binomial()
##   )
##   effects::effect("depth_scaled", fit2)
##   g <- ggeffects::ggeffect(fit2, "depth_scaled [-2.5:2.5, by=.1]")
##   plot(g)
##   expect_s3_class(g, "data.frame")
##
##   fit3 <- sdmTMB(present ~ s(year, k = 3) + depth_scaled + I(depth_scaled^2),
##     data = pcod_2011,
##     mesh = pcod_mesh_2011,
##     family = binomial()
##   )
##   effects::effect("depth_scaled", fit3)
##   g <- ggeffects::ggeffect(fit3, "depth_scaled [-2.5:2.5, by=.1]")
##   plot(g)
##   expect_s3_class(g, "data.frame")
##
##   fit4 <- sdmTMB(present ~ s(depth_scaled, k = 4),
##     data = pcod_2011,
##     mesh = pcod_mesh_2011,
##     family = binomial()
##   )
##   expect_error({
##     effects::effect("depth_scaled", fit4)
##   }, regexp = "missing")
## })

### test_that("ggpredict + sdmTMB", {
###   skip_on_cran()
###   skip_if_not_installed("INLA")
###   skip_if_not_installed("ggeffects")
###   skip_if_not_installed("ggplot2")
###   skip_if_not_installed("sdmTMB")
###
###   d <- pcod_2011
###   d$fyear <- as.factor(d$year)
###
###   # basic quadratic
###   fit <- sdmTMB(present ~ depth_scaled + I(depth_scaled^2) + fyear,
###     data = d,
###     spatial = "off",
###     family = binomial()
###   )
###   g <- ggeffects::ggpredict(fit, "depth_scaled [all]")
###   expect_s3_class(g, "data.frame")
###   plot(g)
###
###   # matches glmmTMB?
###   fit_glmmTMB <- glmmTMB::glmmTMB(
###     present ~ depth_scaled + I(depth_scaled^2) + fyear,
###     data = d,
###     family = binomial()
###   )
###   g_glmmTMB <- ggeffects::ggpredict(fit, "depth_scaled [all]")
###   expect_equal(g, g_glmmTMB, tolerance = 1e-3)
###
###   # with random intercept
###   fit2 <- sdmTMB(
###     present ~ depth_scaled + I(depth_scaled^2) + (1 | fyear),
###     data = d,
###     spatial = "off",
###     family = binomial()
###   )
###   g <- ggeffects::ggpredict(fit2, "depth_scaled [all]")
###   plot(g)
###   expect_s3_class(g, "data.frame")
###
###   # matches glmmTMB?
###   fit_glmmTMB2 <- glmmTMB::glmmTMB(
###     present ~ depth_scaled + I(depth_scaled^2) + (1 | fyear),
###     data = d,
###     family = binomial()
###   )
###   g_glmmTMB <- ggeffects::ggpredict(fit_glmmTMB2, "depth_scaled [all]")
###
###   plot(g)
###   plot(g_glmmTMB)
###   expect_equal(g$predicted, g_glmmTMB$predicted, tolerance = 1e-3)
###   expect_equal(g$std.error, g_glmmTMB$std.error, tolerance = 1e-3)
###
###   expect_error(ggeffects::ggpredict(fit2, "depth_scaled [all]", type = "re"), regexp = "supported")
###
###   # with other smoother terms:
###   fit3 <- sdmTMB(
###     present ~ s(year, k = 3) + depth_scaled + I(depth_scaled^2),
###     data = d,
###     spatial = "off",
###     family = binomial()
###   )
###   g <- ggeffects::ggpredict(fit3, "depth_scaled [all]")
###   expect_s3_class(g, "data.frame")
###
###   # on smoother terms themselves:
###   fit4 <- sdmTMB(
###     present ~ s(depth_scaled, k = 5),
###     data = d,
###     spatial = "off",
###     family = binomial()
###   )
###   g <- ggeffects::ggpredict(fit4, "depth_scaled [all]")
###   plot(g)
###   expect_s3_class(g, "data.frame")
###
###   # similar enough to mgcv?
###   fit_mgcv4 <- mgcv::gam(
###     present ~ s(depth_scaled, k = 5),
###     data = d,
###     family = binomial()
###   )
###   g_mgcv <- ggeffects::ggpredict(fit_mgcv4, "depth_scaled [all]")
###   plot(g_mgcv)
###
###   expect_equal(as.numeric(g$predicted), as.numeric(g_mgcv$predicted), tolerance = 1e-2)
###   expect_equal(as.numeric(g$std.error), as.numeric(g_mgcv$std.error), tolerance = 1e-1)
###
###   # with some random fields:
###   fit_sp <- sdmTMB(
###     present ~ depth_scaled,
###     mesh = pcod_mesh_2011,
###     data = d,
###     spatial = "on",
###     family = binomial()
###   )
###   g <- ggeffects::ggpredict(fit_sp, "depth_scaled [all]")
###   expect_s3_class(g, "data.frame")
###   plot(g)
###
###   gm <- ggeffects::ggemmeans(fit_sp, "depth_scaled [all]")
###   gp <- ggeffects::ggpredict(fit_sp, "depth_scaled [all]")
###   ge <- ggeffects::ggeffect(fit_sp, "depth_scaled [all]")
###   expect_equal(gm$predicted, gp$predicted, tolerance = 1e-3)
###   expect_equal(gm$predicted, ge$predicted, tolerance = 1e-3)
### })
###
