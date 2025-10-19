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

# Test for issue #450: ggpredict() fails with multiple smoothers + an offset
test_that("ggpredict works with multiple smoothers and offset (issue #450)", {
  skip_on_cran()
  skip_if_not_installed("ggeffects")

  # Test 1: Basic model without smoothers
  fit1 <- sdmTMB(
    present ~ depth_scaled,
    data = pcod_2011,
    spatial = "off",
    family = binomial()
  )
  g1 <- ggeffects::ggpredict(fit1, "depth_scaled[all]")
  expect_s3_class(g1, "data.frame")
  expect_true(nrow(g1) > 0)

  # Test 2: Model with single smoother
  fit2 <- sdmTMB(
    present ~ s(depth_scaled, k = 5),
    data = pcod_2011,
    spatial = "off",
    family = binomial()
  )
  g2 <- ggeffects::ggpredict(fit2, "depth_scaled[all]")
  expect_s3_class(g2, "data.frame")
  expect_true(nrow(g2) > 0)

  # Test 3: Model with multiple smoothers
  fit3 <- sdmTMB(
    present ~ s(year, k = 3) + s(depth_scaled, k = 5),
    data = pcod_2011,
    spatial = "off",
    family = binomial()
  )
  g3 <- ggeffects::ggpredict(fit3, "depth_scaled[all]")
  expect_s3_class(g3, "data.frame")
  expect_true(nrow(g3) > 0)

  # Test 4: Model with offset and single smoother
  fit4 <- sdmTMB(
    catch_weight ~ s(depth, k = 3),
    data = dogfish,
    spatial = "off",
    family = tweedie(),
    offset = dogfish$area_swept
  )
  g4 <- suppressMessages(ggeffects::ggpredict(fit4, "depth[all]"))
  expect_s3_class(g4, "data.frame")
  expect_true(nrow(g4) > 0)

  # Test 5: Model with offset and multiple smoothers (the main bug from #450)
  fit5 <- sdmTMB(
    catch_weight ~ s(depth, k = 3) + s(year, k = 3),
    data = dogfish,
    spatial = "off",
    family = tweedie(),
    offset = dogfish$area_swept
  )
  g5 <- suppressMessages(ggeffects::ggpredict(fit5, "depth[all]"))
  expect_s3_class(g5, "data.frame")
  expect_true(nrow(g5) > 0)
  # Check that year is in the "Adjusted for" section
  expect_true("year" %in% names(attributes(g5)$constant.values))

  # Test 6: Model with tensor smoother
  fit6 <- sdmTMB(
    catch_weight ~ t2(depth, year),
    data = dogfish,
    spatial = "off",
    family = tweedie()
  )
  g6 <- ggeffects::ggpredict(fit6, "depth[all]")
  expect_s3_class(g6, "data.frame")
  expect_true(nrow(g6) > 0)

  # Test 7: Verify terms() method works correctly with smoothers
  expect_equal(attr(terms(fit5), "term.labels"), c("depth", "year"))
  expect_equal(attr(terms(fit6), "term.labels"), c("depth", "year"))
})
