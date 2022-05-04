## test_that("Epsilon models work with RW spatiotemporal fields", {
##   skip_on_cran()
##   skip_on_ci()
##   skip_if_not_installed("INLA")
##
##   pcod_spde <- pcod_mesh_2011
##   pcod_2011$year_centered <- pcod_2011$year - mean(pcod_2011$year)
##
##   # Fit model with RW fields, no trend
##   m1 <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
##     data = pcod_2011,
##     time = "year",
##     mesh = pcod_spde,
##     family = tweedie(link = "log"),
##     spatiotemporal = "RW"
##   )
##
##   # The way to check that the models are giving the right results is to
##   # create a new dummy variable, include that as a predictor for the time
##   # varying model. It won't fully converge (without fixing the parameter as we
##   # do below), because it's not identifiable, but parameter estimates for
##   # everything else comparable.
##   pcod_2011$dummy <- 0
##   m2 <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
##     data = pcod_2011,
##     time = "year",
##     mesh = pcod_spde,
##     family = tweedie(link = "log"),
##     epsilon_predictor = "dummy",
##     control = sdmTMBcontrol(
##       lower = list(b_epsilon = -1), upper = list(b_epsilon = 1),
##       map = list(b_epsilon = factor(NA)), start = list(b_epsilon = 0)
##     ),
##     spatiotemporal = "RW"
##   )
##
##   expect_equal(tidy(m1, "ran_par")$estimate, tidy(m2, "ran_par")$estimate, tolerance = 0.001)
##   expect_equal(logLik(m1)[1], logLik(m2)[1])
## })
##
## test_that("Epsilon models work with AR1 spatiotemporal fields", {
##   skip_on_cran()
##   skip_on_ci()
##   skip_if_not_installed("INLA")
##
##   pcod_spde <- pcod_mesh_2011
##   pcod_2011$year_centered <- pcod_2011$year - mean(pcod_2011$year)
##
##   # Fit model with AR1 fields, no trend
##   m1 <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
##     data = pcod_2011,
##     time = "year",
##     mesh = pcod_spde,
##     family = tweedie(link = "log"),
##     spatiotemporal = "AR1"
##   )
##
##   # The way to check that the models are giving the right results is to
##   # create a new dummy variable, include that as a predictor for the time
##   # varying model. It won't fully converge (without fixing the parameter as we
##   # do below), because it's not identifiable, but parameter estimates for
##   # everything else comparable.
##   pcod_2011$dummy <- 0
##   m2 <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
##     data = pcod_2011,
##     time = "year",
##     mesh = pcod_spde,
##     family = tweedie(link = "log"),
##     epsilon_predictor = "dummy",
##     control = sdmTMBcontrol(
##       lower = list(b_epsilon = -1), upper = list(b_epsilon = 1),
##       map = list(b_epsilon = factor(NA)), start = list(b_epsilon = 0)
##     ),
##     spatiotemporal = "AR1"
##   )
##
##   expect_equal(tidy(m1, "ran_par")$estimate, tidy(m2, "ran_par")$estimate, tolerance = 0.001)
##   expect_equal(logLik(m1)[1], logLik(m2)[1])
## })
