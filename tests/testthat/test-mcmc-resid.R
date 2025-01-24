# test_that("Test that MCMC residuals are working", {
#   skip_on_cran()
#   skip_if_not_installed("INLA")
#   skip_if_not_installed("tmbstan")
#
#   set.seed(1)
#   x <- stats::runif(500, -1, 1)
#   y <- stats::runif(500, -1, 1)
#   loc <- data.frame(x = x, y = y)
#   spde <- make_mesh(loc, c("x", "y"), n_knots = 50, type = "kmeans")
#
#   s <- sdmTMB_simulate(
#     ~1,
#     data = loc,
#     mesh = spde,
#     range = 1.4,
#     phi = 0.1,
#     sigma_O = 0.2,
#     seed = 1,
#     B = 0
#   )
#   m1 <- sdmTMB(
#     data = s, time = NULL,
#     formula = observed ~ 1,
#     mesh = spde
#   )
#   resid_mle <- residuals(m1)
#   resid_mcmc <- residuals(m1, type = "mle-mcmc", mcmc_iter = 1000)
#   expect_lt(abs(mean(resid_mcmc)), 0.2)
#   expect_equal(sd(resid_mcmc), 1, tolerance = 0.02)
#   expect_equal(cor(resid_mcmc, resid_mle), 0.955, tolerance = 1e-2)
#
#   # binomial example from scratch/stan-testing
#   set.seed(1)
#   pcod_spde <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
#   m2 <- sdmTMB(present ~ 0 + as.factor(year),
#     data = pcod_2011, mesh = pcod_spde,
#     family = binomial(link = "logit"),
#     priors = sdmTMBpriors(
#       matern_s = pc_matern(range_gt = 15, sigma_lt = 5),
#       b = normal(rep(0, 4), rep(20, 4))
#     )
#   )
#   resid_mle <- residuals(m2)
#   stats::qqnorm(resid_mle)
#   stats::qqline(resid_mle)
#   resid_mcmc <- residuals(m2, type = "mle-mcmc", mcmc_iter = 1000)
#   stats::qqnorm(resid_mcmc)
#   stats::qqline(resid_mcmc)
#
#   expect_lt(abs(mean(resid_mcmc)), 0.1)
#   expect_equal(sd(resid_mcmc), 1.013, tolerance = 1e-2)
#   expect_equal(cor(resid_mcmc, resid_mle), 0.523, tolerance = 1e-2)
# })
