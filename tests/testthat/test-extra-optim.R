context("Extra optimization")

d <- subset(pcod, year >= 2013)
pcod_spde <- make_spde(d$X, d$Y, n_knots = 30)
m <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
  data = d, time = "year", spde = pcod_spde, family = tweedie(link = "log"),
  nlminb_loops = 1)

test_that("Extra optimization runs and reduces gradients", {
  m1 <- run_extra_optimization(m, nlminb_loops = 1, newton_steps = 1)
  expect_lt(max(m1$gradients), max(m$gradients))
})
