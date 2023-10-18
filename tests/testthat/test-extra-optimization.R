test_that("Extra optimization runs and reduces gradients", {
  skip_on_cran()
  skip_on_ci()
  d <- subset(pcod, year >= 2013)
  pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 30)
  m <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
    data = d, time = "year", mesh = pcod_spde, family = tweedie(link = "log"))

  m1 <- run_extra_optimization(m, nlminb_loops = 1, newton_loops = 1)
  expect_lt(max(m1$gradients), max(m$gradients))
})
