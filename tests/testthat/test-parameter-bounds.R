test_that("check_bounds works", {
  check_bounds(c(1, 2, 3), lower = c(0, 0, 0), upper = c(5, 5, 5))
  expect_warning(check_bounds(c(0, 2, 3), lower = c(0, 0, 0), upper = c(5, 5, 5)))
  p <- c(1, 2)
  names(p) <- c("a", "b")
  expect_warning(check_bounds(p, c(1, 1.9), c(5, 5)), regexp = "lower")
  expect_warning(check_bounds(p, c(0, 0), c(1, 5)), regexp = "upper")
})

test_that("lower and upper work", {
  skip_on_cran()
  d <- subset(pcod, year == 2011)
  pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 30)
  expect_warning({
    suppressMessages({
      m <- sdmTMB(density ~ depth_scaled,
        data = d, mesh = pcod_spde, family = tweedie(link = "log"),
        control = sdmTMBcontrol(
          newton_loops = 0,
          lower = list(ln_phi = 0),
          upper = list(ln_phi = 2.5)))
    })},
    regexp = "upper")
  expect_equal(m$model$par[["ln_phi"]], 2.5, tolerance = 1e-6)
  # FIXME NEWTON LOOPS GOING OUTSIDE BOUNDS?

  expect_warning({
    suppressMessages({
      m <- sdmTMB(density ~ depth_scaled,
        data = d, mesh = pcod_spde, family = tweedie(link = "log"),
        control = sdmTMBcontrol(
          newton_loops = 0,
          lower = list(b_j = c(3, -0.46)),
          upper = list(b_j = c(3.5, -0.45))
        ))
    })}, regexp = "bound")
  expect_equal(m$model$par[[2]], -0.45, tolerance = 1e-6)
})
