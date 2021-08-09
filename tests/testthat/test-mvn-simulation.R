test_that("rmvnorm sim prediction works", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 10)
  m <- sdmTMB(data = pcod,
    formula = density ~ 0 + as.factor(year),
    spde = mesh, family = tweedie(link = "log"))
  set.seed(1)
  p <- predict(m, newdata = qcs_grid, sim = 15L)
  p1 <- predict(m, newdata = qcs_grid)
  expect_identical(class(p)[[1]], "matrix")
  expect_identical(ncol(p), 15L)

  expect_equal(p[1:2, 1:10],
    structure(c(1.8699162905496, 1.63531268879039, -0.545589429955865,
      -0.166086444510931, 2.74585833405298, 2.28654710196737, 0.313636663657174,
      0.396008072224339, 2.06241095921304, 1.60963231768441, 1.75521432004763,
      2.07418366348474, 0.333046692746184, 0.555437315814931, 1.97532933260919,
      2.13437800215603, 5.12385511998788, 4.89894219449716, 0.184090550087612,
      0.344108837803646), .Dim = c(2L, 10L), .Dimnames = list(c("0",
        "0"), NULL)),
    tolerance = 1e-5)

  .mean <- apply(p, 1, mean)
  .sd <- apply(p, 1, sd)
  expect_gt(cor(.mean, p1$est), 0.99)
})

test_that("get_index_sims works", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  m <- sdmTMB(density ~ 0 + as.factor(year),
    data = pcod_2011, spde = pcod_mesh_2011, family = tweedie(link = "log"),
    time = "year", spatial_only = TRUE
  )
  qcs_grid_2011 <- subset(qcs_grid, year >= 2011)
  p <- predict(m, newdata = qcs_grid_2011, sims = 3L)
  expect_equal(ncol(p), 3L)
  expect_equal(nrow(p), nrow(qcs_grid_2011))

  x <- get_index_sims(p)
  expect_equal(nrow(x), length(unique(qcs_grid_2011$year)))
  expect_true(sum(is.na(x$se)) == 0L)

  x_sims <- get_index_sims(p, return_sims = TRUE)
  expect_equal(nrow(x_sims), nrow(x) * ncol(p))

  # all areas doubled:
  xa2 <- get_index_sims(p, area = rep(2, nrow(p)))
  expect_equal(xa2$est / x$est, rep(2, 4))

  # 1 year different:
  areas <- rep(1, nrow(p))
  areas[qcs_grid_2011$year == 2011] <- 3.14
  xpi <- get_index_sims(p, area = areas)
  expect_equal(xpi$est / x$est, c(3.14, 1, 1, 1))

  # 5 cells different:
  areas[1:5] <- 0.01
  xpi2 <- get_index_sims(p, area = areas)
  expect_lt(xpi2$est[1], xpi$est[1])
  expect_equal(xpi2$est[-1], xpi$est[-1])

  # arg checking:
  expect_error(get_index_sims(3))
  expect_error(get_index_sims(level = 1.1))
  expect_error(get_index_sims(level = -0.1))
  expect_error(get_index_sims(area = c(1, 2)))
  expect_error(get_index_sims(area = rep(NA, nrow(p))))
  expect_error(get_index_sims(est_function = 3))
  expect_error(get_index_sims(agg_function = 3))
})
