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

  # expect_equal(round(p[1:2, 1:10], 5),
  #   structure(c(1.71569, 1.83575, -1.33492, -1.27293, 0.38908, 0.70163,
  #     1.45686, 1.60475, 2.30503, 2.1425, 1.34876, 1.33194, 4.38547,
  #     4.14214, 1.29596, 1.0981, 0.50995, 0.37118, 1.85081, 1.80129), .Dim = c(2L,
  #       10L), .Dimnames = list(c("0", "0"), NULL)))

  .mean <- apply(p, 1, mean)
  .sd <- apply(p, 1, sd)
  expect_gt(cor(.mean, p1$est), 0.99)
})

test_that("get_index_sims works", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  m <- sdmTMB(density ~ 0 + as.factor(year),
    data = pcod_2011, spde = pcod_mesh_2011, family = tweedie(link = "log"),
    time = "year", spatiotemporal = "off"
  )
  qcs_grid_2011 <- subset(qcs_grid, year >= 2011)
  set.seed(1029)
  p <- predict(m, newdata = qcs_grid_2011, sims = 200L)
  expect_equal(ncol(p), 200L)
  expect_equal(nrow(p), nrow(qcs_grid_2011))

  # library(dplyr)
  # a <- reshape2::melt(p)
  # a <- group_by(a, Var1, Var2) %>% summarize(est = sum(exp(value)))
  # a <- group_by(a, Var1) %>% summarise(est = median(est))

  x <- get_index_sims(p)
  expect_equal(nrow(x), length(unique(qcs_grid_2011$year)))
  expect_true(sum(is.na(x$se)) == 0L)

  p_regular <- predict(m, newdata = qcs_grid_2011, return_tmb_object = TRUE)
  x_regular <- get_index(p_regular)
  # expect_equal(round(x_regular$est/x$est, 5),
  #   c(0.91494, 0.92115, 0.92244, 0.9049))

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
