# basic model fitting and prediction tests

test_that("A model without splines works", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")

  # test regular model w/o spline works
  d <- subset(pcod, year >= 2015)
  pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 25) # a coarse mesh for example speed
  m <- sdmTMB(data = d,
    formula = density ~ depth_scaled,
    spde = pcod_spde,
    family = tweedie(link = "log"))
  expect_identical(class(m), "sdmTMB")

})

test_that("A model with s() spline works", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")

  # test regular model w/o spline works
  d <- subset(pcod, year >= 2015)
  pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 25) # a coarse mesh for example speed
  m <- sdmTMB(data = d,
              formula = density ~ s(depth_scaled),
              spde = pcod_spde,
              family = tweedie(link = "log"))
  expect_identical(class(m), "sdmTMB")

})

test_that("A model with te() spline works", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")

  # test regular model w/o spline works
  d <- subset(pcod, year >= 2015)
  pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 25) # a coarse mesh for example speed
  m <- sdmTMB(data = d,
              formula = density ~ te(depth_scaled),
              spde = pcod_spde,
              family = tweedie(link = "log"))
  expect_identical(class(m), "sdmTMB")

})


test_that("A model with by in spline works", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")

  # test regular model w/o spline works
  d <- subset(pcod, year >= 2015)
  pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 25) # a coarse mesh for example speed
  m <- sdmTMB(data = d,
              formula = density ~ s(depth_scaled,by=year),
              spde = pcod_spde,
              family = tweedie(link = "log"))
  expect_identical(class(m), "sdmTMB")

})
