# basic model fitting and prediction tests

SEED <- 1
set.seed(SEED)

test_that("A model with splines works", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")

  d <- subset(pcod, year >= 2015)
  pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 25) # a coarse mesh for example speed
  m <- sdmTMB(data = d,
    formula = density ~ s(depth_scaled),
    spde = pcod_spde,
    family = tweedie(link = "log"))
  expect_identical(class(m), "sdmTMB")
})
