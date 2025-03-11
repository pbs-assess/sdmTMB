test_that("tidy works", {
  fit <- sdmTMB(
    density ~ depth_scaled,
    data = pcod_2011, mesh = pcod_mesh_2011,
    family = tweedie(link = "log"), time = "year",
  )
  x <- tidy(fit, "ran_pars", conf.int = TRUE)
  x
  expect_true(sum(is.na(x$std.error)) == 0L)

  x <- tidy(fit, conf.int = TRUE)
  expect_true(sum(is.na(x$std.error)) == 0L)

  d <- pcod_2011
  d$year_scaled <- as.numeric(scale(d$year))
  fit <- sdmTMB(
    density ~ year_scaled,
    data = d, mesh = pcod_mesh_2011,
    family = tweedie(link = "log"), time = "year",
    spatiotemporal = "off", spatial_varying = ~ 0 + year_scaled
  )

  x <- tidy(fit, "ran_pars", conf.int = TRUE)
  x
  expect_true(sum(is.na(x$std.error)) == 0L)

  fit <- sdmTMB(
    density ~ depth_scaled,
    data = pcod_2011, mesh = pcod_mesh_2011,
    family = delta_gamma()
  )
  x <- tidy(fit, "ran_pars", conf.int = TRUE)
  x
  expect_true(sum(is.na(x$std.error)) == 0L)

  x <- tidy(fit, "ran_pars", conf.int = TRUE, model = 2)
  x
  expect_true(sum(is.na(x$std.error)) == 0L)

  fit <- sdmTMB(
    density ~ depth_scaled,
    data = pcod_2011, mesh = pcod_mesh_2011, share_range = FALSE,
    family = tweedie(link = "log"), time = "year",
  )
  x <- tidy(fit, "ran_pars", conf.int = TRUE)
  x
  expect_true(sum(is.na(x$std.error)) == 0L)

  # test that parsing of time varying random values works
  mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
  fit <- sdmTMB(
    density ~ 0 + as.factor(year),
    time_varying = ~ 0 + depth_scaled + depth_scaled2, #<
    data = pcod_2011, time = "year", mesh = mesh,
    family = tweedie()
  )
  pars <- tidy(fit, "ran_vals")
  expect_equal(pars$estimate, c(
    -0.87, -0.81, -0.75, -1.11,
    -1.92, -0.92, -1.59, -2.20
  ), tolerance = 0.01)
  expect_equal(pars$term, c(
    "depth_scaled:2011", "depth_scaled:2013", "depth_scaled:2015", "depth_scaled:2017",
    "depth_scaled2:2011", "depth_scaled2:2013", "depth_scaled2:2015", "depth_scaled2:2017"
  ))

  # test smooth handling
  fit <- sdmTMB(
    density ~ s(depth),
    data = pcod_2011, mesh = pcod_mesh_2011,
    family = tweedie(link = "log")
  )
  x <- tidy(fit)
  expect_equal(x$term, c("(Intercept)", "sdepth"))

  # Test with sdmTMB_cv
  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 25)
  m_cv <- sdmTMB_cv(
    density ~ 0 + depth_scaled + depth_scaled2,
    data = pcod, mesh = mesh,
    family = tweedie(link = "log"), k_folds = 2
  )
  model1 <- tidy(m_cv$models[[1]])
  allmodels <- tidy(m_cv)
  expect_equal(allmodels[1:2,1:5], model1[,1:5])
  expect_true("cv_split" %in% names(allmodels))
  expect_equal(allmodels$cv_split, sort(rep(1:2,2)))
})

test_that("printing/tidying works with a delta model that has random intercepts + an AR1 time series #426", {
  skip_on_cran()
  d <- pcod
  d$fake <- rep(c("a", "b", "c"), 9999)[1:nrow(d)]
  fit <- sdmTMB(
    density ~ breakpt(depth_scaled) + (1|fake),
    data = d,
    time = "year",
    extra_time = c(2006, 2008, 2010, 2012, 2014, 2016),
    time_varying = ~1,
    time_varying_type = "ar1",
    spatial = "off",
    spatiotemporal = "off",
    family = delta_gamma(type = "poisson-link")
  )
  b <- tidy(fit, effects = "ran_pars")
  b <- tidy(fit, effects = "ran_vals")
  expect_identical(b$term, c("(Intercept)", "(Intercept)", "(Intercept)", "(Intercept)",
    "(Intercept)", "(Intercept)", "(Intercept):2003", "(Intercept):2004",
    "(Intercept):2005", "(Intercept):2006", "(Intercept):2007", "(Intercept):2008",
    "(Intercept):2009", "(Intercept):2010", "(Intercept):2011", "(Intercept):2012",
    "(Intercept):2013", "(Intercept):2014", "(Intercept):2015", "(Intercept):2016",
    "(Intercept):2017", "(Intercept):2003", "(Intercept):2004", "(Intercept):2005",
    "(Intercept):2006", "(Intercept):2007", "(Intercept):2008", "(Intercept):2009",
    "(Intercept):2010", "(Intercept):2011", "(Intercept):2012", "(Intercept):2013",
    "(Intercept):2014", "(Intercept):2015", "(Intercept):2016", "(Intercept):2017"
  ))
})
