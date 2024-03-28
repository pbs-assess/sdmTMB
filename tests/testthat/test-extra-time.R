test_that("extra time, newdata, and offsets work", {
  # https://github.com/pbs-assess/sdmTMB/issues/270
  skip_on_cran()
  pcod$os <- rep(log(0.01), nrow(pcod)) # offset
  m <- sdmTMB(
    data = pcod,
    formula = density ~ 0,
    time_varying = ~ 1,
    offset = pcod$os,
    family = tweedie(link = "log"),
    spatial = "off",
    time = "year",
    extra_time = c(2006, 2008, 2010, 2012, 2014, 2016),
    spatiotemporal = "off"
  )
  p1 <- predict(m, offset = pcod$os)
  p2 <- predict(m, newdata = pcod, offset = pcod$os)
  p3 <- predict(m, newdata = pcod)
  p4 <- predict(m, newdata = pcod, offset = rep(0, nrow(pcod)))
  expect_equal(nrow(p1), nrow(pcod))
  expect_equal(nrow(p2), nrow(pcod))
  expect_equal(nrow(p3), nrow(pcod))
  expect_equal(nrow(p4), nrow(pcod))
  expect_equal(p1$est, p2$est)
  expect_equal(p3$est, p4$est)

  #273 (with nsim)
  set.seed(1)
  suppressWarnings(p5 <- predict(m, offset = pcod$os, nsim = 2L))
  expect_equal(ncol(p5), 2L)
  expect_equal(nrow(p5), nrow(pcod))

  set.seed(1)
  suppressWarnings(p6 <- predict(m, newdata = pcod, offset = pcod$os, nsim = 2L))
  expect_equal(ncol(p6), 2L)
  expect_equal(nrow(p6), nrow(pcod))
  expect_equal(p6[, 1, drop = TRUE], p5[, 1, drop = TRUE])

  f <- fitted(m)
  expect_equal(length(f), 2143L)
  expect_equal(round(unique(f), 2), c(31.13, 61.93, 64.98, 18.73, 22.76, 42.97, 40.66, 51.65, 26.05))
})

test_that("extra_time, newdata, get_index() work", {
  skip_on_cran()
  m <- sdmTMB(
    density ~ 1,
    time_varying = ~ 1,
    time_varying_type = "ar1",
    data = pcod,
    family = tweedie(link = "log"),
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    extra_time = c(2006, 2008, 2010, 2012, 2014, 2016, 2018) # last real year is 2017
  )

  # missing one extra_time
  nd <- replicate_df(pcod, "year", sort(union(unique(pcod$year), m$extra_time)))
  nd <- subset(nd, year != 2018)
  p <- predict(m, newdata = nd, return_tmb_object = TRUE)
  ind <- get_index(p)
  ind

  # all:
  nd <- replicate_df(pcod, "year", sort(union(unique(pcod$year), m$extra_time)))
  p <- predict(m, newdata = nd, return_tmb_object = TRUE)
  ind2 <- get_index(p)
  ind2
  expect_identical(ind2$year, c(
    2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010,
    2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
  ))

  expect_equal(ind[ind$year %in% pcod$year, "est"], ind2[ind2$year %in% pcod$year, "est"])

  # just original:
  nd <- replicate_df(pcod, "year", unique(pcod$year))
  p <- predict(m, newdata = nd, return_tmb_object = TRUE)
  ind3 <- get_index(p)
  ind3

  expect_equal(ind2[ind2$year %in% pcod$year, "est"], ind3[ind3$year %in% pcod$year, "est"])
  expect_identical(as.numeric(sort(unique(ind3$year))), as.numeric(sort(unique(pcod$year))))

  p$fake_nd <- NULL # mimic old sdmTMB
  expect_error(ind4 <- get_index(p))

  # missing some original time:
  nd <- replicate_df(pcod, "year", unique(pcod$year))
  nd <- subset(nd, year != 2017)
  p <- predict(m, newdata = nd, return_tmb_object = TRUE)
  ind5 <- get_index(p)
  expect_equal(ind2[ind2$year %in% nd$year, "est"], ind5[ind5$year %in% nd$year, "est"])

  # with do_index = TRUE
  nd <- replicate_df(pcod, "year", unique(pcod$year))
  m2 <- sdmTMB(
    density ~ 1,
    time_varying = ~ 1,
    time_varying_type = "ar1",
    data = pcod,
    family = tweedie(link = "log"),
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    do_index = TRUE,
    predict_args = list(newdata = nd),
    index_args = list(area = 1), # used to cause crash b/c extra_time
    extra_time = c(2006, 2008, 2010, 2012, 2014, 2016, 2018) # last real year is 2017
  )
  ind6 <- get_index(m2)
  expect_identical(ind6$year, c(2003, 2004, 2005, 2007, 2009, 2011, 2013, 2015, 2017))
  expect_equal(ind3$est, ind6$est, tolerance = 0.1)
})
test_that("extra time does not affect estimation (without time series estimation)", {
  # there was a bad bug at one point where the likelihood (via the weights)
  # wasn't getting turned off for the extra time data!
  skip_on_cran()
  # adding extra time at beginning or end
  m <- sdmTMB(present ~ depth_scaled,
    family = binomial(),
    data = pcod_2011, spatial = "on", mesh = pcod_mesh_2011
  )
  m1 <- sdmTMB(
    present ~ depth_scaled,
    family = binomial(), data = pcod_2011, spatial = "on",
    mesh = pcod_mesh_2011,
    extra_time = 1990
  )
  m2 <- sdmTMB(
    present ~ depth_scaled,
    family = binomial(), data = pcod_2011, spatial = "on",
    mesh = pcod_mesh_2011,
    extra_time = 3000
  )
  expect_equal(m$model$par, m1$model$par)
  expect_equal(m$model$par, m2$model$par)

  # with weights
  set.seed(1)
  w <- rlnorm(nrow(pcod_2011), meanlog = log(1), sdlog = 0.1)
  m <- sdmTMB(present ~ depth_scaled,
    family = binomial(), weights = w,
    data = pcod_2011, spatial = "on", mesh = pcod_mesh_2011
  )
  m1 <- sdmTMB(
    present ~ depth_scaled, weights = w,
    family = binomial(), data = pcod_2011, spatial = "on",
    mesh = pcod_mesh_2011,
    extra_time = 1990
  )
  expect_equal(m$model$par, m1$model$par)

  # with offset as numeric
  o <- log(w)
  m <- sdmTMB(density ~ depth_scaled,
    family = tweedie(), offset = o,
    data = pcod_2011, mesh = pcod_mesh_2011
  )
  m1 <- sdmTMB(density ~ depth_scaled,
    family = tweedie(), offset = o,
    data = pcod_2011, mesh = pcod_mesh_2011,
    extra_time = 1990
  )
  expect_equal(m$model$par, m1$model$par)

  # with offset as character
  pcod_2011$off <- o
  m <- sdmTMB(density ~ depth_scaled,
    family = tweedie(), offset = "off",
    data = pcod_2011, mesh = pcod_mesh_2011
  )
  m1 <- sdmTMB(density ~ depth_scaled,
    family = tweedie(), offset = "off",
    data = pcod_2011, mesh = pcod_mesh_2011,
    extra_time = 1990
  )
  expect_equal(m$model$par, m1$model$par)
})

test_that("factor year extra time clash is detected and warned about", {
  skip_on_cran()
  mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
  expect_warning({fit <- sdmTMB(
    density ~ 0 + as.factor(year),
    time = "year", extra_time = 2030, do_fit = FALSE,
    data = pcod_2011, mesh = mesh,
    family = tweedie(link = "log")
  )}, regexp = "rename")
  pcod_2011$year_f <- factor(pcod_2011$year)
  expect_warning({fit <- sdmTMB(
    density ~ 0 + year_f,
    time = "year_f", do_fit = FALSE, extra_time = factor(2030),
    data = pcod_2011, mesh = mesh,
    family = tweedie(link = "log")
  )})
  pcod_2011$year_f2 <- pcod_2011$year_f
  fit <- sdmTMB(
    density ~ 0 + year_f,
    time = "year_f2", do_fit = FALSE, extra_time = factor(2030),
    data = pcod_2011, mesh = mesh,
    family = tweedie(link = "log")
  )
})

test_that("update() works on an extra_time model", {
  skip_on_cran()
  pcod$os <- rep(log(0.01), nrow(pcod)) # offset check
  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 30)
  m <- sdmTMB(
    data = pcod,
    formula = density ~ 0,
    time_varying = ~ 1,
    mesh = mesh,
    offset = pcod$os,
    family = tweedie(link = "log"),
    spatial = "on",
    time = "year",
    extra_time = c(2012),
    spatiotemporal = "off"
  )
  m2 <- update(m, time_varying_type = "ar1")
  expect_s3_class(m2, "sdmTMB")

  m2 <- update(m, time_varying_type = "ar1", mesh = m$spde)
  expect_s3_class(m2, "sdmTMB")

  m2 <- update(m, time_varying_type = "ar1", extra_time = m$extra_time)
  expect_s3_class(m2, "sdmTMB")

  m2 <- update(m, time_varying_type = "ar1", extra_time = m$extra_time, mesh = m$spde)
  expect_s3_class(m2, "sdmTMB")
})
