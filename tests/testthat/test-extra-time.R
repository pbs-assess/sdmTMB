test_that("extra time, newdata, and offsets work", {
  # https://github.com/pbs-assess/sdmTMB/issues/270
  skip_on_cran()
  skip_on_ci()
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
