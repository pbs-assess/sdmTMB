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
})
