test_that("Prediction works with missing time", {
  skip_on_cran()
  fit <- sdmTMB(
    density ~ 1,
    data = pcod_2011, mesh = pcod_mesh_2011, time = "year",
    family = tweedie(link = "log")
  )
  nd <- pcod_2011[pcod_2011$year %in% c(2013, 2017), ]
  p1 <- predict(fit, newdata = nd)

  p2 <- predict(fit, newdata = pcod_2011)
  p2 <- p2[p2$year %in% c(2013, 2017), ]

  expect_equal(nrow(p1), nrow(p2))
  expect_equal(p1$est, p2$est)
  expect_equal(p1$year, p2$year)
  expect_equal(p1, p2)

  expect_warning(p3 <- predict(fit, newdata = nd, return_tmb_object = TRUE), regexp = "time")
  expect_error(get_index(p3), regexp = "time")
})
