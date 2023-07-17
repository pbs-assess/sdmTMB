test_that("coef and vcov and confint work", {
  skip_if_not_installed("INLA")
  fit <- sdmTMB(
    density ~ depth,
    data = pcod_2011, spatial = "off",
    family = tweedie(link = "log")
  )
  x <- coef(fit)
  expect_equal(round(unname(x), 3), c(5.347, -0.010))
  expect_equal(names(x), c("(Intercept)", "depth"))

  x <- vcov(fit)
  expect_equal(nrow(x), 2L)
  expect_equal(ncol(x), 2L)
  expect_equal(colnames(x)[1], "(Intercept)")
  expect_equal(rownames(x)[1], "(Intercept)")

  x <- vcov(fit, complete = TRUE)
  expect_equal(nrow(x), 4L)
  expect_equal(ncol(x), 4L)

  x <- confint(fit)
  expect_equal(nrow(x), 2L)
  expect_equal(ncol(x), 3L)
  expect_true(grepl("2\\.5", colnames(x))[1])
  expect_true(grepl("97\\.5", colnames(x))[2])
  expect_true(grepl("Estimate", colnames(x))[3])
})

