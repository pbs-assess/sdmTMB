test_that("coef and vcov and confint work", {
  skip_on_ci()
  skip_on_cran()
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

test_that("various methods work", {
  skip_on_ci()
  skip_on_cran()
  fit <- sdmTMB(
    density ~ depth,
    data = pcod_2011, spatial = "off",
    family = tweedie(link = "log")
  )
  f <- fitted(fit)
  expect_equal(nrow(pcod_2011), length(f))

  fit <- sdmTMB(
    density ~ depth,
    data = pcod_2011, spatial = "off",
    family = delta_gamma()
  )
  f <- fitted(fit)
  expect_equal(nrow(pcod_2011), length(f))

  a <- AIC(fit)
  expect_equal(round(a, 3), 6062.726)

  f <- fixef(fit)
  expect_length(f, 2L)

  f <- family(fit)
  expect_identical(f$family, c("binomial", "Gamma"))

  x <- terms(fit)
  expect_identical(as.character(x), c("~", "density", "depth"))

  pcod_2011$fyear <- as.factor(pcod_2011$year)
  fit <- sdmTMB(
    density ~ (1 | fyear),
    data = pcod_2011, spatial = "off",
    family = tweedie(link = "log")
  )
  x <- ranef(fit)
  expect_length(x$cond$fyear[,1], 4L)
})
