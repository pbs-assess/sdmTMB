test_that("reloading a model fit with another sdmTMB version warns", {
  withr::local_options(rlib_warning_verbosity = "verbose")
  fit <- sdmTMB(density ~ 1, data = pcod_2011, spatial = "off")
  fit$version <- "0.0.1"
  reloaded <- unserialize(serialize(fit, NULL))
  expect_warning(capture.output(print(reloaded)), "fit with sdmTMB 0.0.1")

  fit$version <- utils::packageVersion("sdmTMB")
  reloaded <- unserialize(serialize(fit, NULL))
  expect_no_warning(capture.output(print(reloaded)))
})

test_that("fitted() on a very old model suggests refitting", {
  fit <- sdmTMB(density ~ 1, data = pcod_2011, spatial = "off")
  fit$offset <- NULL
  expect_error(fitted(fit), "refit the model")
})
