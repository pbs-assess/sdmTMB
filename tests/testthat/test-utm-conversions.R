test_that("UTM conversion works", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("sf")

  d <- data.frame(lat = c(52.1, 53.4), lon = c(-130.0, -131.4))
  x <- add_utm_columns(d, c("lon", "lat"))
  expect_equal(class(x), "data.frame")
  expect_true("X" %in% names(x))
  expect_true("Y" %in% names(x))

  expect_error(add_utm_columns(d, c("xx", "lat"))) # column missing
  d$X <- 1
  expect_error(add_utm_columns(d, c("lon", "lat"))) # X already there
  d$X <- NULL

  dd <- d
  names(dd) <- c("lon", "lat")
  expect_warning(add_utm_columns(dd, c("lat", "lon"))) # reversed
  expect_error(add_utm_columns(dd, c("lat", "lon", "x"))) # too many

  # CRS:
  expect_identical(get_crs(d, c("lon", "lat")), 32609)

  # multiple utm zones:
  d <- data.frame(lat = c(52.1, 53.4, 53), lon = c(-130.0, -131.4, -120))
  expect_warning(get_crs(d, c("lon", "lat")), regexp = "zones")

  # N and S
  d <- data.frame(lat = c(52.1, 53.4, -53), lon = c(-130.0, -131.4, -130))
  expect_warning(get_crs(d, c("lon", "lat")), regexp = "North")
})
