test_that("get_kappa_map() works", {
  x <- get_kappa_map(
    n_m = 2,
    spatial = c("off", "off"),
    spatiotemporal = c("off", "off"),
    share_range = c(FALSE, FALSE)
  )
  expect_identical(x, factor(c(NA, NA, NA, NA)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("off", "off"),
    spatiotemporal = c("off", "off"),
    share_range = c(TRUE, TRUE)
  )
  expect_identical(x, factor(c(NA, NA, NA, NA)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("on", "on"),
    spatiotemporal = c("off", "off"),
    share_range = c(TRUE, TRUE)
  )
  expect_identical(x, factor(c(1, 1, 2, 2)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("on", "on"),
    spatiotemporal = c("off", "off"),
    share_range = c(FALSE, FALSE)
  )
  expect_identical(x, factor(c(1, 1, 2, 2)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("off", "off"),
    spatiotemporal = c("on", "on"),
    share_range = c(TRUE, TRUE)
  )
  expect_identical(x, factor(c(1, 1, 2, 2)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("off", "off"),
    spatiotemporal = c("on", "on"),
    share_range = c(FALSE, FALSE)
  )
  expect_identical(x, factor(c(1, 1, 2, 2)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("off", "on"),
    spatiotemporal = c("on", "on"),
    share_range = c(TRUE, TRUE)
  )
  expect_identical(x, factor(c(1, 1, 2, 2)))

  # x <- get_kappa_map(
  #   n_m = 2,
  #   spatial = c("off", "on"),
  #   spatiotemporal = c("on", "on"),
  #   share_range = c(FALSE, FALSE)
  # )
  # expect_identical(x, factor(c(1, 1, 2, 3)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("on", "on"),
    spatiotemporal = c("off", "on"),
    share_range = c(TRUE, TRUE)
  )
  expect_identical(x, factor(c(1, 1, 2, 2)))

  # x <- get_kappa_map(
  #   n_m = 2,
  #   spatial = c("on", "on"),
  #   spatiotemporal = c("off", "on"),
  #   share_range = c(FALSE, FALSE)
  # )
  # expect_identical(x, factor(c(1, 1, 2, 3)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("on", "on"),
    spatiotemporal = c("on", "off"),
    share_range = c(TRUE, TRUE)
  )
  expect_identical(x, factor(c(1, 1, 2, 2)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("on", "on"),
    spatiotemporal = c("on", "off"),
    share_range = c(FALSE, FALSE)
  )
  expect_identical(x, factor(c(1, 2, 3, 3)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("on", "off"),
    spatiotemporal = c("on", "on"),
    share_range = c(TRUE, TRUE)
  )
  expect_identical(x, factor(c(1, 1, 2, 2)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("on", "off"),
    spatiotemporal = c("on", "on"),
    share_range = c(FALSE, FALSE)
  )
  expect_identical(x, factor(c(1, 2, 3, 3)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("on", "on"),
    spatiotemporal = c("on", "on"),
    share_range = c(FALSE, FALSE)
  )
  expect_identical(x, factor(c(1, 2, 3, 4)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("on", "on"),
    spatiotemporal = c("on", "on"),
    share_range = c(TRUE, TRUE)
  )
  expect_identical(x, factor(c(1, 1, 2, 2)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("on", "off"),
    spatiotemporal = c("on", "on"),
    share_range = c(FALSE, TRUE)
  )
  expect_identical(x, factor(c(1, 2, 3, 3)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("on", "off"),
    spatiotemporal = c("off", "off"),
    share_range = c(TRUE, TRUE)
  )
  expect_identical(x, factor(c(1, 1, NA, NA)))

  # x <- get_kappa_map(
  #   n_m = 2,
  #   spatial = c("on", "on"),
  #   spatiotemporal = c("on", "on"),
  #   share_range = c(FALSE, TRUE)
  # )
  # expect_identical(x, factor(c(1, 2, 3, 3)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("on", "on"),
    spatiotemporal = c("on", "on"),
    share_range = c(TRUE, FALSE)
  )
  expect_identical(x, factor(c(1, 1, 2, 3)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("on", "off"),
    spatiotemporal = c("off", "on"),
    share_range = c(FALSE, FALSE)
  )
  expect_identical(x, factor(c(1, 1, 2, 2)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("off", "on"),
    spatiotemporal = c("off", "on"),
    share_range = c(FALSE, FALSE)
  )
  expect_identical(x, factor(c(NA, NA, 1, 2)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("on", "off"),
    spatiotemporal = c("on", "off"),
    share_range = c(FALSE, FALSE)
  )
  expect_identical(x, factor(c(1, 2, NA, NA)))

  # non-delta:
  x <- get_kappa_map(
    n_m = 1,
    spatial = "on",
    spatiotemporal = "on",
    share_range = TRUE
  )
  expect_identical(x, factor(c(1, 1)))

  x <- get_kappa_map(
    n_m = 1,
    spatial = "off",
    spatiotemporal = "on",
    share_range = TRUE
  )
  expect_identical(x, factor(c(1, 1)))

  x <- get_kappa_map(
    n_m = 1,
    spatial = "on",
    spatiotemporal = "off",
    share_range = TRUE
  )
  expect_identical(x, factor(c(1, 1)))

  x <- get_kappa_map(
    n_m = 1,
    spatial = "off",
    spatiotemporal = "off",
    share_range = TRUE
  )
  expect_identical(x, factor(c(NA, NA)))

  #######

  x <- get_kappa_map(
    n_m = 1,
    spatial = "on",
    spatiotemporal = "on",
    share_range = FALSE
  )
  expect_identical(x, factor(c(1, 2)))

  x <- get_kappa_map(
    n_m = 1,
    spatial = "off",
    spatiotemporal = "on",
    share_range = FALSE
  )
  expect_identical(x, factor(c(1, 1)))

  x <- get_kappa_map(
    n_m = 1,
    spatial = "on",
    spatiotemporal = "off",
    share_range = FALSE
  )
  expect_identical(x, factor(c(1, 1)))

  x <- get_kappa_map(
    n_m = 1,
    spatial = "off",
    spatiotemporal = "off",
    share_range = FALSE
  )
  expect_identical(x, factor(c(NA, NA)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("on", "on"),
    spatiotemporal = c("off", "on"),
    share_range = c(FALSE, FALSE)
  )
  expect_identical(x, factor(c(1, 1, 2, 3)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("off", "off"),
    spatiotemporal = c("on", "on"),
    share_range = c(FALSE, TRUE)
  )
  expect_identical(x, factor(c(1, 1, 2, 2)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("off", "off"),
    spatiotemporal = c("on", "on"),
    share_range = c(FALSE, FALSE)
  )
  expect_identical(x, factor(c(1, 1, 2, 2)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("on", "on"),
    spatiotemporal = c("on", "on"),
    share_range = c(FALSE, TRUE)
  )
  expect_identical(x, factor(c(1, 2, 3, 3)))
})
