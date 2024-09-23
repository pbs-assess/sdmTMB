test_that("make_mesh returns informative error if some coordinates are NA", {
  d <- data.frame(x = c(NA, runif(10)), y = c(NA, runif(10)))
  expect_error({make_mesh(d, xy_cols = c("x", "y"), cutoff = 1)}, regexp = "NA")

  d <- data.frame(x = c(1, runif(10)), y = c(NA, runif(10)))
  expect_error({make_mesh(d, xy_cols = c("x", "y"), cutoff = 1)}, regexp = "NA")

  d <- data.frame(x = c(NA, runif(10)), y = c(1, runif(10)))
  expect_error({make_mesh(d, xy_cols = c("x", "y"), cutoff = 1)}, regexp = "NA")

  d <- data.frame(x = c(1, runif(10)), y = c(1, runif(10)))
  mesh <- make_mesh(d, xy_cols = c("x", "y"), cutoff = 1)
})
