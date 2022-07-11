# basic plot tests

SEED <- 123

test_that("PC prior plotting works", {
  local_edition(2)
  p <- plot_pc_matern(range_gt = 0.5,
                      sigma_lt = 0.1,
                      range_prob = 0.05,
                      sigma_prob = 0.05, plot = FALSE)
  expect_true(is.matrix(p))
  expect_true(is.array(p))
  expect_equal(dim(p)[1], 201)
  expect_equal(dim(p)[2], 200)
  x <- seq(0.5 * 0.1, 0.5 * 10, length.out=200)
  expect_lt(max(abs(as.numeric(colnames(p)) - x)), 1e-12)
  y <- seq(0, 0.1 * 2, length.out=201)
  expect_lt(max(abs(as.numeric(rownames(p)) - y)), 1e-12)

  expect_equal(p[12,37], 2.002172, tolerance = 1e-3)
  expect_equal(p[37,75], 0.6592268, tolerance = 1e-3)
  expect_equal(p[75,114], -1.038814, tolerance = 1e-3)
  expect_equal(p[114,165], -2.780357, tolerance = 1e-3)
})
