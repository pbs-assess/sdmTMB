test_that("rmvnorm sim prediction works", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 10)
  m <- sdmTMB(data = pcod,
    formula = density ~ 0 + as.factor(year),
    spde = mesh, family = tweedie(link = "log"))
  set.seed(1)
  p <- predict(m, newdata = qcs_grid, sim = 100L)
  p1 <- predict(m, newdata = qcs_grid)
  expect_identical(class(p)[[1]], "matrix")
  expect_identical(ncol(p), 100L)

  expect_equal(p[1:2, 1:10],
    structure(c(1.71568612013519, 1.83575278171995, -1.33492461718702,
      -1.27293303646145, 0.389082530561698, 0.701626838525375, 1.4568565871509,
      1.60475341525115, 2.30503024445254, 2.14250496100078, 1.3487550058187,
      1.33193607741541, 4.38547481705569, 4.14213965225322, 1.29596339134799,
      1.09810432208885, 0.509954819233504, 0.371179632960271, 1.85081423956647,
      1.80129275148604), .Dim = c(2L, 10L), .Dimnames = list(c("0",
        "0"), NULL)),
    tolerance = 1e-6)

  .mean <- apply(p, 1, mean)
  .sd <- apply(p, 1, sd)
  expect_gt(cor(.mean, p1$est), 0.999)
  # qcs_grid$se <- .sd
  # library(ggplot2)
  # ggplot(qcs_grid, aes(X, Y, fill = se)) + geom_raster()
})
