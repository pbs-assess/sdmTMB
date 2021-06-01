test_that("rmvnorm sim prediction works", {
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
    structure(c(1.74529288326404, 1.84540275192375, 0.803440491931372,
    1.11365791215555, 2.13224380123888, 2.33087978879283, 2.38945020735315,
    2.5273650857693, -0.386991701267006, -0.240176121269049, -0.311765593772021,
    -0.339317094353514, 1.38186046416679, 1.65844041727968, 3.06651731693835,
    3.01683843796059, 1.7315609398691, 1.59050905137807, 3.10402659778617,
    3.11058531772459), .Dim = c(2L, 10L)), tolerance = 1e-6)
  .mean <- apply(p, 1, mean)
  .sd <- apply(p, 1, sd)
  expect_gt(cor(.mean, p1$est), 0.999)
  # qcs_grid$se <- .sd
  # library(ggplot2)
  # ggplot(qcs_grid, aes(X, Y, fill = se)) + geom_raster()
})
