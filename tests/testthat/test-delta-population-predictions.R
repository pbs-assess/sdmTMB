test_that("Standard errors on overall predictions from delta models work", {
  skip_on_cran()
  skip_on_ci()

  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 12)
  fit <- sdmTMB(
    density ~ depth_scaled + I(depth_scaled^2),
    data = pcod,
    mesh = mesh,
    family = delta_gamma(),
  )
  fit

  nd <- data.frame(depth_scaled = seq(min(pcod$depth_scaled), max(pcod$depth_scaled), length.out = 50))

  # link
  p <- predict(fit, model = NA, re_form = NA, newdata = nd, type = "link")
  head(p)

  p <- predict(fit, model = 1, re_form = NA, newdata = nd, type = "link")
  head(p)

  p <- predict(fit, model = 2, re_form = NA, newdata = nd, type = "link")
  head(p)


  # response
  p <- predict(fit, model = NA, re_form = NA, newdata = nd, type = "response")
  head(p)

  p <- predict(fit, model = 1, re_form = NA, newdata = nd, type = "response")
  head(p)

  p <- predict(fit, model = 2, re_form = NA, newdata = nd, type = "response")
  head(p)


  # re_form = NULL, se_fit = FALSE
  nd$X <- mean(pcod$X)
  nd$Y <- mean(pcod$Y)
  p <- predict(fit, model = NA, newdata = nd, type = "response")
  head(p)

  p <- predict(fit, model = 1, newdata = nd, type = "response")
  head(p)

  p <- predict(fit, model = 2, newdata = nd, type = "response")
  head(p)


  # link with se_fit = TRUE, re_form = NA
  p <- predict(fit, model = NA, re_form = NA, newdata = nd, type = "link", se_fit = TRUE)
  head(p)
  expect_true("est" %in% names(p))
  expect_true("est_se" %in% names(p))
  plot(p$depth_scaled, p$est)
  lines(p$depth_scaled, p$est - 2 * p$est_se)
  lines(p$depth_scaled, p$est + 2 * p$est_se)
  expect_equal(as.numeric(p$est), as.numeric(log(plogis(p$est1) * exp(p$est2))))

  p <- predict(fit, model = 1, re_form = NA, newdata = nd, type = "link", se_fit = TRUE)
  head(p)
  plot(p$depth_scaled, p$est)
  lines(p$depth_scaled, p$est - 2 * p$est_se)
  lines(p$depth_scaled, p$est + 2 * p$est_se)
  expect_equal(p$est, p$est1)

  p <- predict(fit, model = 2, re_form = NA, newdata = nd, type = "link", se_fit = TRUE)
  head(p)
  plot(p$depth_scaled, p$est)
  lines(p$depth_scaled, p$est - 2 * p$est_se)
  lines(p$depth_scaled, p$est + 2 * p$est_se)
  expect_equal(p$est, p$est2)

  visreg_delta(fit, xvar = "depth_scaled", nn = 10, model = 1)
  visreg_delta(fit, xvar = "depth_scaled", nn = 10, model = 2)

  # *response* with se_fit = TRUE, re_form = NA
  # should go link anyways with warning
  expect_warning({
    p <- predict(fit, model = NA, re_form = NA, newdata = nd, type = "response", se_fit = TRUE)
  }, regexp = "link")

})
