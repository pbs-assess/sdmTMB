test_that("Link/response type works", {

  # https://github.com/pbs-assess/sdmTMB/issues/110
  skip_on_cran()
  skip_if_not_installed("INLA")

  fit <- sdmTMB(
    density ~ 1,
    family = tweedie(),
    data = pcod_2011,
    mesh = pcod_mesh_2011
  )

  p <- predict(fit, type = "link")
  expect_lt(mean(p$est), 5)

  p <- predict(fit, type = "link", re_form = NA)
  expect_lt(mean(p$est), 5)

  p <- predict(fit, type = "response")
  expect_gt(mean(p$est), 20)

  p <- predict(fit, type = "response", re_form = NA)
  expect_gt(mean(p$est), 20)

  fit_delt <- sdmTMB(
    density ~ 1,
    family = delta_gamma(),
    data = pcod_2011,
    mesh = pcod_mesh_2011
  )

  p <- predict(fit_delt, type = "link")
  expect_false("est" %in% names(p))
  expect_lt(mean(p$est1), 1)
  expect_lt(mean(p$est2), 5)

  p <- predict(fit_delt, type = "link", re_form = NA)
  expect_false("est" %in% names(p))
  expect_lt(mean(p$est1), 1)
  expect_lt(mean(p$est2), 5)

  p <- predict(fit_delt, type = "response")
  expect_gt(mean(p$est), 20)
  expect_gt(mean(p$est1), 0)
  expect_gt(mean(p$est2), 30)

  p <- predict(fit_delt, type = "response", re_form = NA)
  expect_gt(mean(p$est), 20)
  expect_gt(mean(p$est1), 0)
  expect_gt(mean(p$est2), 30)

  # with std. error
  # p <- predict(fit, type = "link", re_form = NA, se_fit = TRUE)
  # mean(p$est)
  #
  # p <- predict(fit_delt, type = "link", re_form = NA, se_fit = TRUE)
  # mean(p$est)
})
