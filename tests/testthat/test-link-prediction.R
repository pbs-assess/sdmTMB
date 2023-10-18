test_that("Link/response type works", {

  # see also test-delta-population-predictions.R
  # https://github.com/pbs-assess/sdmTMB/issues/110
  skip_on_cran()
  skip_on_ci()

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
  # expect_false("est" %in% names(p))
  expect_lt(mean(p$est1), 1)
  expect_lt(mean(p$est2), 5)

  p <- predict(fit_delt, type = "link", re_form = NA)
  # expect_false("est" %in% names(p))
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
  expect_message(p <- predict(fit, type = "link", re_form = NULL, se_fit = TRUE), regexp = "slow")
  mean(p$est)

  p <- predict(fit_delt, type = "link", re_form = NA, se_fit = TRUE)
  mean(p$est)

  expect_warning(
    p <- predict(fit_delt, type = "response", re_form = NA, se_fit = TRUE), regexp = "link"
  )
})

test_that("Response prediction works as reported in https://github.com/pbs-assess/sdmTMB/issues/160#issuecomment-1380920333", {
  skip_on_cran()
  skip_on_ci()

  fit_dg <- sdmTMB(density ~ 1 + s(depth),
    data = pcod_2011,
    mesh = pcod_mesh_2011,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = delta_gamma()
  )

  nd <- data.frame(
    depth = seq(min(pcod$depth),
      max(pcod$depth),
      length.out = 20
    ),
    X = 0,
    Y = 0,
    year = 2015L # a chosen year
  )

  p1 <- predict(fit_dg,
    newdata = nd,
    re_form = NULL,
    re_form_iid = NULL,
    se_fit = TRUE,
    type = "link"
  )

  expect_warning(
    p2 <- predict(fit_dg,
      newdata = nd,
      re_form = NULL,
      re_form_iid = NULL,
      se_fit = TRUE,
      type = "response"
    ),
    regexp = "link"
  )

  expect_equal(p1$est1, p2$est1)
  expect_equal(p1$est2, p2$est2)
  expect_equal(p1$est, p2$est)

  # without se_fit = TRUE
  p1 <- predict(fit_dg,
    newdata = nd,
    re_form = NULL,
    re_form_iid = NULL,
    se_fit = FALSE,
    type = "link"
  )

  p2 <- predict(fit_dg,
    newdata = nd,
    re_form = NULL,
    re_form_iid = NULL,
    se_fit = FALSE,
    type = "response"
  )

  expect_equal(plogis(p1$est1), p2$est1)
  expect_equal(exp(p1$est2), p2$est2)
  expect_equal(plogis(p1$est1) * exp(p1$est2), p2$est)
})
