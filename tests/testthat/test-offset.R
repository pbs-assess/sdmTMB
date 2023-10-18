test_that("Offset works", {
  skip_on_cran()
  skip_on_ci()

  pcod$offset <- rnorm(nrow(pcod))
  fit2 <- sdmTMB(density ~ 1,
    offset = pcod$offset,
    data = pcod, spatial = "off",
    family = tweedie()
  )
  expect_error(fit2 <- sdmTMB(density ~ 1,
    offset = year,
    data = pcod, spatial = "off",
    family = tweedie()
  ), regexp = "year")
})


test_that("Offset matches glmmTMB", {
  skip_on_cran()
  skip_on_ci()

  set.seed(1)
  pcod$offset <- rnorm(nrow(pcod))
  pcod_pos <- subset(pcod, density > 0)

  fit1 <- glmmTMB::glmmTMB(density ~ 1,
    data = pcod_pos,
    family = Gamma(link = "log")
  )

  fit1_off <- glmmTMB::glmmTMB(density ~ 1,
    offset = pcod_pos$offset,
    data = pcod_pos,
    family = Gamma(link = "log")
  )

  pcod_pos$offset2 <- log(1)
  fit1_off0 <- glmmTMB::glmmTMB(density ~ 1,
    offset = pcod_pos$offset2,
    data = pcod_pos,
    family = Gamma(link = "log")
  )

  fit2 <- sdmTMB(density ~ 1,
    data = pcod_pos, spatial = "off",
    family = Gamma(link = "log")
  )

  fit2_off <- sdmTMB(density ~ 1,
    offset = pcod_pos$offset,
    data = pcod_pos, spatial = "off",
    family = Gamma(link = "log")
  )

  fit2_off0 <- sdmTMB(density ~ 1,
    offset = pcod_pos$offset2,
    data = pcod_pos, spatial = "off",
    family = Gamma(link = "log")
  )

  b1 <- summary(fit1)$coefficients$cond[1]
  b1_offset <- summary(fit1_off)$coefficients$cond[1]
  b1_offset0 <- summary(fit1_off0)$coefficients$cond[1]

  b2 <- tidy(fit2)$estimate[1]
  b2_offset <- tidy(fit2_off)$estimate[1]
  b2_offset0 <- tidy(fit2_off0)$estimate[1]

  # test that glmmTMB and sdmTMB agree
  expect_equal(b2, b1, tolerance = 1e-4)
  expect_equal(b2_offset, b1_offset, tolerance = 1e-4)
  expect_equal(b2_offset0, b1_offset0, tolerance = 1e-4)

  # the offset of 0 is same as no offset
  expect_equal(b2_offset0, b2, tolerance = 1e-8)

  # the offset is doing something
  expect_false(((b2_offset - b2) == 0))
})

test_that("Offset works with extra_time", {
  skip_on_cran()
  skip_on_ci()
  set.seed(1)
  pcod$offset <- rnorm(nrow(pcod))
  mesh <- make_mesh(pcod, xy_cols = c("X", "Y"), n_knots = 80)
  fit <- sdmTMB(density ~ 1,
    offset = pcod$offset,
    mesh = mesh,
    time = "year",
    extra_time = c(2006L, 2008L, 2010L, 2012L, 2014L, 2016L),
    data = pcod, spatial = "off",
    spatiotemporal = "ar1",
    family = tweedie()
  )
  expect_true(inherits(fit, "sdmTMB"))
  b <- tidy(fit, "ran_pars")
  expect_equal(round(b$estimate[b$term == "rho"], 2), 0.91)
})

test_that("Offset prediction matches glm()", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  set.seed(1)
  pcod$offset <- rnorm(nrow(pcod))
  fit <- sdmTMB(
    present ~ 1,
    offset = pcod$offset,
    data = pcod, spatial = "off",
    family = binomial()
  )
  fit_glm <- glm(
    present ~ 1,
    offset = pcod$offset,
    data = pcod,
    family = binomial()
  )
  fit_glmmTMB <- glmmTMB::glmmTMB(
    present ~ 1,
    offset = pcod$offset,
    data = pcod,
    family = binomial()
  )

  p <- predict(fit)
  p_glm <- predict(fit_glm)
  p_glmmTMB <- predict(fit_glmmTMB)

  expect_equal(p$est, unname(p_glm))
  expect_equal(p$est, p_glmmTMB)

  set.seed(1)
  p <- predict(fit, nsim = 1000)
  mu <- apply(p, 1, mean)
  plot(mu, p_glm)
  expect_equal(unname(mu), unname(p_glm), tolerance = 0.01)

  # sdmTMB ignores offset here (but not glm() or glmmTMB()!)
  # p <- predict(fit, newdata = pcod)
  # p_glmmTMB <- predict(fit_glmmTMB, newdata = pcod)
  # expect_equal(p$est, unname(p_glmmTMB))
})
