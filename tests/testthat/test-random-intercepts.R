test_that("RE group factor levels are properly checked.", {
  expect_error(check_valid_factor_levels(c(1, 2, 3), "test"))
  expect_error(check_valid_factor_levels(c("A", "B")))
  x <- factor(c("a", "b", "c"))
  expect_true(check_valid_factor_levels(x))
  x <- x[-1]
  expect_error(check_valid_factor_levels(x, "test"))
})

test_that("Model with random intercepts fits appropriately.", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("glmmTMB")
  set.seed(1)
  x <- stats::runif(500, -1, 1)
  y <- stats::runif(500, -1, 1)
  loc <- data.frame(x = x, y = y)
  spde <- make_mesh(loc, c("x", "y"), n_knots = 50, type = "kmeans")

  s <- sdmTMB_simulate(
    ~ 1,
    data = loc,
    mesh = spde,
    range = 1.4,
    phi = 0.1,
    sigma_O = 0.2,
    seed = 1,
    B = 0
  )

  g <- rep(gl(30, 10), 999)
  set.seed(134)
  RE_vals <- rnorm(30, 0, 0.4)
  h <- rep(gl(40, 10), 999)
  set.seed(1283)
  RE_vals2 <- rnorm(40, 0, 0.2)
  s$g <- g[seq_len(nrow(s))]
  s$h <- h[seq_len(nrow(s))]
  s$observed <- s$observed + RE_vals[s$g] + RE_vals2[s$h]

  # ignore RE:
  m1 <- sdmTMB(data = s, formula = observed ~ 1, mesh = spde)
  tidy(m1, "fixed", conf.int = TRUE)
  .t1 <- tidy(m1, "ran_pars", conf.int = TRUE)

  # with RE:
  m <- sdmTMB(data = s, time = NULL,
    formula = observed ~ 1 + (1 | g) + (1 | h), mesh = spde)
  tidy(m, "fixed", conf.int = TRUE)
  .t <- tidy(m, "ran_pars", conf.int = TRUE)
  print(m)

  expect_gt(.t1$estimate[.t1$term == "phi"], .t$estimate[.t$term == "phi"])
  expect_gt(.t1$estimate[.t1$term == "sigma_O"], .t$estimate[.t$term == "sigma_O"])
  expect_lt(nrow(.t1), nrow(.t))

  b <- as.list(m$sd_report, "Estimate")
  .cor <- cor(c(RE_vals, RE_vals2), b$RE[,1])
  expect_equal(round(.cor, 5), 0.8313)
  expect_equal(round(b$RE[seq_len(5)], 5),
    c(-0.28645, 0.68619, 0.10028, -0.31436, -0.61168), tolerance = 1e-5)

  # missing a factor level:
  s_drop <- s[s$g != 1, , drop = FALSE]
  spde_drop <- make_mesh(s_drop, c("x", "y"), n_knots = 10, type = "kmeans")
  expect_error(
    sdmTMB(data = s_drop,
      formula = observed ~ 1 + (1 | g) + (1 | h), mesh = spde_drop),
    regexp = "levels"
  )

  p <- predict(m)
  p.nd <- predict(m, newdata = s)
  # newdata is not the same as fitted data:
  p.nd2 <- predict(m, newdata = s[1:3, , drop = FALSE])

  expect_equal(p.nd2$est[1:3], p$est[1:3], tolerance = 1e-4)
  expect_equal(p.nd2$est_non_rf[1:3], p$est_non_rf[1:3], tolerance = 1e-4)
  expect_equal(p.nd2$est[1:3], p.nd$est[1:3], tolerance = 1e-9)
  expect_equal(p$est, p.nd$est, tolerance = 1e-4)
  expect_equal(p$est_rf, p.nd$est_rf, tolerance = 1e-4)
  expect_equal(p$est_non_rf, p.nd$est_non_rf, tolerance = 1e-4)

  # prediction with missing level in `newdata` works:
  s_drop <- s[s$g != 1, , drop = FALSE]
  p.nd <- predict(m, newdata = s_drop)
  p <- p[s$g != 1, , drop = FALSE]
  expect_equal(p$est, p.nd$est, tolerance = 1e-4)

  # prediction without random intercepts included:
  p.nd.null <- predict(m, newdata = s, re_form_iid = NULL)
  p.nd.na <- predict(m, newdata = s, re_form_iid = NA)
  p.nd.0 <- predict(m, newdata = s, re_form_iid = ~ 0)
  expect_identical(p.nd.na, p.nd.0)
  expect_false(identical(p.nd.null$est, p.nd.0$est))

  # random ints match glmmTMB exactly:
  m <- sdmTMB(data = s,
    formula = observed ~ 1 + (1 | g) + (1 | h), mesh = spde, spatial = "off")
  .t <- tidy(m, "ran_pars")
  m.glmmTMB <- glmmTMB::glmmTMB(data = s,
    formula = observed ~ 1 + (1 | g) + (1 | h))
  .v <- glmmTMB::VarCorr(m.glmmTMB)
  expect_equal(.t$estimate[.t$term == "sigma_G"][1],
    sqrt(as.numeric(.v$cond$g)), tolerance = 1e-5)
  expect_equal(.t$estimate[.t$term == "sigma_G"][2],
    sqrt(as.numeric(.v$cond$h)), tolerance = 1e-5)

  sdmTMB_re <- as.list(m$sd_report, "Estimate")
  glmmTMB_re <- glmmTMB::ranef(m.glmmTMB)$cond
  expect_equal(c(glmmTMB_re$g$`(Intercept)`, glmmTMB_re$h$`(Intercept)`),
    sdmTMB_re$RE[,1], tolerance = 1e-5)
})


test_that("Tidy returns random intercepts appropriately.", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("glmmTMB")
  set.seed(1)
  x <- stats::runif(500, -1, 1)
  y <- stats::runif(500, -1, 1)
  loc <- data.frame(x = x, y = y)
  spde <- make_mesh(loc, c("x", "y"), n_knots = 50, type = "kmeans")

  s <- sdmTMB_simulate(
    ~ 1,
    data = loc,
    mesh = spde,
    range = 1.4,
    phi = 0.1,
    sigma_O = 0.2,
    seed = 1,
    B = 0
  )

  g <- rep(gl(30, 10), 999)
  set.seed(134)
  RE_vals <- rnorm(30, 0, 0.4)
  h <- rep(gl(40, 10), 999)
  set.seed(1283)
  RE_vals2 <- rnorm(40, 0, 0.2)
  s$g <- g[seq_len(nrow(s))]
  s$h <- h[seq_len(nrow(s))]
  s$observed <- s$observed + RE_vals[s$g] + RE_vals2[s$h]

  # with RE:
  m <- sdmTMB(data = s, time = NULL,
              formula = observed ~ 1 + (1 | g) + (1 | h),
              mesh = spde,
              spatial = "off")

  ranpars <- tidy(m, "ran_pars", conf.int = TRUE)
  expect_equal(ranpars$estimate,
               c(0.1934564, 0.4422655, 0.1960184), tolerance = 1e-5)
  expect_equal(ranpars$conf.low,
               c(0.1812356, 0.3367532, 0.1414036), tolerance = 1e-5)
  ranint <- tidy(m, "ran_vals", conf.int = TRUE)
  expect_equal(ranint$estimate[1:5],
               c(-0.2281940, 0.6663989, 0.1411399, -0.3220671, -0.6363942), tolerance = 1e-5)
  expect_equal(ranint$conf.low[1:5],
               c(-0.5029686,0.3915682,-0.1338101,-0.5979386,-0.9114315), tolerance = 1e-5)

  # Test against same model estimated from glmmTMB
  fit_glmmtmb <- glmmTMB::glmmTMB(data = s,
                                formula = observed ~ 1 + (1 | g) + (1 | h))
  expect_equal(ranef(fit_glmmtmb)$cond$g[[1]],
               ranint$estimate[1:30], tolerance = 1e-5)

  # also check that ranef returns the same thing with same names
  expect_equal(names(ranef(fit_glmmtmb)$cond), names(ranef(m)$cond))

  # and check that they return the same values
  expect_equal(ranef(fit_glmmtmb)$cond$g[[1]], ranef(m)$cond$g[[1]], tolerance = 1e-5)
})

test_that("random slopes throw an error", {
  pcod_2011$fyear <- as.factor(pcod_2011$year)
  expect_error({
    fit <- sdmTMB(
      density ~ 1 + (1 + depth | fyear),
      data = pcod_2011, mesh = pcod_mesh_2011,
      family = tweedie(link = "log")
    )
  }, regexp = "slope")
  expect_error({
    fit <- sdmTMB(
      density ~ 1 + (1 + depth | fyear) + (Y | fyear),
      data = pcod_2011, mesh = pcod_mesh_2011,
      family = tweedie(link = "log")
    )
  }, regexp = "slope")
  expect_error({
    fit <- sdmTMB(
      density ~ 1 + (depth | fyear),
      data = pcod_2011, mesh = pcod_mesh_2011,
      family = tweedie(link = "log")
    )
  }, regexp = "slope")
  fit <- sdmTMB( # but random intercepts still work
    density ~ 1 + (1 | fyear),
    data = pcod_2011, mesh = pcod_mesh_2011,
    family = tweedie(link = "log")
  )
  expect_s3_class(fit, "sdmTMB")
})

