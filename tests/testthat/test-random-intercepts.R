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
  s <- sdmTMB_sim(x = x, y = y,
    betas = 0, time = 1L,
    phi = 0.1, range = 1.4,
    sigma_O = 0.2, sigma_E = 0,
    seed = 1, mesh = spde
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
  expect_equal(round(.cor, 6), 0.827327)
  expect_equal(round(b$RE[seq_len(5)], 6),
    c(-0.328542, 0.680869, 0.106166, -0.369494, -0.63912), tolerance = 1e-5)

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
    formula = observed ~ 1 + (1 | g) + (1 | h), mesh = spde, control = sdmTMBcontrol(map_rf = TRUE))
  .t <- tidy(m, "ran_pars")
  m.glmmTMB <- glmmTMB::glmmTMB(data = s, formula = observed ~ 1 + (1 | g) + (1 | h))
  .v <- glmmTMB::VarCorr(m.glmmTMB)
  expect_equal(.t$estimate[.t$term == "sigma_G"][1], sqrt(as.numeric(.v$cond$g)), tolerance = 1e-5)
  expect_equal(.t$estimate[.t$term == "sigma_G"][2], sqrt(as.numeric(.v$cond$h)), tolerance = 1e-5)

  sdmTMB_re <- as.list(m$sd_report, "Estimate")
  glmmTMB_re <- glmmTMB::ranef(m.glmmTMB)$cond
  expect_equal(c(glmmTMB_re$g$`(Intercept)`, glmmTMB_re$h$`(Intercept)`),
    sdmTMB_re$RE[,1], tolerance = 1e-7)
})
