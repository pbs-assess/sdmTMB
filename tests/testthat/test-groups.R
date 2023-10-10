test_that("groups work", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  # skip_if_not_installed("ggplot2")
  # library(ggplot2)

  expect_error(make_groups(c(1, 2)), regexp = "not a factor")
  expect_identical(make_groups(factor(c('a', 'b', 'c'))), c(0L, 1L, 2L))
  expect_error(make_groups(factor(c('a', 'b'), levels = c('a', 'b', 'c'))), "Extra factor")
  expect_identical(make_groups(factor(c('a', 'b', 'c')), prev_levels = c('a', 'b', 'c')), c(0L, 1L, 2L))
  expect_identical(make_groups(factor(c('a', 'b'), levels = c('a', 'b', 'c')), prev_levels = c('a', 'b', 'c')), c(0L, 1L))
  expect_identical(make_groups(factor(c('a', 'c'), levels = c('a', 'b', 'c')), prev_levels = c('a', 'b', 'c')), c(0L, 2L))
  expect_error(make_groups(factor(c('a', 'c'), levels = c('a', 'b', 'c', 'd')), prev_levels = c('a', 'b', 'c')), regexp = "Extra")

  s <- lapply(1:3, function(i) {
    predictor_dat <- data.frame(
      X = runif(300), Y = runif(300),
      year = rep(1:6, each = 50)
    )
    mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)
    sdmTMB_simulate(
      formula = ~1,
      data = predictor_dat,
      time = "year",
      mesh = mesh,
      family = gaussian(),
      range = 0.5,
      sigma_E = 0.2,
      phi = 0.1,
      sigma_O = 0,
      spatiotemporal = "rw",
      seed = i,
      B = 0
    )
  })
  s <- do.call("rbind", s)
  s$category <- rep(1:3, each = 300L)
  s$category <- as.factor(s$category)
  head(s)

  # ggplot(s, aes(X, Y, colour = epsilon_st)) +
  #   geom_point() +
  #   facet_grid(category ~ year) +
  #   scale_colour_viridis_c()

  mesh <- make_mesh(s, c("X", "Y"), cutoff = 0.1)
  fit <- sdmTMB(
    observed ~ 0 + as.factor(category),
    data = s,
    time = "year",
    group = "category",
    spatial = "off",
    spatiotemporal = "rw",
    silent = FALSE,
    mesh = mesh
  )

  # fit$tmb_obj$par
  # fit$tmb_params$upsilon_stc
  # dim(fit$tmb_params$upsilon_stc)
  #
  # fit$tmb_map$ln_kappa
  # fit$tmb_map$epsilon_st
  # fit$tmb_map$epsilon_re
  # fit$tmb_map$upsilon_stc
  # fit$tmb_map$ln_tau_E
  # fit$tmb_map$ln_tau_O
  # fit$tmb_data$n_c
  # fit$tmb_data$mvrw_cat_i
  # fit$tmb_data$year_i
  # fit$tmb_data$rw_fields
  # fit$tmb_data$ar1_fields

  nd <- expand.grid(X = seq(0, 1, length.out = 50),
    Y = seq(0, 1, length.out = 50), category = unique(s$category),
    year = unique(s$year))
  p1 <- predict(fit, newdata = NULL)
  p2 <- predict(fit, newdata = s)
  plot(p1$est, p2$est)
  expect_equal(p1$est, p2$est)
  plot(p2$est, s$mu)
  expect_gt(stats::cor(p2$est, s$mu), 0.95)

  p_temp <- predict(fit, newdata = nd, return_tmb_report = TRUE)

  p <- predict(fit, newdata = nd)
  p$upsilon_stc <- p_temp$proj_upsilon_st_A_vec

  est <- as.list(fit$sd_report, "Estimate", report = TRUE)
  se <- as.list(fit$sd_report, "Std. Error", report = TRUE)
  est$log_sigma_U
  est$sigma_U
  se$log_sigma_U

  fit <- sdmTMB(
    observed ~ 0 + as.factor(category),
    data = s,
    time = "year",
    group = "category",
    spatial = "off",
    spatiotemporal = "rw",
    silent = FALSE,
    mesh = mesh
  )

  nd2 <- nd
  nd2 <- subset(nd2, category %in% c(2, 3))
  p2 <- predict(fit, newdata = nd2)

  pp <- subset(p, category %in% c(2, 3))
  plot(pp$est, p2$est)
  expect_equal(pp$est, p2$est)

  nd3 <- nd
  nd3$category <- factor(nd3$category, levels = c("1", "2", "3", "4"))
  expect_error(p3 <- predict(fit, newdata = nd3), "Extra")

  # ggplot(p, aes(X, Y, colour = upsilon_stc)) +
  #   geom_point() +
  #   facet_grid(category ~ year) +
  #   scale_colour_gradient2()
  #
  # ggplot(s, aes(X, Y, colour = mu)) +
  #   geom_point() +
  #   facet_grid(category ~ year) +
  #   scale_colour_viridis_c()
  #
  # ggplot(p, aes(X, Y, fill = est)) +
  #   geom_raster() +
  #   facet_grid(category ~ year) +
  #   scale_fill_viridis_c()

})

