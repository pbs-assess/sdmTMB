sim_breakpt <- function(n = 300, cut = 0.3, sd = 0.4, seed = 1) {
  set.seed(seed)
  x <- rnorm(n)
  data.frame(x = x, y = 1 + pmin(x, cut) + rnorm(n, sd = sd))
}

test_that("profile_breakpt() recovers the cutpoint and fits the final model", {
  skip_on_cran()
  d <- sim_breakpt()
  p <- profile_breakpt(y ~ breakpt(x), data = d, spatial = "off")
  e <- p$estimate
  expect_s3_class(p, "breakpt_profile")
  expect_true(e$lwr < 0.3 && e$upr > 0.3)
  expect_true(e$lwr < e$cut && e$cut < e$upr)
  expect_false(e$open_interval)
  expect_false(e$few_below)
  expect_true(e$refine_complete)
  expect_equal(e$n_failed, 0)
  expect_equal(min(p$profile$dnll, na.rm = TRUE), 0)
  expect_false(is.unsorted(p$profile$cut))

  # the final fit holds the cutpoint at the best value and matches the profile
  b <- p$fit$model$par[names(p$fit$model$par) == "b_threshold"]
  expect_length(b, 1L) # only the slope is estimated
  expect_equal(p$fit$tmb_obj$env$parList()$b_threshold[2, 1], e$cut)
  expect_equal(e$final_nll, min(p$profile$nll, na.rm = TRUE), tolerance = 1e-4)
  expect_equal(e$aic, 2 * e$nll + 2 * e$k)
  expect_equal(e$k, length(p$fit$model$par) + 1L)

  # the profile minimum is no worse than fitting the cutpoint directly from the truth
  m <- suppressMessages(suppressWarnings(sdmTMB(y ~ breakpt(x), data = d, spatial = "off",
    control = sdmTMBcontrol(start = list(b_threshold = matrix(c(1, 0.3), 2L, 1L))))))
  expect_lt(e$nll, m$model$objective + 0.01) # within the grid's resolution
  expect_output(print(p), "Cutpoint")
})

test_that("profile_breakpt() warns when the cutpoint is poorly identified", {
  skip_on_cran()
  set.seed(2)
  d <- data.frame(x = rnorm(200))
  d$y <- 1 + d$x + rnorm(200, sd = 0.5) # linear: no bend within the data
  w <- character()
  p <- withCallingHandlers(
    profile_breakpt(y ~ breakpt(x), data = d, spatial = "off", final = FALSE),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(p$estimate$open_interval)
  expect_true(any(grepl("end of the grid", w)))
  expect_null(p$fit)
  expect_output(print(p), "reaches the end of the grid")

  # an open interval: the shaded band runs to the end of the grid
  skip_if_not_installed("ggplot2")
  g <- plot(p)
  expect_s3_class(g, "ggplot")
  band <- ggplot2::layer_data(g, 1)
  expect_equal(band$xmax, max(p$profile$cut))
  expect_equal(band$xmin, p$estimate$lwr)
})

test_that("profile_breakpt() checks its inputs", {
  skip_on_cran()
  d <- sim_breakpt(n = 100)
  expect_error(profile_breakpt(y ~ x, data = d, spatial = "off"), "breakpt")
  expect_error(profile_breakpt(y ~ logistic(x), data = d, spatial = "off"), "breakpt")
  expect_error(profile_breakpt(y ~ breakpt(x), data = d, spatial = "off",
    control = sdmTMBcontrol(map = list(b_threshold = factor(c(1, NA))))), "map")
  expect_error(profile_breakpt(y ~ breakpt(x), data = d, spatial = "off", tol = 0), "positive")
  expect_error(profile_breakpt(y ~ breakpt(x), data = d, spatial = "off", max_fits = Inf), "finite")
  expect_error(profile_breakpt(y ~ breakpt(x), data = d, spatial = "off", cuts = c(0, 1)), "3")
  d$y <- ifelse(d$y > 1, d$y, 0)
  expect_error(profile_breakpt(y ~ breakpt(x), data = d, spatial = "off",
    family = delta_gamma()), "Delta")
})
