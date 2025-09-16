test_that("Zero inflation works and matches glmmTMB", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")

  set.seed(123)

  # Create test data with zero inflation
  n <- 100
  data <- data.frame(
    x = runif(n, -1, 1),
    y = runif(n, -1, 1),
    depth = rnorm(n)
  )

  # Simulate zero-inflated Poisson data
  zi_prob <- 0.3
  lambda <- exp(1 + 0.5 * data$depth)
  zi_draws <- rbinom(n, 1, zi_prob)
  pois_draws <- rpois(n, lambda)
  data$count <- ifelse(zi_draws == 1, 0, pois_draws)

  # Fit with glmmTMB (reference)
  fit_glmmTMB <- glmmTMB::glmmTMB(count ~ depth,
    data = data,
    ziformula = ~1,
    family = poisson()
  )

  # Fit with sdmTMB
  fit_sdmTMB <- sdmTMB(count ~ depth,
    data = data,
    zi = TRUE,
    spatial = "off",
    family = poisson()
  )

  # Extract zero inflation probabilities
  glmmTMB_zi <- as.numeric(plogis(glmmTMB::fixef(fit_glmmTMB)$zi[1]))
  sdmTMB_zi <- tidy(fit_sdmTMB, "fixed")$estimate[tidy(fit_sdmTMB, "fixed")$term == "zi_prob"]

  # Extract log-likelihoods
  glmmTMB_loglik <- as.numeric(logLik(fit_glmmTMB))
  sdmTMB_loglik <- -fit_sdmTMB$model$objective

  # Test that estimates match
  expect_equal(sdmTMB_zi, glmmTMB_zi, tolerance = 0.01)
  expect_equal(sdmTMB_loglik, glmmTMB_loglik, tolerance = 0.01)
})

test_that("Zero inflation parameter is not estimated when zi = FALSE", {
  set.seed(123)
  data <- data.frame(
    x = runif(50, -1, 1),
    y = runif(50, -1, 1),
    count = rpois(50, 2)
  )
  fit <- sdmTMB(count ~ 1,
    data = data,
    zi = FALSE,
    spatial = "off",
    family = poisson()
  )
  expect_false("zi_p" %in% names(fit$sd_report$value))
})

test_that("Zero inflation cannot be combined with delta models", {
  data <- data.frame(
    x = runif(50, -1, 1),
    y = runif(50, -1, 1),
    count = rpois(50, 2)
  )
  expect_error(
    sdmTMB(count ~ 1,
      data = data,
      spatial = "off",
      zi = TRUE,
      family = delta_gamma()
    ),
    "Zero inflation.*cannot be combined with delta"
  )
})

test_that("Zero inflation parameter appears in tidy fixed effects", {
  set.seed(123)
  data <- data.frame(
    x = runif(50, -1, 1),
    y = runif(50, -1, 1),
    count = c(rep(0, 15), rpois(35, 3))  # some zeros plus Poisson data
  )

  fit <- sdmTMB(count ~ 1,
    data = data,
    zi = TRUE,
    spatial = "off",
    family = poisson()
  )

  # Check that zi_prob is in fixed effects
  tidy_fixed <- tidy(fit, "fixed")
  expect_true("zi_prob" %in% tidy_fixed$term)

  # Check that zi_prob is NOT in ran_pars
  tidy_ran <- tidy(fit, "ran_pars")
  expect_false("zi_prob" %in% tidy_ran$term)

  # Check that zi_prob has reasonable values
  zi_row <- tidy_fixed[tidy_fixed$term == "zi_prob", ]
  expect_true(zi_row$estimate > 0 & zi_row$estimate < 1)  # probability
  expect_true(zi_row$std.error > 0)  # positive standard error
})

test_that("Zero inflation parameter appears in print output", {
  set.seed(123)
  data <- data.frame(
    x = runif(50, -1, 1),
    y = runif(50, -1, 1),
    count = c(rep(0, 15), rpois(35, 3))  # some zeros plus Poisson data
  )

  fit <- sdmTMB(count ~ 1,
    data = data,
    zi = TRUE,
    spatial = "off",
    family = poisson()
  )

  # Capture print output
  print_output <- capture.output(print(fit))
  print_text <- paste(print_output, collapse = " ")

  # Check that zi_prob appears in the print output
  expect_true(grepl("zi_prob", print_text))
})
