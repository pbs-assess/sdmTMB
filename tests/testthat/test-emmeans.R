test_that("emmeans works with parametric models", {
  skip_if_not_installed("emmeans")
  library(emmeans)
  
  # Basic parametric model without spatial component  
  fit <- sdmTMB(
    present ~ as.factor(year),
    data = pcod_2011,
    spatial = "off",
    family = binomial()
  )
  
  # Test emmeans creation
  emm <- emmeans::emmeans(fit, ~ year)
  expect_s4_class(emm, "emmGrid")
  
  # Test pairwise comparisons
  pairs_result <- pairs(emm)  # pairs is a method, not namespaced
  expect_s4_class(pairs_result, "emmGrid")
  
  # Test response scale
  emm_resp <- emmeans::emmeans(fit, ~ year, type = "response")
  expect_s4_class(emm_resp, "emmGrid")
  
  # Check that results have correct dimensions
  emm_summary <- summary(emm)
  expect_equal(nrow(emm_summary), length(unique(pcod_2011$year)))
  expect_true(all(c("year", "emmean", "SE", "df", "lower.CL", "upper.CL") %in% names(emm_summary)))
})

test_that("emmeans works with continuous covariates", {
  skip_if_not_installed("emmeans")
  library(emmeans)
  
  fit <- sdmTMB(
    present ~ as.factor(year) + depth_scaled,
    data = pcod_2011,
    spatial = "off",
    family = binomial()
  )
  
  # Test emmeans with continuous covariate
  emm <- emmeans::emmeans(fit, ~ year)
  expect_s4_class(emm, "emmGrid")
  
  # Test emtrends for slopes
  emtrends_result <- emmeans::emtrends(fit, ~ year, var = "depth_scaled")
  expect_s4_class(emtrends_result, "emmGrid")
  
  # Test interaction
  fit_int <- sdmTMB(
    present ~ as.factor(year) * depth_scaled,
    data = pcod_2011,
    spatial = "off",
    family = binomial()
  )
  
  emm_int <- emmeans::emmeans(fit_int, ~ year)
  expect_s4_class(emm_int, "emmGrid")
})

test_that("emmeans works with smoothers", {
  skip_if_not_installed("emmeans")
  library(emmeans)
  
  # Model with smoother - this was previously broken
  fit_smooth <- sdmTMB(
    present ~ as.factor(year) + s(depth_scaled),
    data = pcod_2011,
    spatial = "off",
    family = binomial()
  )
  
  # Test that emmeans works with smoothers
  emm_smooth <- emmeans::emmeans(fit_smooth, ~ year)
  expect_s4_class(emm_smooth, "emmGrid")
  
  # Test pairwise comparisons with smoothers
  pairs_smooth <- pairs(emm_smooth)
  expect_s4_class(pairs_smooth, "emmGrid")
  
  # Test response scale with smoothers  
  emm_smooth_resp <- emmeans::emmeans(fit_smooth, ~ year, type = "response")
  expect_s4_class(emm_smooth_resp, "emmGrid")
  
  # Check that results are reasonable
  emm_summary <- summary(emm_smooth)
  expect_equal(nrow(emm_summary), length(unique(pcod_2011$year)))
  expect_true(all(is.finite(emm_summary$emmean)))
  expect_true(all(emm_summary$SE > 0))
})

test_that("emmeans works with multiple smoothers", {
  skip_if_not_installed("emmeans")
  library(emmeans)
  
  # Model with multiple smoothers
  fit_multi_smooth <- sdmTMB(
    present ~ as.factor(year) + s(depth_scaled) + s(X),
    data = pcod_2011,
    spatial = "off",
    family = binomial()
  )
  
  # Test that emmeans still works with multiple smoothers
  emm_multi <- emmeans::emmeans(fit_multi_smooth, ~ year)
  expect_s4_class(emm_multi, "emmGrid")
  
  # Test pairwise comparisons
  pairs_multi <- pairs(emm_multi)
  expect_s4_class(pairs_multi, "emmGrid")
  
  # Results should be sensible
  emm_summary <- summary(emm_multi)
  expect_true(all(is.finite(emm_summary$emmean)))
})

test_that("emmeans results are consistent between parametric and smooth models", {
  skip_if_not_installed("emmeans")
  library(emmeans)
  
  # Parametric model  
  fit_param <- sdmTMB(
    present ~ as.factor(year) + depth_scaled,
    data = pcod_2011,
    spatial = "off",
    family = binomial()
  )
  
  # Smooth model
  fit_smooth <- sdmTMB(
    present ~ as.factor(year) + s(depth_scaled),
    data = pcod_2011,
    spatial = "off", 
    family = binomial()
  )
  
  # Get emmeans for both
  emm_param <- emmeans::emmeans(fit_param, ~ year)
  emm_smooth <- emmeans::emmeans(fit_smooth, ~ year)
  
  # Results should be reasonably similar (smooth vs linear depth effect)
  param_est <- summary(emm_param)$emmean
  smooth_est <- summary(emm_smooth)$emmean
  
  # They shouldn't be identical but should be correlated
  expect_true(cor(param_est, smooth_est) > 0.8)
})

test_that("emmeans works with delta models", {
  skip_if_not_installed("emmeans")
  library(emmeans)

  # Delta model
  fit_delta <- sdmTMB(
    density ~ as.factor(year),
    data = pcod_2011,
    spatial = "off",
    family = delta_gamma()
  )

  # Test binomial component (model = 1)
  emm1 <- emmeans::emmeans(fit_delta, ~ year, model = 1)
  expect_s4_class(emm1, "emmGrid")

  # Test positive component (model = 2)
  emm2 <- emmeans::emmeans(fit_delta, ~ year, model = 2)
  expect_s4_class(emm2, "emmGrid")

  # Check that results have correct dimensions
  emm1_summary <- summary(emm1)
  emm2_summary <- summary(emm2)
  expect_equal(nrow(emm1_summary), length(unique(pcod_2011$year)))
  expect_equal(nrow(emm2_summary), length(unique(pcod_2011$year)))

  # Test pairwise comparisons work
  pairs1 <- pairs(emm1)
  pairs2 <- pairs(emm2)
  expect_s4_class(pairs1, "emmGrid")
  expect_s4_class(pairs2, "emmGrid")

  # Test response scale works
  emm1_resp <- emmeans::emmeans(fit_delta, ~ year, model = 1, type = "response")
  emm2_resp <- emmeans::emmeans(fit_delta, ~ year, model = 2, type = "response")
  expect_s4_class(emm1_resp, "emmGrid")
  expect_s4_class(emm2_resp, "emmGrid")
})

test_that("emmeans works with binomial family and smoothers", {
  skip_if_not_installed("emmeans")
  library(emmeans)
  
  # Test binomial with smoother - this is the main use case that was broken
  fit <- sdmTMB(
    present ~ as.factor(year) + s(depth_scaled),
    data = pcod_2011,
    spatial = "off",
    family = binomial()
  )
  
  # Test emmeans works
  emm <- emmeans::emmeans(fit, ~ year)
  expect_s4_class(emm, "emmGrid")
  
  # Test that we can extract meaningful results
  emm_summary <- summary(emm)
  expect_true(nrow(emm_summary) == 4)  # 4 years
  expect_true(all(c("year", "emmean") %in% names(emm_summary)))
})