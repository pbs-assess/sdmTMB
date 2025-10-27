test_that("build_family_arrays works for single Gaussian family", {
  n_obs <- 10
  data <- data.frame(y = rnorm(n_obs), x = rnorm(n_obs))

  result <- build_family_arrays(
    family = gaussian(),
    distribution_column = NULL,
    data = data,
    n_obs = n_obs
  )

  # Check dimensions
  expect_equal(dim(result$family), c(n_obs, 1))
  expect_equal(dim(result$link), c(n_obs, 1))
  expect_equal(length(result$e_i), n_obs)

  # Check all observations have same family/link
  expect_true(all(result$family[, 1] == .valid_family["gaussian"]))
  expect_true(all(result$link[, 1] == .valid_link["identity"]))

  # Check e_i maps all to group 0
  expect_true(all(result$e_i == 0))

  # Check component usage
  expect_true(all(result$component_usage[, 1] == 1))

  # Check metadata
  expect_equal(result$n_m_effective, 1)
  expect_equal(result$n_families, 1)
  expect_equal(result$components, 1)
})

test_that("build_family_arrays works for single delta family", {
  n_obs <- 10
  data <- data.frame(y = c(0, 0, 1, 2, 0, 3, 0, 1, 0, 2), x = rnorm(n_obs))

  result <- build_family_arrays(
    family = delta_gamma(),
    distribution_column = NULL,
    data = data,
    n_obs = n_obs
  )

  # Check dimensions for delta model
  expect_equal(dim(result$family), c(n_obs, 2))
  expect_equal(dim(result$link), c(n_obs, 2))

  # Check family codes for delta
  expect_true(all(result$family[, 1] == .valid_family["binomial"]))
  expect_true(all(result$family[, 2] == .valid_family["Gamma"]))

  # Check link codes
  expect_true(all(result$link[, 1] == .valid_link["logit"]))
  expect_true(all(result$link[, 2] == .valid_link["log"]))

  # Check component usage - both components active
  expect_true(all(result$component_usage[, 1] == 1))
  expect_true(all(result$component_usage[, 2] == 1))

  # Check metadata
  expect_equal(result$n_m_effective, 2)
  expect_equal(result$components, 2)
})

test_that("build_family_arrays works for integrated model with non-delta families", {
  n_obs <- 20
  data <- data.frame(
    y = rpois(n_obs, 5),
    x = rnorm(n_obs),
    data_type = rep(c("binomial_data", "poisson_data"), each = 10)
  )

  family_list <- list(
    "binomial_data" = binomial(link = "cloglog"),
    "poisson_data" = poisson(link = "log")
  )

  result <- build_family_arrays(
    family = family_list,
    distribution_column = "data_type",
    data = data,
    n_obs = n_obs
  )

  # Check dimensions
  expect_equal(dim(result$family), c(n_obs, 1))

  # Check first 10 are binomial
  expect_true(all(result$family[1:10, 1] == .valid_family["binomial"]))
  expect_true(all(result$link[1:10, 1] == .valid_link["cloglog"]))
  expect_true(all(result$e_i[1:10] == 0))

  # Check last 10 are poisson
  expect_true(all(result$family[11:20, 1] == .valid_family["poisson"]))
  expect_true(all(result$link[11:20, 1] == .valid_link["log"]))
  expect_true(all(result$e_i[11:20] == 1))

  # Check metadata
  expect_equal(result$n_families, 2)
  expect_equal(result$n_m_effective, 1)
  expect_equal(as.integer(result$components), c(1, 1))
})

test_that("build_family_arrays works for integrated model with mixed delta/non-delta", {
  n_obs <- 30
  data <- data.frame(
    y = c(rbinom(10, 1, 0.5), rgamma(20, 2, 1)),
    x = rnorm(n_obs),
    data_type = c(rep("presence_absence", 10), rep("biomass", 20))
  )

  family_list <- list(
    "presence_absence" = binomial(link = "cloglog"),
    "biomass" = delta_gamma()
  )

  result <- build_family_arrays(
    family = family_list,
    distribution_column = "data_type",
    data = data,
    n_obs = n_obs
  )

  # Check dimensions - should be n_m = 2 (max components)
  expect_equal(dim(result$family), c(n_obs, 2))
  expect_equal(result$n_m_effective, 2)

  # Check first 10 (binomial) - only component 1 active
  expect_true(all(result$family[1:10, 1] == .valid_family["binomial"]))
  expect_true(all(is.na(result$family[1:10, 2])))
  expect_true(all(result$component_usage[1:10, 1] == 1))
  expect_true(all(result$component_usage[1:10, 2] == 0))

  # Check last 20 (delta_gamma) - both components active
  expect_true(all(result$family[11:30, 1] == .valid_family["binomial"]))
  expect_true(all(result$family[11:30, 2] == .valid_family["Gamma"]))
  expect_true(all(result$component_usage[11:30, 1] == 1))
  expect_true(all(result$component_usage[11:30, 2] == 1))

  # Check e_i mapping
  expect_true(all(result$e_i[1:10] == 0))  # binomial -> group 0
  expect_true(all(result$e_i[11:30] == 1)) # delta_gamma -> group 1

  # Check components vector
  expect_equal(as.integer(result$components), c(1, 2))
})

test_that("build_family_arrays detects invalid distribution_column", {
  data <- data.frame(y = rnorm(10), x = rnorm(10), type = rep("A", 10))

  family_list <- list("B" = gaussian(), "C" = tweedie())

  expect_error(
    build_family_arrays(family_list, "type", data, 10),
    "don't match"
  )
})

test_that("build_family_arrays detects missing distribution_column", {
  data <- data.frame(y = rnorm(10), x = rnorm(10))

  family_list <- list("A" = gaussian())

  expect_error(
    build_family_arrays(family_list, "nonexistent", data, 10),
    "not found in data"
  )
})

test_that("build_family_arrays requires named list for integrated models", {
  data <- data.frame(y = rnorm(10), type = rep("A", 10))

  # Unnamed list should error
  expect_error(
    build_family_arrays(list(gaussian(), tweedie()), "type", data, 10),
    "must be a named list"
  )
})

test_that("build_family_arrays handles poisson_link_delta correctly", {
  n_obs <- 10
  data <- data.frame(y = rpois(n_obs, 5), x = rnorm(n_obs))

  result <- build_family_arrays(
    family = delta_gamma(type = "poisson-link"),
    distribution_column = NULL,
    data = data,
    n_obs = n_obs
  )

  # All observations should have poisson_link_delta = 1
  expect_true(all(result$poisson_link_delta == 1))

  # Compare to regular delta_gamma
  result2 <- build_family_arrays(
    family = delta_gamma(),
    distribution_column = NULL,
    data = data,
    n_obs = n_obs
  )

  # Regular delta should have poisson_link_delta = 0
  expect_true(all(result2$poisson_link_delta == 0))
})

test_that("build_family_arrays handles three different families", {
  n_obs <- 30
  data <- data.frame(
    y = rnorm(n_obs),
    data_type = rep(c("type1", "type2", "type3"), each = 10)
  )

  family_list <- list(
    "type1" = gaussian(),
    "type2" = poisson(),
    "type3" = nbinom2()
  )

  result <- build_family_arrays(
    family = family_list,
    distribution_column = "data_type",
    data = data,
    n_obs = n_obs
  )

  # Check we have 3 unique families
  expect_equal(result$n_families, 3)

  # Check e_i maps correctly (0, 1, 2)
  expect_true(all(result$e_i[1:10] == 0))
  expect_true(all(result$e_i[11:20] == 1))
  expect_true(all(result$e_i[21:30] == 2))

  # Check family codes
  expect_true(all(result$family[1:10, 1] == .valid_family["gaussian"]))
  expect_true(all(result$family[11:20, 1] == .valid_family["poisson"]))
  expect_true(all(result$family[21:30, 1] == .valid_family["nbinom2"]))
})
