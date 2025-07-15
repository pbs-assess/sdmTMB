if (requireNamespace("tinyVAST", quietly = TRUE)) {
  TOL <- 1e-3 # Tolerance for parameter comparisons

  test_that("distribution_column functionality works for mixed presence-absence and zero-inflated biomass data", {
    skip_on_cran()

    set.seed(123)

    # Create mixed dataset with presence-absence and biomass data
    n_obs <- 200
    mixed_data <- data.frame(
      x = runif(n_obs, 0, 10),
      y = runif(n_obs, 0, 10),
      depth = rnorm(n_obs, 50, 15),
      year = sample(2015:2020, n_obs, replace = TRUE),
      data_type = rep(c("presence_absence", "biomass"), each = n_obs / 2)
    )

    # Generate response data
    # Presence-absence: binomial with cloglog link
    # Biomass: zero-inflated gamma (delta_gamma)
    mixed_data$response <- ifelse(
      mixed_data$data_type == "presence_absence",
      rbinom(n_obs, 1, plogis(-2 + 0.05 * mixed_data$depth)),
      ifelse(rbinom(n_obs, 1, 0.5) == 1,
        rgamma(n_obs, shape = 2, rate = 1 / exp(1.5 + 0.02 * mixed_data$depth)), 0
      )
    )

    # Add variable column for tinyVAST compatibility
    mixed_data$var <- "logdens"

    # Define family list for mixed distributions
    family_list <- list(
      "presence_absence" = binomial(link = "cloglog"),
      "biomass" = delta_gamma(type = "poisson-link")
    )

    # Fit model with sdmTMB using distribution_column
    fit_sdmTMB <- sdmTMB(
      formula = response ~ depth + factor(year),
      data = mixed_data,
      family = family_list,
      distribution_column = "data_type",
      spatial = "off"
    )

    # Check that model converged successfully
    expect_true(fit_sdmTMB$model$convergence == 0)
    expect_true(fit_sdmTMB$pos_def_hessian)

    # Fit equivalent model with tinyVAST for comparison
    suppressWarnings({
      fit_tinyVAST <- tinyVAST::tinyVAST(
        data = mixed_data,
        formula = response ~ depth + factor(year),
        delta_options = list(formula = ~ depth + factor(year)),
        distribution_column = "data_type",
        family = family_list,
        variable_column = "var"
      )
    })

    # Extract parameters for comparison
    sdmTMB_pars <- fit_sdmTMB$model$par
    tinyVAST_pars <- fit_tinyVAST$opt$par

    # Check that models have reasonable parameter estimates
    expect_true(all(is.finite(sdmTMB_pars)))

    # Compare log-likelihoods
    ll_sdmTMB <- logLik(fit_sdmTMB)
    ll_tinyVAST <- logLik(fit_tinyVAST)
    expect_equal(as.numeric(ll_sdmTMB), as.numeric(ll_tinyVAST), tolerance = TOL)

    # Extract and compare coefficients using tidy methods
    coef_sdmTMB <- tidy(fit_sdmTMB, silent = TRUE)

    # Test that coefficients are reasonable
    expect_true(nrow(coef_sdmTMB) > 0)
    expect_true(all(is.finite(coef_sdmTMB$estimate)))
    expect_true(all(is.finite(coef_sdmTMB$std.error)))

    # Test prediction functionality
    predict_grid <- expand.grid(
      depth = seq(20, 80, length.out = 10),
      year = 2018
    )

    # Make predictions with sdmTMB
    pred_sdmTMB <- predict(fit_sdmTMB, newdata = predict_grid)

    # Check predictions are reasonable
    expect_true(all(is.finite(pred_sdmTMB$est)))
    expect_true(nrow(pred_sdmTMB) == nrow(predict_grid))

    print(fit_sdmTMB)
  })

  test_that("distribution_column functionality works for mixed binomial, delta_gamma, and Poisson families", {
    skip_on_cran()

    set.seed(789)

    # Create mixed dataset with three family types
    n_obs <- 20
    mixed_data_1 <- data.frame(
      depth = rnorm(n_obs, 50, 15),
      # year = sample(2015:2020, n_obs, replace = TRUE),
      # data_type = rep(c("presence_absence", "biomass", "count_data"), each = n_obs / 3)
      data_type = rep(c("biomass", "count_data"), each = n_obs / 2)
    )

    # Generate response data for three families
    # Presence-absence: binomial with cloglog link
    # Biomass: zero-inflated gamma (delta_gamma)
    # Count data: Poisson with log link
    mixed_data_1$response <- ifelse(
      mixed_data_1$data_type == "presence_absence",
      rbinom(n_obs, 1, plogis(-2 + 0.05 * mixed_data_1$depth)),
      ifelse(mixed_data_1$data_type == "biomass",
        ifelse(rbinom(n_obs, 1, 0.5) == 1,
          rgamma(n_obs, shape = 2, rate = 1 / exp(1.5 + 0.02 * mixed_data_1$depth)), 0
        ),
        rpois(n_obs, lambda = exp(0.5 + 0.02 * mixed_data_1$depth))
      )
    )

    # Define family list for three distributions
    family_list_1 <- list(
      "biomass" = delta_gamma(type = "standard"),
      "count_data" = poisson(link = "log")
    )

    # Fit model with sdmTMB using distribution_column
    fit_1_families <- sdmTMB(
      formula = response ~ 1,
      data = mixed_data_1,
      family = family_list_1,
      distribution_column = "data_type",
      spatial = "off"
    )

    suppressWarnings({
      fit_tinyVAST <- tinyVAST::tinyVAST(
        data = mixed_data_1,
        formula = response ~ 1,
        delta_options = list(formula = ~1),
        distribution_column = "data_type",
        family = family_list_1,
        variable_column = "var"
      )
    })

    ll_sdmTMB <- logLik(fit_1_families)
    ll_tinyVAST <- logLik(fit_tinyVAST)
    expect_equal(as.numeric(ll_sdmTMB), as.numeric(ll_tinyVAST), tolerance = TOL)


    rs <- fit_1_families$tmb_obj$report()
    rs$eta_i

    rt <- fit_tinyVAST$obj$report()
    rt$p_i
    rt$p2_i

    p <- fit_tinyVAST$opt$par
    rs <- fit_1_families$tmb_obj$report(par = p)
    rs$eta_i
    rt$p_i
    rt$p2_i

    # now try with binomial + delta_gamma:

    mixed_data_2 <- data.frame(
      depth = rnorm(n_obs, 50, 15),
      # year = sample(2015:2020, n_obs, replace = TRUE),
      # data_type = rep(c("presence_absence", "biomass", "count_data"), each = n_obs / 3)
      data_type = rep(c("biomass", "presence_absence"), each = n_obs / 2)
    )

    mixed_data_2$response <- ifelse(
      mixed_data_2$data_type == "presence_absence",
      rbinom(n_obs, 1, plogis(-2 + 0.05 * mixed_data_1$depth)),
      ifelse(mixed_data_1$data_type == "biomass",
        ifelse(rbinom(n_obs, 1, 0.5) == 1,
          rgamma(n_obs, shape = 2, rate = 1 / exp(1.5 + 0.02 * mixed_data_1$depth)), 0
        ),
        rpois(n_obs, lambda = exp(0.5 + 0.02 * mixed_data_1$depth))
      )
    )

    family_list_2 <- list(
      "biomass" = delta_gamma(type = "standard"),
      "presence_absence" = binomial(link = "cloglog")
    )

    fit_2_families <- sdmTMB(
      formula = response ~ 1,
      data = mixed_data_2,
      family = family_list_2,
      distribution_column = "data_type",
      spatial = "off"
    )

    fit_tinyVAST <- tinyVAST::tinyVAST(
      data = mixed_data_2,
      formula = response ~ 1,
      delta_options = list(formula = ~1),
      distribution_column = "data_type",
      family = family_list_2,
      variable_column = "var"
    )

    ll_sdmTMB <- logLik(fit_2_families)
    ll_tinyVAST <- logLik(fit_tinyVAST)
    expect_equal(as.numeric(ll_sdmTMB), as.numeric(ll_tinyVAST), tolerance = TOL)

    rs <- fit_2_families$tmb_obj$report()
    rs$eta_i

    rt <- fit_tinyVAST$obj$report()
    rt$p_i
    rt$p2_i

    # just poisson:
    mixed_data_3 <- data.frame(
      depth = rnorm(n_obs, 50, 15),
      data_type = rep(c("count"), each = n_obs / 2)
    )
    mixed_data_3$response <- ifelse(
      mixed_data_3$data_type == "presence_absence",
      rbinom(n_obs, 1, plogis(-2 + 0.05 * mixed_data_3$depth)),
      ifelse(mixed_data_3$data_type == "biomass",
        ifelse(rbinom(n_obs, 1, 0.5) == 1,
          rgamma(n_obs, shape = 2, rate = 1 / exp(1.5 + 0.02 * mixed_data_3$depth)), 0
        ),
        rpois(n_obs, lambda = exp(0.5 + 0.02 * mixed_data_3$depth))
      )
    )

    family_list_3 <- list(
      "count" = poisson()
    )

    fit_3_families <- sdmTMB(
      formula = response ~ 1,
      data = mixed_data_3,
      family = family_list_3,
      distribution_column = "data_type",
      spatial = "off"
    )

    fit_tinyVAST <- tinyVAST::tinyVAST(
      data = mixed_data_3,
      formula = response ~ 1,
      delta_options = list(formula = ~1),
      distribution_column = "data_type",
      family = family_list_3,
      variable_column = "var"
    )

    ll_sdmTMB <- logLik(fit_3_families)
    ll_tinyVAST <- logLik(fit_tinyVAST)
    expect_equal(as.numeric(ll_sdmTMB), as.numeric(ll_tinyVAST), tolerance = TOL)
  })
}

