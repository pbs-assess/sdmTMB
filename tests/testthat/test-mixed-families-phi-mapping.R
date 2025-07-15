test_that("phi mapping works correctly for mixed delta/non-delta families", {
  skip_on_cran()
  
  # Create test data with mixed delta and non-delta families
  set.seed(42)
  n_obs <- 200
  
  # Presence-absence data (binomial - no phi needed)
  pa_data <- data.frame(
    x = runif(n_obs/2, 0, 10),
    y = runif(n_obs/2, 0, 10),
    depth = rnorm(n_obs/2, 50, 10),
    response = rbinom(n_obs/2, 1, 0.3),
    data_type = "presence_absence"
  )
  
  # Biomass data (delta_gamma - phi needed for gamma component)
  biomass_data <- data.frame(
    x = runif(n_obs/2, 0, 10),
    y = runif(n_obs/2, 0, 10),
    depth = rnorm(n_obs/2, 50, 10),
    response = c(rep(0, n_obs/4), rgamma(n_obs/4, shape = 2, rate = 1)),
    data_type = "biomass"
  )
  
  test_data <- rbind(pa_data, biomass_data)
  mesh <- make_mesh(test_data, xy_cols = c("x", "y"), cutoff = 2)
  
  # Family specification
  family_list <- list(
    "presence_absence" = binomial(link = "logit"),
    "biomass" = delta_gamma()
  )
  
  # Fit model without optimization to examine mapping
  fit <- sdmTMB(
    formula = response ~ depth,
    data = test_data,
    mesh = mesh,
    family = family_list,
    distribution_column = "data_type",
    spatial = "off",
    do_fit = FALSE
  )
  
  # Examine the phi mapping
  phi_map <- fit$tmb_map$ln_phi
  
  # Debug output
  cat("\\nPhi mapping:", as.character(phi_map), "\\n")
  cat("Data types:", table(test_data$data_type), "\\n")
  cat("d_i values:", unique(fit$tmb_data$d_i), "\\n")
  
  # Expected behavior (corrected):
  # - d_i assignment: presence_absence=1, biomass=0 (based on factor ordering)
  # - n_m = 2 (forced for mixed case)
  # - phi_map matrix: [1,3; 2,4] then flattened by columns
  # 
  # Mapping should be:
  # Group 1 (d_i=0, biomass):     [1,3] -> [NA,3] (binomial->NA, gamma->keep)  
  # Group 2 (d_i=1, presence):    [2,4] -> [NA,NA] (binomial->NA, unused->NA)
  # Flattened: [NA, NA, 3, NA]
  
  expect_equal(length(phi_map), 4)
  expect_true(is.na(phi_map[1]))   # Group 1, component 1: binomial part of delta (no phi)
  expect_true(is.na(phi_map[2]))   # Group 2, component 1: binomial (no phi)
  expect_false(is.na(phi_map[3]))  # Group 1, component 2: gamma part of delta (needs phi)
  expect_true(is.na(phi_map[4]))   # Group 2, component 2: unused (no phi)
})