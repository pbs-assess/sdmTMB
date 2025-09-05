test_that("make_category_svc works with basic data", {
  # Create test data
  set.seed(123)
  data <- data.frame(
    age = factor(rep(1:3, each = 20)),
    year = rep(2020:2022, 20), 
    abundance = rnorm(60),
    x = runif(60), 
    y = runif(60)
  )
  
  # Basic setup with shared SDs
  setup <- make_category_svc(
    data = data,
    category_column = "age", 
    time_column = "year",
    share_spatial_sd = TRUE,
    share_spatiotemporal_sd = TRUE
  )
  
  # Check return structure
  expect_type(setup, "list")
  expect_named(setup, c("data_expanded", "svc_formula", "svc_map", "info"))
  
  # Check data expansion
  expect_true(nrow(setup$data_expanded) == nrow(data))
  expect_true(ncol(setup$data_expanded) > ncol(data))
  
  # Check formula
  expect_s3_class(setup$svc_formula, "formula")
  
  # Check map structure  
  expect_type(setup$svc_map, "list")
  expect_named(setup$svc_map, "ln_tau_Z")
  expect_s3_class(setup$svc_map$ln_tau_Z, "factor")
  
  # Check info content
  expect_type(setup$info, "list")
  expect_equal(setup$info$n_categories, 3)
  expect_equal(setup$info$n_times, 3)
  expect_equal(setup$info$n_variance_parameters, 2) # shared spatial + shared spatiotemporal
})

test_that("make_category_svc works with different sharing options", {
  # Create small test data
  data <- data.frame(
    age = factor(rep(1:2, each = 4)),
    year = rep(2020:2021, 4), 
    abundance = rnorm(8),
    x = runif(8), 
    y = runif(8)
  )
  
  # Test all sharing combinations
  setup1 <- make_category_svc(data, "age", "year", TRUE, TRUE)   # both shared
  setup2 <- make_category_svc(data, "age", "year", TRUE, FALSE)  # spatial shared only
  setup3 <- make_category_svc(data, "age", "year", FALSE, TRUE)  # spatiotemporal shared only
  setup4 <- make_category_svc(data, "age", "year", FALSE, FALSE) # none shared
  
  # Check number of variance parameters
  expect_equal(setup1$info$n_variance_parameters, 2)  # 1 spatial + 1 spatiotemporal
  expect_equal(setup2$info$n_variance_parameters, 5)  # 1 spatial + 4 spatiotemporal
  expect_equal(setup3$info$n_variance_parameters, 3)  # 2 spatial + 1 spatiotemporal
  expect_equal(setup4$info$n_variance_parameters, 6)  # 2 spatial + 4 spatiotemporal
  
  # Check map factor levels
  expect_equal(max(as.numeric(setup1$svc_map$ln_tau_Z)), 2)
  expect_equal(max(as.numeric(setup2$svc_map$ln_tau_Z)), 5)
  expect_equal(max(as.numeric(setup3$svc_map$ln_tau_Z)), 3)
  expect_equal(max(as.numeric(setup4$svc_map$ln_tau_Z)), 6)
})

test_that("make_category_svc validates inputs", {
  data <- data.frame(
    age = factor(1:2),
    year = 2020:2021,
    x = 1:2
  )
  
  # Valid case should work
  expect_no_error(
    make_category_svc(data, "age", "year")
  )
  
  # Invalid data frame
  expect_error(
    make_category_svc("not_a_df", "age", "year"),
    "data is not a data frame"
  )
  
  # Missing column
  expect_error(
    make_category_svc(data, "missing_col", "year"),
    "Column 'missing_col' not found in data"
  )
  
  expect_error(
    make_category_svc(data, "age", "missing_col"),
    "Column 'missing_col' not found in data"
  )
  
  # Invalid logical arguments
  expect_error(
    make_category_svc(data, "age", "year", share_spatial_sd = "not_logical"),
    "is.logical\\(share_spatial_sd\\) is not TRUE"
  )
  
  expect_error(
    make_category_svc(data, "age", "year", share_spatiotemporal_sd = 1),
    "is.logical\\(share_spatiotemporal_sd\\) is not TRUE"
  )
})

test_that("make_category_svc handles different category types", {
  # Test with different category column types
  data1 <- data.frame(
    cat = factor(rep(letters[1:3], 2)),  # factor categories
    year = rep(2020:2021, 3),
    x = 1:6
  )
  
  data2 <- data.frame(
    cat = rep(letters[1:3], 2),  # character categories  
    year = rep(2020:2021, 3),
    x = 1:6
  )
  
  setup1 <- make_category_svc(data1, "cat", "year")
  setup2 <- make_category_svc(data2, "cat", "year")
  
  # Both should work and produce valid outputs
  expect_type(setup1, "list")
  expect_type(setup2, "list")
  expect_equal(setup1$info$n_categories, 3)
  expect_equal(setup2$info$n_categories, 3)
})

test_that("make_category_svc creates correct model matrices", {
  data <- data.frame(
    age = factor(rep(1:2, 2)),
    year = rep(2020:2021, each = 2),
    x = 1:4
  )
  
  setup <- make_category_svc(data, "age", "year")
  
  # Check that spatial terms match expected structure
  expected_spatial_cols <- paste0("age", 1:2)
  spatial_cols <- setup$info$spatial_terms
  expect_equal(sort(spatial_cols), sort(expected_spatial_cols))
  
  # Check that spatiotemporal terms include year-age interactions
  spatiotemporal_cols <- setup$info$spatiotemporal_terms
  expect_true(length(spatiotemporal_cols) == 4) # 2 years × 2 ages
  expect_true(all(grepl(":", spatiotemporal_cols))) # should contain interactions
})