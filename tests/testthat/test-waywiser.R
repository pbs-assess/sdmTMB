test_that("cv_to_waywiser() works for basic models", {
  skip_on_cran()
  skip_if_not_installed("sf")

  d <- pcod_2011
  mesh <- make_mesh(d, c("X", "Y"), cutoff = 10)

  set.seed(123)
  m_cv <- sdmTMB_cv(
    density ~ s(depth_scaled),
    data = d,
    mesh = mesh,
    family = tweedie(),
    spatial = "off",
    k_folds = 2
  )

  # Convert to sf
  cv_sf <- cv_to_waywiser(m_cv, ll_names = c("lon", "lat"))

  # Check that it's an sf object
  expect_s3_class(cv_sf, "sf")

  # Check that it has the right columns
  expect_true("truth" %in% names(cv_sf))
  expect_true("estimate" %in% names(cv_sf))
  expect_true("geometry" %in% names(cv_sf))

  # Check dimensions
  expect_equal(nrow(cv_sf), nrow(d))

  # Check that truth and estimate are numeric
  expect_type(cv_sf$truth, "double")
  expect_type(cv_sf$estimate, "double")

  # Check that geometry is POINT
  expect_true(inherits(sf::st_geometry(cv_sf), "sfc_POINT"))
})

test_that("cv_to_waywiser() works for delta models", {
  skip_on_cran()
  skip_if_not_installed("sf")

  d <- pcod_2011
  mesh <- make_mesh(d, c("X", "Y"), cutoff = 10)

  set.seed(456)
  m_cv_delta <- sdmTMB_cv(
    density ~ s(depth_scaled),
    data = d,
    mesh = mesh,
    family = delta_gamma(),
    spatial = "off",
    k_folds = 2
  )

  # Convert delta model (returns combined predictions)
  cv_sf <- cv_to_waywiser(m_cv_delta, ll_names = c("lon", "lat"))
  expect_s3_class(cv_sf, "sf")
  expect_true("truth" %in% names(cv_sf))
  expect_true("estimate" %in% names(cv_sf))
  expect_equal(nrow(cv_sf), nrow(d))

  # Truth should be the original response variable
  expect_equal(cv_sf$truth, d$density)
  # Estimate should be the combined CV predictions
  expect_equal(cv_sf$estimate, m_cv_delta$data$cv_predicted)
})

test_that("cv_to_waywiser() validates column names", {
  skip_on_cran()
  skip_if_not_installed("sf")

  d <- pcod_2011
  mesh <- make_mesh(d, c("X", "Y"), cutoff = 10)

  set.seed(202)
  m_cv <- sdmTMB_cv(
    density ~ s(depth_scaled),
    data = d,
    mesh = mesh,
    family = tweedie(),
    spatial = "off",
    k_folds = 2
  )

  # Invalid ll_names should error
  expect_error(
    cv_to_waywiser(m_cv, ll_names = c("longitude", "latitude")),
    "not found in CV data"
  )
})

test_that("cv_to_waywiser() accepts manual UTM CRS specification", {
  skip_on_cran()
  skip_if_not_installed("sf")

  d <- pcod_2011
  mesh <- make_mesh(d, c("X", "Y"), cutoff = 10)

  set.seed(303)
  m_cv <- sdmTMB_cv(
    density ~ s(depth_scaled),
    data = d,
    mesh = mesh,
    family = tweedie(),
    spatial = "off",
    k_folds = 2
  )

  # Manually specify UTM CRS
  cv_sf <- cv_to_waywiser(m_cv, ll_names = c("lon", "lat"), utm_crs = 32609)

  expect_s3_class(cv_sf, "sf")
  expect_equal(sf::st_crs(cv_sf)$epsg, 32609)
})

