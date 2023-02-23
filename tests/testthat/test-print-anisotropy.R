test_that("Print anisotropy prints correctly", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_on_ci()

  # No anisotropy
  fit1 <- sdmTMB(
    data = pcod_2011,
    formula = density ~ 1,
    mesh = pcod_mesh_2011,
    family = tweedie(),
    anisotropy = FALSE
  )

  expect_output(print(fit1), regexp = "range: 33.23")

  # Anisotropy when not shared across random fields
  set.seed(1)
  fit2 <- sdmTMB(
    data = pcod_2011,
    formula = density ~ 1,
    mesh = pcod_mesh_2011,
    family = tweedie(),
    share_range = FALSE,
    time = "year",
    anisotropy = TRUE
  )

  expect_output(cat(print_anisotropy(fit2)), regexp = "\\(spatial\\): 20.11 to 65.79 at 117")
  expect_output(cat(print_anisotropy(fit2)), regexp = "\\(spatiotemporal\\): 0.01 to 0.02 at 117")

  # Anisotropy when shared across random fields in delta model
  fit_dg_shared <- sdmTMB(
    data = pcod_2011,
    formula = density ~ 1,
    mesh = pcod_mesh_2011,
    family = delta_gamma(),
    share_range = TRUE,
    time = "year",
    anisotropy = TRUE,
    control = sdmTMBcontrol(newton_loops = 1L)
  )

  expect_output(cat(print_anisotropy(fit_dg_shared, m = 1)), regexp = "\\(spatial\\): 36.42 to 80.86 at 34")
  expect_output(cat(print_anisotropy(fit_dg_shared, m = 2)), regexp = "\\(spatial\\): 2.09 to 4.65 at 34")

  # Anisotropy when not shared across random fields in delta model
  test_mesh <- make_mesh(data = pcod, xy_cols = c("X", "Y"), cutoff = 20)

  fit_dg_not_shared <- sdmTMB(
    data = pcod,
    formula = density ~ 1,
    mesh = test_mesh,
    family = delta_gamma(),
    share_range = FALSE,
    time = "year",
    anisotropy = TRUE,
    control = sdmTMBcontrol(newton_loops = 1L)
  )

  expect_output(cat(print_anisotropy(fit_dg_not_shared, m = 1)), regexp = "\\(spatial\\): 19.46 to 64.60 at 59")
  expect_output(cat(print_anisotropy(fit_dg_not_shared, m = 1)), regexp = "\\(spatiotemporal\\): 122.96 to 408.11 at 59")
  expect_output(cat(print_anisotropy(fit_dg_not_shared, m = 2)), regexp = "\\(spatial\\): 0.01 to 0.04 at 59")
  expect_output(cat(print_anisotropy(fit_dg_not_shared, m = 2)), regexp = "\\(spatiotemporal\\): 9.10 to 30.22 at 59")
})
