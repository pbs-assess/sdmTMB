test_that("Print anisotropy prints correctly", {
  # skip_on_cran()
  skip_if_not_installed("INLA")

  # No anisotropy
  fit1 <- sdmTMB(
    data = pcod_2011,
    formula = density ~ 1,
    mesh = pcod_mesh_2011,
    family = tweedie(),
    anisotropy = FALSE
  )

  expect_output(print(fit1), regexp = "Matern range: 33.23")

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
  a_df <- plot_anisotropy(fit2, return_data = TRUE)
  a_df$degree <- a_df$angle * 180 / pi

  expect_output(cat(print_anisotropy(fit2)), regexp = "\\(spatial\\): 20.11 to 65.78 at 116.98")
  expect_output(cat(print_anisotropy(fit2)), regexp = "\\(spatiotemporal\\): 0.01 to 0.04 at 116.98")

  # Anisotropy when shared across random fields in delta model
  fit_dg_shared <- sdmTMB(
    data = pcod_2011,
    formula = density ~ 1,
    mesh = pcod_mesh_2011,
    family = delta_gamma(),
    share_range = TRUE,
    time = "year",
    anisotropy = TRUE
  )
  a_df <- plot_anisotropy(fit_dg_shared, return_data = TRUE)
  a_df$degree <- a_df$angle * 180 / pi

  expect_output(cat(print_anisotropy(fit_dg_shared, m = 1)), regexp = "\\(spatial\\): 36.42 to 80.86 at 33.52")
  expect_output(cat(print_anisotropy(fit_dg_shared, m = 2)), regexp = "\\(spatial\\): 36.42 to 80.86 at 33.52")

  # Anisotropy when not shared across random fields in delta model
  pcod_test <- pcod
  test_mesh <- make_mesh(data = pcod, xy_cols = c("X", "Y"), cutoff = 20)

  fit_dg_shared <- sdmTMB(
    data = pcod_test,
    formula = density ~ 1,
    mesh = test_mesh,
    family = delta_gamma(),
    share_range = FALSE,
    time = "year",
    anisotropy = TRUE
  )
  a_df <- plot_anisotropy(fit_dg_shared, return_data = TRUE)
  a_df$degree <- a_df$angle * 180 / pi

  expect_output(cat(print_anisotropy(fit_dg_not_shared, m = 1)), regexp = "\\(spatial\\): 19.47 to 64.60 at 58.72")
  expect_output(cat(print_anisotropy(fit_dg_not_shared, m = 1)), regexp = "\\(spatiotemporal\\): 122.96 to 408.09 at 58.72")
  expect_output(cat(print_anisotropy(fit_dg_not_shared, m = 2)), regexp = "\\(spatial\\): 0.02 to 0.06 at 58.72")
  expect_output(cat(print_anisotropy(fit_dg_not_shared, m = 2)), regexp = "\\(spatiotemporal\\): 9.10 to 30.22 at 58.72")
})
