test_that("Dispersion formula works", {
  # 1 phi
  m1 <- sdmTMB(
    density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
    data = pcod_2011, time = "year", mesh = pcod_mesh_2011,
    family = tweedie(link = "log"), spatial = "off", spatiotemporal = "off"
  )

  # 1 phi delta
  m1d <- sdmTMB(
    density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
    data = pcod_2011, time = "year", mesh = pcod_mesh_2011,
    family = delta_gamma(), spatial = "off", spatiotemporal = "off"
  )

  # phi by year
  m2 <- sdmTMB(
    density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
    data = pcod_2011, time = "year", mesh = pcod_mesh_2011, family = tweedie(link = "log"),
    dispformula = ~ 0 + as.factor(year),
    spatial = "off", spatiotemporal = "off"
  )
  # phi by year delta
  m2d <- sdmTMB(
    density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
    data = pcod_2011, time = "year", mesh = pcod_mesh_2011, family = delta_gamma(),
    dispformula = ~ 0 + as.factor(year),
    spatial = "off", spatiotemporal = "off"
  )
  # 1 phi
  m3 <- glmmTMB::glmmTMB(
    density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
    data = pcod_2011, family = glmmTMB::tweedie(link = "log"),
    dispformula = ~ 1
  )
  # phi by year
  m4 <- glmmTMB::glmmTMB(
    density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
    data = pcod_2011, family = glmmTMB::tweedie(link = "log"),
    dispformula = ~ 0 + as.factor(year)
  )

  x1 <- summary(m1$sd_report)
  p1 <- as.numeric(x1[row.names(x1) == "b_disp_k", 1])
  x2 <- summary(m2$sd_report)
  p2 <- as.numeric(x2[row.names(x2) == "b_disp_k", 1])
  x3 <- summary(m3$sdr)
  p3 <- as.numeric(x3[row.names(x3) == "betad", 1])
  x4 <- summary(m4$sdr)
  p4 <- as.numeric(x4[row.names(x4) == "betad", 1])

  expect_equal(p1, p3, tolerance = 1e-5)
  expect_equal(p2, p4, tolerance = 1e-5)

  expect_output(print(m2), "Dispersion model")
  expect_output(print(m2d), "Dispersion model")

  expect_warning({
    m <- sdmTMB(
      present ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
      data = pcod_2011, time = "year", mesh = pcod_mesh_2011, family = binomial(),
      dispformula = ~ 0 + as.factor(year),
      spatial = "off", spatiotemporal = "off"
    )
  }, regexp = "ignored")
})
