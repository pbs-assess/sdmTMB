test_that("r2 function works", {
  .compare <- function(.p1, .p2, tol = 1e-3) {
    expect_equal(.p1$R2[p1$component == "conditional"], unname(.p2$R2_conditional), tolerance = tol)
    expect_equal(.p1$R2[p1$component == "marginal"], unname(.p2$R2_marginal), tolerance = tol)
  }

  d <- dogfish
  d$fyear <- as.factor(dogfish$year)
  d$density <- d$catch_weight / d$area_swept
  d$scaled_log_depth <- as.numeric(scale(log(d$depth)))
  dpos <- d[d$density > 0, ]

  # gaussian
  fit <- sdmTMB(log(density) ~ scaled_log_depth + (1 | fyear), dpos, spatial = "off")
  fit2 <- glmmTMB::glmmTMB(log(density) ~ scaled_log_depth + (1 | fyear), dpos)
  p1 <- r2(fit)
  p2 <- performance::r2_nakagawa(fit2)
  .compare(p1, p2)

  # tweedie
  fit <- sdmTMB(density ~ scaled_log_depth + (1 | fyear), data = d, spatial = "off", family = tweedie())
  fit2 <- glmmTMB::glmmTMB(density ~ scaled_log_depth + (1 | fyear), data = d, family = glmmTMB::tweedie())
  p1 <- r2(fit)
  p2 <- performance::r2_nakagawa(fit2)
  .compare(p1, p2)

  # multiple random intercepts:
  set.seed(1)
  d$g <- as.factor(sample(letters[1:5], size = nrow(d), replace = TRUE))
  fit <- sdmTMB(density ~ scaled_log_depth + (1 | g) + (1 | fyear), data = d, spatial = "off", family = tweedie())
  fit2 <- glmmTMB::glmmTMB(density ~ scaled_log_depth + (1 | g) + (1 | fyear), data = d, family = glmmTMB::tweedie())
  p1 <- r2(fit)
  p2 <- r2(fit2)
  .compare(p1, p2)

  # no (1 | g)
  fit <- sdmTMB(
    density ~ s(depth),
    data = pcod_2011, mesh = pcod_mesh_2011,
    family = tweedie(link = "log")
  )
  p1 <- r2(fit)
  p1
  expect_equal(p1$R2[p1$component == "conditional"], expected = 0.9142605, tolerance = 1e-3)

  # binomial
  fit <- sdmTMB(present ~ scaled_log_depth + (1 | fyear), data = d, spatial = "off", family = binomial())
  fit2 <- glmmTMB::glmmTMB(present ~ scaled_log_depth + (1 | fyear), data = d, family = binomial())
  p1 <- r2(fit)
  p2 <- performance::r2_nakagawa(fit2)
  .compare(p1, p2)

  # gamma
  fit <- sdmTMB(density ~ scaled_log_depth + (1 | fyear), data = dpos, spatial = "off", family = Gamma(link = "log"))
  fit2 <- glmmTMB::glmmTMB(density ~ scaled_log_depth + (1 | fyear), data = dpos, family = Gamma(link = "log"))
  p1 <- r2(fit)
  p2 <- performance::r2_nakagawa(fit2)
  .compare(p1, p2)

  # poisson
  d$count <- round(d$density)
  fit <- sdmTMB(count ~ scaled_log_depth + (1 | fyear), data = d, spatial = "off", family = poisson(link = "log"))
  fit2 <- glmmTMB::glmmTMB(count ~ scaled_log_depth + (1 | fyear), data = d, family = poisson(link = "log"))
  p1 <- r2.sdmTMB(fit)
  p2 <- performance::r2_nakagawa(fit2)
  .compare(p1, p2)

  # spatial
  mesh_pos <- make_mesh(dpos, c("X", "Y"), cutoff = 15)
  fit <- sdmTMB(log(density) ~ scaled_log_depth + (1 | fyear), dpos, spatial = "on", mesh = mesh_pos)
  r2(fit)
  expect_true("partial_spatial" %in% r2(fit)$component)

  # ST
  fit <- sdmTMB(log(density) ~ scaled_log_depth + (1 | fyear), dpos, time = "year", mesh = mesh_pos)
  r2(fit)
  expect_true("partial_spatiotemporal" %in% r2(fit)$component)

  # ST + smooths
  fit <- sdmTMB(log(density) ~ s(scaled_log_depth) + (1 | fyear), dpos, time = "year", mesh = mesh_pos)
  r2(fit)
  expect_true("partial_smoothers" %in% r2(fit)$component)
})
