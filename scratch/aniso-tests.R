mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 6)
fit <- sdmTMB(
  density ~ 1,
  data = pcod, mesh = mesh,
  anisotropy = TRUE,
  control = sdmTMBcontrol(map = list(ln_H_input = factor(c(1, 2, 3, 4)))),
  family = delta_gamma()
)
plot_anisotropy2(fit)
plot_anisotropy2(fit, return_data = TRUE)

fit <- sdmTMB(
  density ~ 1,
  data = pcod, mesh = mesh,
  anisotropy = TRUE,
  control = sdmTMBcontrol(map = list(ln_H_input = factor(c(1, 2, 1, 2)))),
  family = delta_gamma()
)
plot_anisotropy2(fit)
plot_anisotropy2(fit, return_data = TRUE)

fit <- sdmTMB(
  density ~ 1,
  data = pcod, mesh = mesh,
  anisotropy = TRUE,
  family = tweedie()
)
plot_anisotropy2(fit)
plot_anisotropy2(fit, return_data = TRUE)

mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 9)
fit <- sdmTMB(
  density ~ 1,
  data = pcod, mesh = mesh, spatiotemporal = "iid",
  time = "year", silent = FALSE,
  anisotropy = TRUE, share_range = FALSE,
  family = tweedie()
)
plot_anisotropy2(fit)
plot_anisotropy2(fit, return_data = TRUE)

mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 9)
fit <- sdmTMB(
  density ~ 1,
  data = pcod, mesh = mesh,
  anisotropy = TRUE,
  time = "year",
  share_range = FALSE,
  family = tweedie()
)
plot_anisotropy2(fit)
plot_anisotropy2(fit, return_data = TRUE)

mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 8)
fit <- sdmTMB(
  density ~ as.factor(year),
  data = pcod,
  mesh = mesh,
  anisotropy = TRUE,
  time = "year",
  share_range = FALSE,
  silent = FALSE,
  family = delta_gamma()
)
plot_anisotropy2(fit)
plot_anisotropy2(fit, return_data = TRUE)

mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 10)
fit <- sdmTMB(
  density ~ as.factor(year),
  data = pcod,
  mesh = mesh,
  anisotropy = TRUE,
  time = "year",
  spatial = "off",
  silent = FALSE,
  family = delta_gamma()
)
plot_anisotropy2(fit)
plot_anisotropy2(fit, return_data = TRUE)

fit <- sdmTMB(
  density ~ as.factor(year),
  data = pcod,
  mesh = mesh,
  anisotropy = TRUE,
  time = "year",
  spatiotemporal = "off",
  silent = FALSE,
  family = delta_gamma()
)
plot_anisotropy2(fit)
plot_anisotropy2(fit, return_data = TRUE)

fit <- sdmTMB(
  density ~ as.factor(year),
  data = pcod,
  mesh = mesh,
  anisotropy = TRUE,
  time = "year",
  share_range = FALSE,
  spatiotemporal = list("off", "iid"),
  silent = FALSE,
  family = delta_gamma()
)
plot_anisotropy2(fit)
plot_anisotropy2(fit, return_data = TRUE)
