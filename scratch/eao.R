library(TMB)
library(VAST)
library(sp)
library(sdmTMB)
library(here)

FieldConfig <- matrix(c("0", "0", "IID", "Identity", "IID", "IID", "IID", "Identity"),
  ncol = 2, nrow = 4,
  dimnames = list(
    c("Omega", "Epsilon", "Beta", "Epsilon_year"),
    c("Component_1", "Component_2")
  )
)
FieldConfig

RhoConfig <- c("Beta1" = 3, "Beta2" = 0, "Epsilon1" = 0, "Epsilon2" = 0)

qcs_grid_ll <- replicate_df(qcs_grid, "year", unique(pcod$year))
qcs_grid_ll$Y <- qcs_grid_ll$Y * 1000
qcs_grid_ll$X <- qcs_grid_ll$X * 1000

coordinates(qcs_grid_ll) <- ~ X + Y
proj4string(qcs_grid_ll) <- CRS("+proj=utm +zone=9")
qcs_grid_ll <- as.data.frame(spTransform(qcs_grid_ll, CRS("+proj=longlat +datum=WGS84")))

qcs_grid_ll <- subset(qcs_grid_ll, year == min(qcs_grid_ll$year))
input_grid <- cbind(Lat = qcs_grid_ll$coords.x2, Lon = qcs_grid_ll$coords.x1, Area_km2 = 4)

settings <- make_settings(
  n_x = 205, # number of vertices in the SPDE mesh
  Region = "User",
  purpose = "index2", # use recommended defaults for an index of abundance fine_scale = TRUE, # use bilinear interpolation from the INLA 'A' matrix zone = 9,
  FieldConfig = FieldConfig,
  RhoConfig = RhoConfig,
  ObsModel = c(10, 2), # use the Tweedie distribution as the observation model bias.correct = FALSE,
  use_anisotropy = FALSE,
  Options = list(Calculate_effective_area = TRUE),
  max_cells = Inf, # use all grid cells from the extrapolation grid knot_method = "grid" # or "samples"
)

dir.create(paste0(here("scratch", "VAST-EAO-test")), showWarnings = FALSE)

pcod$effort <- 1
pcod <- as.data.frame(pcod) # ensure not a tibble
fit <- fit_model(
  settings = settings,
  Lat_i = pcod[, "lat"],
  Lon_i = pcod[, "lon"],
  t_i = pcod[, "year"],
  b_i = pcod[, "density"],
  a_i = pcod[, "effort"],
  input_grid = input_grid,
  working_dir = paste0(here("scratch", "VAST-EAO-test"), "/")
)

fit$parameter_estimates$SD

mesh <- sdmTMB::make_mesh(pcod,
  xy_cols = c("X", "Y"),
  mesh = fit$spatial_list$MeshList$isotropic_mesh
)

fit_sdmTMB <- sdmTMB(density ~ 0 + as.factor(year),
  data = pcod, mesh = mesh, family = tweedie(link = "log"),
  time = "year", spatiotemporal = "iid", spatial = "on"
)
nd <- replicate_df(qcs_grid, "year", unique(pcod$year))
p <- predict(fit_sdmTMB, newdata = nd, return_tmb_object = TRUE)
ind <- get_index(p, area = 4, bias_correct = FALSE)

s1 <- fit$ParHat$beta1_ft[fit$ParHat$beta2_ft != 0]
s2 <- fit$ParHat$beta2_ft[fit$ParHat$beta2_ft != 0]
b_VAST <- as.numeric(s1 + s2)
b_sdmTMB <- tidy(fit_sdmTMB)
plot(b_sdmTMB$estimate, b_VAST)
abline(0, 1)

s <- fit$parameter_estimates$SD
b <- as.list(s, what = "Estimate", report = TRUE)

b$mean_D_ctl

b$effective_area_ctl
b$log_effective_area_ctl

b.se <- as.list(s, what = "Std. Error", report = TRUE)
b.se$log_effective_area_ctl

vast_est1 <- as.list(s, "Estimate", report = FALSE)
vast_est2 <- as.list(s, "Estimate", report = TRUE)
vast_sd1 <- as.list(s, "Std. Error", report = FALSE)
vast_sd2 <- as.list(s, "Std. Error", report = TRUE)

exp(vast_est2$ln_Index_ctl)
ind$est

# now try in sdmTMB...

library(dplyr)
group_by(p$data, year) |>
  summarise(total = sum(exp(est) * 4))

b$mean_D_ctl

group_by(p$data, year) |>
  summarise(mean_dens = weighted.mean(exp(est), w = exp(est)))

group_by(p$data, year) |>
  summarise(mean_dens = weighted.mean(exp(est), w = exp(est)))

# supp. eqn. version:
group_by(p$data, year) |>
  mutate(d_jt = exp(est), b_t = sum(d_jt * 4), a_j = 4) |>
  summarise(mean_dens = sum((a_j * d_jt)/b_t[1] * d_jt))

group_by(p$data, year) |>
  mutate(d_jt = exp(est), b_t = sum(d_jt), a_j = 4) |>
  summarise(
    weighted_mean = weighted.mean(exp(est), w = exp(est)),
    mean1 = mean(d_jt),
    mean2 = sum(d_jt * d_jt/sum(d_jt)) #< simplest
  )

exp(b$log_effective_area_ctl)
exp(b$ln_Index_ctl)
exp(b$log_mean_D_ctl)

s <- fit$parameter_estimates$SD
b <- as.list(s, what = "Estimate", report = TRUE)

b$mean_D_ctl

b$effective_area_ctl
b$log_effective_area_ctl

exp(vast_est2$ln_Index_ctl)
ind$est

# now try in sdmTMB...

group_by(p$data, year) |>
  summarise(total = sum(exp(est) * 4))

b$mean_D_ctl

group_by(p$data, year) |>
  summarise(mean_dens = weighted.mean(exp(est), w = exp(est)))

group_by(p$data, year) |>
  summarise(mean_dens = weighted.mean(exp(est), w = exp(est)))

# supp. eqn. version:
group_by(p$data, year) |>
  mutate(d_jt = exp(est), b_t = sum(d_jt * 4), a_j = 4) |>
  summarise(mean_dens = sum((a_j * d_jt)/b_t[1] * d_jt))

group_by(p$data, year) |>
  mutate(d_jt = exp(est), b_t = sum(d_jt), a_j = 4) |>
  summarise(
    weighted_mean = weighted.mean(exp(est), w = exp(est)),
    mean1 = mean(d_jt),
    mean2 = sum(d_jt * d_jt/sum(d_jt)) #< simplest
  )

exp(b$log_effective_area_ctl)
exp(b$ln_Index_ctl)
exp(b$log_mean_D_ctl)

group_by(p$data, year) |>
  mutate(d_jt = exp(est), b_t = sum(d_jt * 4), a_j = 4) |>
  summarise(
    b_t = b_t[1],
    weighted_mean = weighted.mean(exp(est), w = exp(est)),
    mean1 = mean(d_jt),
    mean2 = sum(d_jt * d_jt/sum(d_jt)),
    h_t_supp = b_t[1]^2 / sum(a_j * d_jt^2),
    h_t = b_t[1] / mean2 #< simplest
  )

eao0 <- get_eao(p, area = 4, bias_correct = FALSE)
eao0

no_data <- vast_est1$beta2_ft[1,] == 0
ve <- b$log_effective_area_ctl[,,1]

plot(eao0$log_est, ve[!no_data]);abline(0, 1)

eao <- get_eao(p, area = 4, bias_correct = TRUE)
eao

