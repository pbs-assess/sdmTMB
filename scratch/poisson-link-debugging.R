library(VAST)
library(sdmTMB)
library(units)
library(dplyr)
library(ggplot2)

sp <- "arrowtooth-flounder"

dat <- readRDS(paste0("../gfsynopsis-2021/report/data-cache-feb-2023/", sp, ".rds"))$survey_sets
dat <- as.data.frame(dat) |>
  dplyr::filter(survey_abbrev %in% "SYN QCS")
dat$area_swept1 <- dat$doorspread_m * (dat$speed_mpm * dat$duration_min)
dat$area_swept2 <- dat$tow_length_m * dat$doorspread_m
dat$area_swept <- ifelse(!is.na(dat$area_swept2), dat$area_swept2, dat$area_swept1)

# VAST ----------------------------------------------------------------

FieldConfig <- matrix(c("IID", "IID", "IID", "IID", "IID", "IID"),
  ncol = 2, nrow = 3,
  dimnames = list(
    c("Omega", "Epsilon", "Beta"),
    c("Component_1", "Component_2")
  )
)
FieldConfig
RhoConfig <- c("Beta1" = 0, "Beta2" = 0, "Epsilon1" = 0, "Epsilon2" = 0)

# get coordinates in geographic coordinates from UTM projection,
# after changing units to m
qcs_grid_ll <- replicate_df(sdmTMB::qcs_grid, "year", unique(dat$year))
qcs_grid_ll$Y <- qcs_grid_ll$Y * 1000
qcs_grid_ll$X <- qcs_grid_ll$X * 1000

qcs_grid_ll <- sf::st_as_sf(
  x = qcs_grid_ll,
  coords = c("X", "Y"),
  crs = "+proj=utm +zone=9"
)
qcs_grid_ll <- sf::st_transform(qcs_grid_ll, crs = "+proj=longlat +datum=WGS84")
co <- sf::st_coordinates(qcs_grid_ll)
qcs_grid_ll$X <- as.numeric(co[,1])
qcs_grid_ll$Y <- as.numeric(co[,2])

# remove replicate locations for each year and format for VAST
qcs_grid_ll <- subset(qcs_grid_ll, year == min(qcs_grid_ll$year))
input_grid <- cbind(Lat = qcs_grid_ll$Y, Lon = qcs_grid_ll$X, Area_km2 = 4)

settings <- make_settings(
  n_x = 135, # number of vertices in the SPDE mesh
  Region = "User",
  purpose = "index2", # index of abundance with Gamma for positive catches
  fine_scale = TRUE, # use bilinear interpolation from the INLA 'A' matrix
  zone = 9,
  FieldConfig = FieldConfig,
  RhoConfig = RhoConfig,
  ObsModel = c(2, 1), # conventional logit-linked delta-Gamma; c(10, 2) for Tweedie
  bias.correct = FALSE,
  use_anisotropy = FALSE,
  max_cells = Inf, # use all grid cells from the extrapolation grid
  knot_method = "grid" # or "samples"
)

input_grid <- cbind(Lat = qcs_grid_ll$Y, Lon = qcs_grid_ll$X, Area_km2 = 4)

b_i <- dat$catch_weight
units(b_i) <- "kg"
a_i <- dat$area_swept/10000
# a_i <- rep(1.0, nrow(dat))
units(a_i) <- "km^2"

dir.create("scratch/vast-test", showWarnings = FALSE)
tictoc::tic()
fv <- fit_model(
  settings = settings,
  Lat_i = dat[["latitude"]],
  Lon_i = dat[["longitude"]],
  t_i = dat[["year"]],
  b_i = b_i,
  a_i = a_i,
  input_grid = input_grid,
  working_dir = paste0(here::here("scratch", "vast-test"), "/")
)
tictoc::toc()

fv$parameter_estimates$SD

plot(
  fv,
  check_residuals = FALSE,
  working_dir = paste0(here::here("scratch", "vast-test"), "/")
)

# sdmTMB --------------------------------------------------------------

dat <- sdmTMB::add_utm_columns(dat)
mesh <- make_mesh(dat, xy_cols = c("X", "Y"),
  mesh = fv$spatial_list$MeshList$isotropic_mesh)

dat$offset <- log(dat$area_swept/10000)
# dat$offset <- 0
mean(dat$offset)

nd <- replicate_df(sdmTMB::qcs_grid, "year", unique(dat$year))

tictoc::tic()
fs <- sdmTMB(
  catch_weight ~ 0 + as.factor(year),
  data = dat, mesh = mesh,
  offset = "offset",
  family = delta_poisson_link_gamma(),
  time = "year", silent = FALSE,
  # do_fit = F
  do_index = TRUE,
  index_args = list(area = rep(4, nrow(nd)), bias_correct = TRUE),
  predict_args = list(newdata = nd),
  control = sdmTMBcontrol(
    # start = list(
    #   ln_tau_O = c(0.9191302, 0.5868656),
    #   ln_tau_E = c(2.2391263, 0.6580646),
    #   ln_kappa = rbind(
    #     c(-2.4824288, -1.9137340),
    #     c(-2.4824288, -1.9137340)
    #   ),
    #   ln_phi = c(0, 0.3032190)
    # ),
  newton_loops = 1L
  )
)
tictoc::toc()

tictoc::tic()
ind <- get_index(fs, bias_correct = FALSE)
tictoc::toc()

vast_i <- read.csv("scratch/vast-test/Index.csv") %>%
  mutate(index = "VAST", year = as.numeric(Time), est = Estimate,
    se = Std..Error.for.ln.Estimate.) %>%
  select(index, year, est, se) %>%
  mutate(lwr = exp(log(est) + qnorm(0.025) * se)) %>%
  mutate(upr = exp(log(est) + qnorm(0.975) * se))
sdm_i <- ind %>% mutate(index = "sdmTMB")

both_i <- bind_rows(sdm_i, vast_i) %>% filter(est > 0)
ggplot(both_i, aes(x = year, y = est, ymin = lwr, ymax = upr, colour = index)) +
  geom_ribbon(alpha = 0.1) +
  geom_line(alpha = 0.8) +
  ylim(0, max(both_i$upr)) +
  coord_cartesian(expand = FALSE) +
  facet_wrap(~index)

# tictoc::tic()
# fs2 <- sdmTMB(
#   catch_weight ~ 0 + as.factor(year),
#   data = dat, mesh = mesh,
#   offset = "offset",
#   family = delta_gamma(),
#   time = "year", silent = FALSE
# )
# tictoc::toc()

fs
sanity(fs)
fs$sd_report

plot_betas <- function(vast_model, sdmTMB_model, vast_par = "beta1_ft", sdmTMB_pars = 1) {
  s <- vast_model$parameter_estimates$SD
  vast_est1 <- as.list(s, "Estimate", report = FALSE)
  vast_est2 <- as.list(s, "Estimate", report = TRUE)
  vast_sd1 <- as.list(s, "Std. Error", report = FALSE)
  vast_sd2 <- as.list(s, "Std. Error", report = TRUE)
  sdmTMB_est <- as.list(sdmTMB_model$sd_report, "Estimate", report = FALSE)
  sdmTMB_sd <- as.list(sdmTMB_model$sd_report, "Std. Error", report = FALSE)
  b_year_vast <- vast_est1[[vast_par]][!is.na(vast_sd1[[vast_par]])]
  b_year_vast_se <- vast_sd1[[vast_par]][!is.na(vast_sd1[[vast_par]])]
  years <- sort(unique(dat$year))
  lwr_vast <- b_year_vast - 2 * b_year_vast_se
  upr_vast <- b_year_vast + 2 * b_year_vast_se
  plot(years, b_year_vast, ylim = range(c(lwr_vast, upr_vast)))
  segments(years, lwr_vast, years, upr_vast)
  years <- years + 0.05
  if (sdmTMB_pars == 1) {
    points(years, sdmTMB_est$b_j)
    segments(years, sdmTMB_est$b_j - 2 * sdmTMB_sd$b_j,
      years, sdmTMB_est$b_j + 2 * sdmTMB_sd$b_j,
      col = "red")
  } else {
    points(years, sdmTMB_est$b_j2)
    segments(years, sdmTMB_est$b_j2 - 2 * sdmTMB_sd$b_j2,
      years, sdmTMB_est$b_j2 + 2 * sdmTMB_sd$b_j2,
      col = "red")
  }
  legend("topright", legend = c("VAST", "sdmTMB"),
    col = c("black", "red"), bty = "n", lty = c(1, 1))
}

par(mfrow = c(2, 1), cex = 0.8, mar = c(1.5, 1, 1, 1), oma = c(2, 3, 1, 1))
plot_betas(fv, fs, "beta1_ft", sdmTMB_pars = 1)
plot_betas(fv, fs, "beta2_ft", sdmTMB_pars = 2)

