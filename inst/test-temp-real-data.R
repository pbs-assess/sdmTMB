pcod_spde <- make_spde(pcod$X, pcod$Y, n_knots = 100)
m <- sdmTMB(
  pcod, density ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
  time = "year", spde = pcod_spde, family = tweedie(link = "log"),
  silent = FALSE
)

# Predictions onto new data:
# expand for time units:
newdata <- qcs_grid
original_time <- sort(unique(m$data[[m$time]]))
nd <- do.call("rbind",
  replicate(length(original_time), newdata, simplify = FALSE))
nd[[m$time]] <- rep(original_time, each = nrow(newdata))

p <- predict(m, newdata = nd)
cog <- get_cog(p)
cog

library(ggplot2)
ggplot(cog, aes(year, est, ymin = lwr, ymax = upr)) +
  geom_ribbon(alpha = 0.5) + geom_line() + facet_wrap(~coord, scales = "free_y")


# get_index(p, value_name = "link_total", bias_correct = FALSE)

library(dplyr)
data.frame(X = p$data$X, est = exp(p$data$est), year = p$data$year) %>%
  group_by(year) %>% summarize(cog = sum(X * est) / sum(est))

predictions <- p$data
# A short function for plotting our predictions:
plot_map <- function(dat, column = "est") {
  ggplot(dat, aes_string("X", "Y", fill = column)) +
    geom_raster() +
    facet_wrap(~year) +
    coord_fixed()
}

plot_map(predictions, "exp(est)") +
  scale_fill_viridis_c(trans = "sqrt") +
  ggtitle("Prediction (fixed effects + all random effects)")


# Simulation
library(ggplot2)
set.seed(183228)

# simulate data
dat <- sim(
  time_steps = 9, ar1_fields = TRUE, ar1_phi = 0.5,
  sigma_O = 0.001, sigma_E = 0.3, phi = 0.05
)


spde <- make_spde(x = dat$x, y = dat$y, n_knots = 200)

m <- sdmTMB(silent = FALSE, ar1_fields = TRUE, include_spatial = FALSE,
  data = dat, formula = z ~ 1, time = "time",
  family = gaussian(link = "identity"), spde = spde
)


# create fine-scale grid to predict on
newdata <- expand.grid(x = seq(min(dat$x), max(dat$x), length.out = 100),
  y = seq(min(dat$y), max(dat$y), length.out = 100))

# replicate for each time period
original_time <- sort(unique(m$data[[m$time]]))
nd <- do.call("rbind",
  replicate(length(original_time), newdata, simplify = FALSE))
nd[[m$time]] <- rep(original_time, each = nrow(newdata))

# run TMB with prediction turned on
predictions <- predict(m, newdata = nd, xy_cols = c("x", "y"))$data


plot_map <- function(dat, fill = "est") {
  ggplot(dat, aes_string("x", "y", fill = fill)) +
    geom_raster() +
    facet_wrap(~time) +
    coord_fixed()
}

plot_map(predictions, "est") +
  scale_fill_viridis_c() +
  ggtitle("Prediction (fixed effects + all random effects)")

plot_map(predictions, "est") +
  scale_fill_gradient2() +
  ggtitle("Prediction (fixed effects + all random effects)")

plot_map(predictions, "est_re_st") +
  ggtitle("Spatiotemporal random effects only") +
  scale_fill_gradient2()



set.seed(1838)
dat <- sim(
  time_steps = 9, ar1_fields = TRUE, ar1_phi = 0.5,
  sigma_O = 0.01, sigma_E = 0.3, phi = 0.01
)

# create subset to fit to


spde <- make_spde(x = dat$x, y = dat$y, n_knots = 120)
m <- sdmTMB(silent = FALSE, ar1_fields = TRUE, include_spatial = FALSE,
  data = dat, formula = z ~ 1, time = "time",
  family = gaussian(link = "identity"), spde = spde
)

# create fine-scale grid to predict on
newdata <- expand.grid(x = seq(min(dat$x), max(dat$x), length.out = 100),
  y = seq(min(dat$y), max(dat$y), length.out = 100))

# replicate for each time period
original_time <- sort(unique(m$data[[m$time]]))
nd <- do.call("rbind",
  replicate(length(original_time), newdata, simplify = FALSE))
nd[[m$time]] <- rep(original_time, each = nrow(newdata))

# run TMB with prediction turned on
predictions <- predict(m, newdata = nd, xy_cols = c("x", "y"))$data

# contrast predictions with simulated true data


r <- m$tmb_obj$report()
r$sigma_E # spatio-temporal
exp(r$ln_kappa) # Matern for both space and space-time
2 * plogis(m$model$par[['ar1_phi']]) - 1 # back transform

# -----------------------------------------

library(dplyr)
survey <- "SYN QCS"
d <- gfplot::get_survey_sets("pacific cod", ssid = 1)

d <- readRDS("~/src/gfsynopsis/report/data-cache/pacific-ocean-perch.rds")
d <- readRDS("~/src/gfsynopsis/report/data-cache/quillback-rockfish.rds")
d <- readRDS("~/src/gfsynopsis/report/data-cache/rougheye-blackspotted-rockfish-complex.rds")
d <- readRDS("~/src/gfsynopsis/report/data-cache/petrale-sole.rds")
d <- readRDS("~/src/gfsynopsis/report/data-cache/pacific-cod.rds")

d <- d$survey_sets
sd <- gfplot::get_sensor_data_trawl(ssid = 1, attribute = 1)

# d2 <- left_join(d, sd)
# d2 <- mutate(d2, depth_m = avg_value) # Temporary hack
col <- if (grepl("SYN", survey)) "density_kgpm2" else "density_ppkm2"
dat <- gfplot:::tidy_survey_sets(d, survey, years = seq(1, 1e6),
  density_column = col)

dat <- left_join(dat, select(sd, avg_value, year, fishing_event_id)) %>%
  rename(temperature = avg_value)
# if (any(is.na(dat$depth)))
# dat <- gfplot:::interp_survey_bathymetry(dat)$data
nrow(dat)
dat <- filter(dat, !is.na(temperature), !is.na(depth))
nrow(dat)
# dat <- gfplot:::scale_survey_predictors(dat)

spde <- sdmTMB::make_spde(dat$X, dat$Y, n_knots = 150)
plot_spde(spde)

dat <- mutate(dat,
  temperature_mean = mean(temperature),
  temperature_sd = sd(temperature),
  depth_mean = mean(depth),
  depth_sd = sd(depth),
  depth_scaled = (depth - mean(depth)) / sd(depth),
  temperature_scaled = (temperature - mean(temperature)) / sd(temperature),
  temperature_scaled2 = temperature_scaled^2,
  depth_scaled2 = depth_scaled^2,
)

m <- sdmTMB::sdmTMB(
  formula = density ~ 0 + as.factor(year),# + temperature_scaled + temperature_scaled2,
  time_varying =
    ~ 0 +
    # depth_scaled + depth_scaled2,
    temperature_scaled + temperature_scaled2,
  data = dat, time = "year",
  ar1_fields = FALSE,
  spde = spde, include_spatial = TRUE,
  family = sdmTMB::tweedie(link = "log"),
  anisotropy = FALSE,
  silent = FALSE)
# 2 * plogis(m$model$par[['ar1_phi']]) - 1 # back transform

get_y_hat <- function(b0, b1, b2, year, predictor, mean_column, sd_column) {
  x_pred <- seq(min(dat[[predictor]]), max(dat[[predictor]]), length.out = 300)
  data.frame(
    x = x_pred * dat[[sd_column]][[1]] + dat[[mean_column]][[1]],
    # x = -exp((x_pred * dat$depth_sd[1] + dat$depth_mean[1])),
    y_hat = exp(b0 + b1 * x_pred + b2 * x_pred^2),
    year = year)
}

# r <- m$tmb_obj$report()
# r$b_rw_t
# n_t <- nrow(r$b_rw_t)
# yrs <- sort(unique(dat$year))
# pred_depth <- purrr::map_df(seq_len(n_t), function(.t) {
#   get_y_hat(r$b_j[.t], r$b_rw_t[.t,1], r$b_rw_t[.t,2], yrs[.t], sd_column = "depth_sd",
#     mean_column = "depth_mean", predictor = "depth_scaled")
# })
#
# library(ggplot2)
# ggplot(pred_depth, aes(x = x, ymax = 10000*y_hat + year*0.0, ymin = 0 + year*0.0,
#   group = year, colour = year, fill = year)) +
#   geom_ribbon(lwd = 0.5, alpha = 0.2) +
#   xlab("Depth") +
#   scale_color_viridis_c(option = "C") +
#   scale_fill_viridis_c(option = "C") +
#   ylab("Predicted density in some units") +
#   gfplot::theme_pbs() +
#   coord_cartesian(expand = FALSE)

pred_temperature <- purrr::map_df(seq_len(n_t), function(.t) {
  get_y_hat(r$b_j[.t], r$b_rw_t[.t,1], r$b_rw_t[.t,2], yrs[.t], sd_column = "temperature_sd",
    mean_column = "temperature_mean", predictor = "temperature_scaled")
})

ggplot(pred_temperature, aes(x = x, ymax = 10000*y_hat + year*0.1, ymin = 0 + year*0.1,
  group = year, colour = year, fill = year)) +
  geom_ribbon(lwd = 0.5, alpha = 0.2) +
  xlab("Temperature") +
  scale_color_viridis_c(option = "C") +
  scale_fill_viridis_c(option = "C") +
  ylab("Predicted density in some units") +
  gfplot::theme_pbs() +
  coord_cartesian(expand = FALSE)

# -----------------------------------
# functionize it

get_quants <- function(i) {
  f <- function(x) {
    exp(r$b_j[i] + r$b_rw_t[i,1] * x + r$b_rw_t[i,2] * x^2)
  }

  depth2scaled <- function(x) {
    (x - dat$temperature_mean[1]) / dat$temperature_sd[1]
  }

  total <- integrate(f,
    upper = depth2scaled(2),
    lower = depth2scaled(12))[[1]]

  thresh <- seq(2, 12, length.out = 1000)
  areas <- sapply(thresh, function(.x)
    integrate(f, upper = depth2scaled(.x), lower = depth2scaled(-1))[[1]])

  data.frame(
    lwr = thresh[max(which(areas/total < 0.25))],
    mid = thresh[max(which(areas/total < 0.50))],
    upr = thresh[max(which(areas/total < 0.75))],
    year = yrs[i])
}

out <- purrr::map_df(seq_along(unique(dat$year)), get_quants)

ggplot(out, aes(x = year, ymax = upr, ymin = lwr, y = mid)) +
  geom_ribbon(fill = "grey50", alpha = 0.4) +
  geom_line() +
  xlab("Year") +
  ylab("Depth (m)") +
  gfplot::theme_pbs()
