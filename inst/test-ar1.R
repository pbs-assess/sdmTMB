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


# get_index(p, value_name = "log_total", bias_correct = FALSE)

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
  plot = TRUE, sigma_O = 0.001, sigma_E = 0.3, phi = 0.05
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
  plot = TRUE, sigma_O = 0.01, sigma_E = 0.3, phi = 0.01
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
