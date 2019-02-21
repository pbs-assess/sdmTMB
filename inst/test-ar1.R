pcod_spde <- make_spde(pcod$X, pcod$Y, n_knots = 100)
m <- sdmTMB(
 pcod, density ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
 time = "year", spde = pcod_spde, family = tweedie(link = "log"),
 silent = FALSE
)

# Predictions at original data locations:
predictions <- predict(m)$data
cols <- c("year", "X", "Y", "est", "est_fe",
  "est_re_s", "est_re_st", "s_i")
head(predictions[,cols])

predictions$resids <- residuals(m) # randomized quantile residuals
library(ggplot2)
ggplot(predictions, aes(X, Y, col = resids)) + scale_colour_gradient2() +
  geom_point() + facet_wrap(~year)
hist(predictions$resids)
qqnorm(predictions$resids);abline(a = 0, b = 1)

# Predictions onto new data:
    # expand for time units:
newdata <- qcs_grid
original_time <- sort(unique(m$data[[m$time]]))
nd <- do.call("rbind",
  replicate(length(original_time), newdata, simplify = FALSE))
nd[[m$time]] <- rep(original_time, each = nrow(newdata))
predictions <- predict(m, newdata = nd)$data

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

set.seed(183228)
dat <- sim(
  time_steps = 9, ar1_fields = TRUE, ar1_phi = 0.8,
  plot = TRUE, sigma_O = 0.001, sigma_E = 0.3, phi = 0.05
)
spde <- make_spde(x = dat$x, y = dat$y, n_knots = 200)
m <- sdmTMB(silent = FALSE, ar1_fields = TRUE, include_spatial = FALSE,
  data = dat, formula = z ~ 1, time = "time",
  family = gaussian(link = "identity"), spde = spde
)
newdata <- expand.grid(x = seq(min(dat$x), max(dat$x), length.out = 100),
  y = seq(min(dat$y), max(dat$y), length.out = 100))
original_time <- sort(unique(m$data[[m$time]]))
nd <- do.call("rbind",
  replicate(length(original_time), newdata, simplify = FALSE))
nd[[m$time]] <- rep(original_time, each = nrow(newdata))
predictions <- predict(m, newdata = nd, xy_cols = c("x", "y"))$data

plot_map <- function(dat, column = "est") {
  ggplot(dat, aes_string("x", "y", fill = column)) +
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
