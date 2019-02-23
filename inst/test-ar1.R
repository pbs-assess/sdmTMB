
pcod_spde <- make_spde(pcod$X, pcod$Y, n_knots = 100)

m <- sdmTMB(
 pcod, density ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
 time = "year", spde = pcod_spde, family = tweedie(link = "log"),
 silent = FALSE
)

# Predictions onto new data:
    # expand for time units:

newdata <- qcs_grid # use bathymetry grid of Queen Charlotte Sound
ggplot(qcs_grid, aes(X, Y, colour=log(depth))) + geom_point() + scale_color_viridis_c(direction= -1)

original_time <- sort(unique(m$data[[m$time]])) # select and order all unique time periods

nd <- do.call("rbind",
  replicate(length(original_time), newdata, simplify = FALSE)) # repeat all rows of 'newdata' (grid of points in bathymetry layer) once for each year
nd[[m$time]] <- rep(original_time, each = nrow(newdata)) # label rows of newdata (now 'nd') with year ('time')

p <- predict(m, newdata = nd) # use model 'm' to predict values for each point
cog <- get_cog(p) # calculate centre of gravity for each data point
cog

library(ggplot2)
ggplot(cog, aes(year, est, ymin = lwr, ymax = upr)) +
  geom_ribbon(alpha = 0.5) + geom_line() + facet_wrap(~coord, scales = "free_y")

# get_index(p, value_name = "log_total", bias_correct = FALSE)

library(dplyr)
data.frame(Y = p$data$Y, est = exp(p$data$est), year = p$data$year) %>%
  group_by(year) %>% summarize(cog = sum(Y * est) / sum(est))

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



######################
# Simulation
######################
library(ggplot2)
set.seed(183228)


# simulate 'true' data
######################

## create fine-scale square grid to predict on
grid <- expand.grid(X = seq(1:100), Y= seq(1:100))

## or use the boundaries of the Queen Charlotte Sound
# grid <- qcs_grid

simdat <- sim( x=grid$X, y=grid$Y,
  time_steps = 9, ar1_fields = FALSE, ar1_phi = 0.5,
  plot = TRUE, sigma_O = 0.3, sigma_E = 0.3, kappa = 0.05, phi = 0.05
)


# sub-sample from 'true' data
######################

library(dplyr)
dat <- simdat %>% group_by(time) %>% sample_n(9000) %>% ungroup()
# dat <- simdat %>% filter((x <25|x>50) & (y<25|y>50)) # try removing whole chunks of data

spde <- make_spde(x = dat$x, y = dat$y, n_knots = 200)

plot_spde(spde)

m <- sdmTMB(silent = FALSE,
  ar1_fields = TRUE,
  include_spatial = TRUE, # no fixed spatial random field
  data = dat, formula = z ~ 1, time = "time",
  family = gaussian(link = "identity"), spde = spde
)


# replicate for each time period
original_time <- sort(unique(m$data[[m$time]]))
nd <- do.call("rbind",
  replicate(length(original_time), grid, simplify = FALSE))
nd[[m$time]] <- rep(original_time, each = nrow(grid))
attributes(nd)
# run TMB with prediction turned on but replace sample 'dat' with new grid 'nd'
predictions <- predict(m, newdata = nd, xy_cols = c("X", "Y"))$data
head(predictions)


# model output
r <- m$tmb_obj$report()
names(r)
r$sigma_E # spatio-temporal random field SD
exp(r$ln_kappa) # Matern for both space and space-time where distance at which ~%10 correlation = sqrt(Type(8.0)) / exp(ln_kappa)
r$range
2 * plogis(m$model$par[['ar1_phi']]) - 1 # back transformed temporal correlation coefficent
exp(r$ln_phi) # SD of observation error (aka sigma)




plot_map <- function(dat, fill = "est") {
  ggplot(dat, aes_string("X", "Y", fill = fill)) +
    geom_raster() +
    facet_wrap(~time) +
    coord_fixed()
}


ggplot(simdat, aes_string("x", "y", fill = "real_z")) +
  geom_raster() +
  facet_wrap(~time) +
  coord_fixed() +
  scale_fill_viridis_c() +
  ggtitle("Simulated 'real' values")


plot_map(predictions, "est") +
  scale_fill_viridis_c() +
  ggtitle("Prediction (fixed effects + all random effects)")




######################
######################
######################
######################



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
newdata <- expand.grid(x = seq(min(grid$x), max(grid$x), length.out = 100),
  y = seq(min(grid$y), max(grid$y), length.out = 100))

# replicate for each time period
original_time <- sort(unique(m$data[[m$time]]))
nd <- do.call("rbind",
  replicate(length(original_time), newdata, simplify = FALSE))
nd[[m$time]] <- rep(original_time, each = nrow(newdata))

# run TMB with prediction turned on
predictions <- predict(m, newdata = nd, xy_cols = c("x", "y"))$data
?predict

# contrast predictions with simulated true data


r <- m$tmb_obj$report()
r$sigma_E # spatio-temporal
exp(r$ln_kappa) # Matern for both space and space-time
2 * plogis(m$model$par[['ar1_phi']]) - 1 # back transform
