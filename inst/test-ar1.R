######################
# Estimate Centre of Gravity
######################

library(sdmTMB)
library(ggplot2)
library(dplyr)


pcod_spde <- make_spde(pcod$X, pcod$Y, n_knots = 100)

m <- sdmTMB(
 pcod, density ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
 time = "year", spde = pcod_spde, family = tweedie(link = "log"),
 silent = FALSE
)


# Predictions onto new data:
    # expand for time units:

newdata <- qcs_grid # use bathymetry grid of Queen Charlotte Sound
# ggplot(qcs_grid, aes(X, Y, colour=log(depth))) + geom_point() + scale_color_viridis_c(direction= -1)

original_time <- sort(unique(m$data[[m$time]])) # select and order all unique time periods

nd <- do.call("rbind",
  replicate(length(original_time), newdata, simplify = FALSE)) # repeat all rows of 'newdata' (grid of points in bathymetry layer) once for each year
nd[[m$time]] <- rep(original_time, each = nrow(newdata)) # label rows of newdata (now 'nd') with year ('time')

p <- predict(m, newdata = nd) # use model 'm' to predict values for each point



# get_index(p, value_name = "log_total", bias_correct = FALSE)


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



# Calculate centre of gravity for latitude and longitude

cog <- get_cog(p) # calculate centre of gravity for each data point
cog

ggplot(cog, aes(year, est, ymin = lwr, ymax = upr)) +
  geom_ribbon(alpha = 0.5) + geom_line() + facet_wrap(~coord, scales = "free_y") # what is free_y?


# table of COG by latitude
data.frame(Y = p$data$Y, est = exp(p$data$est), year = p$data$year) %>%
  group_by(year) %>% summarize(cog = sum(Y * est) / sum(est))







######################
# Simulation test of "real" vs. predicted values
######################

library(sdmTMB)
library(ggplot2)
library(dplyr)

set.seed(183228)


# simulate 'true' data
######################

## create fine-scale square 100 x 100 grid to predict on
grid <- expand.grid(X = seq(1:100), Y= seq(1:100))

## or use the boundaries of the Queen Charlotte Sound
# grid <- qcs_grid

simdat <- sim_args_vec( x=grid$X, y=grid$Y, time_steps = 9, plot = TRUE,
              ar1_fields = TRUE,
              ar1_phi = 0.5,
              sigma_O = 0.3,
              sigma_E = 0.3,
              kappa = 0.05,
              phi = 0.05
              )


# sub-sample from 'true' data
######################

dat <- simdat[[1]] %>% group_by(time) %>% sample_n(200) %>% ungroup()
# dat <- simdat %>% filter((x <25|x>50) & (y<25|y>50)) # try removing whole chunks of data


# model from sub-sample
######################

spde <- make_spde(x = dat$x, y = dat$y, n_knots = 200)

plot_spde(spde)

m <- sdmTMB(silent = FALSE,
  ar1_fields = TRUE,
  include_spatial = TRUE, # no fixed spatial random field
  data = dat, formula = z ~ 1, time = "time",
  family = gaussian(link = "identity"), spde = spde
)

# parameter estimates
######################

# if sim args are in list form
# inputs <- simdat[[2]]
#
# sigma_E_diff <- inputs$sigma_E - r$sigma_E # spatio-temporal random field SD
# kappa_diff <- inputs$kappa - exp(r$ln_kappa) # Matern for both space and space-time where distance at which ~%10 correlation = sqrt(Type(8.0)) / exp(ln_kappa)
# ar1_diff <- inputs$ar1_phi - (2 * plogis(m$model$par[['ar1_phi']]) - 1) # back transformed temporal correlation coefficent
# phi_diff <- inputs$phi - exp(r$ln_phi) # SD of observation error (aka sigma)
#


# if vector
inputs <- as.data.frame(reshape::melt (simdat[[2]]))
row.names(inputs)

r <- m$tmb_obj$report()
# r$range # distance at which ~%10 correlation

estimates <- c( ar1_phi = (2 * plogis(m$model$par[['ar1_phi']]) - 1), # back transformed temporal correlation coefficent
  sigma_O = r$sigma_O,
  sigma_E = r$sigma_E, # spatio-temporal random field SD
  kappa = exp(r$ln_kappa), # Matern for both space and space-time where distance at which ~%10 correlation = sqrt(Type(8.0)) / exp(ln_kappa)
  phi = exp(r$ln_phi) # SD of observation error (aka sigma)
)
estimates <- as.data.frame(reshape::melt(estimates))


parameter <- row.names(inputs)[1:5]
true <- inputs[1:5,]
par_diff <- cbind(parameter, true, estimates) %>% mutate (diff = value - true)

par_diff



# predict for all points in simulated data
######################

# replicate for each time period
original_time <- sort(unique(m$data[[m$time]]))
nd <- do.call("rbind",
  replicate(length(original_time), grid, simplify = FALSE))
nd[[m$time]] <- rep(original_time, each = nrow(grid))
# attributes(nd)

# run TMB with prediction turned on but replace sample 'dat' with new grid 'nd'
predictions <- predict(m, newdata = nd, xy_cols = c("X", "Y"))$data
head(predictions)



# combine true and predicted values for each point in space and time

spatial_bias_dat <- full_join(simdat, predictions, by=c("x"="X", "y"="Y", "time"="time"), suffix = c("", "_est")) %>% mutate (diff = est - real_z)
head(spatial_bias_dat)

hist(spatial_bias_dat$diff, breaks = 40)


# plot spatial distribution of model estimates beside simulated true values and the difference between them

plot_map_diff <- function(dat, id = c("x", "y", "time"), values = c("real_z","est", "diff"), time_periods = c("1","5","9") ) {

  melted <- reshape::melt(dat, id=id) %>% filter(variable %in% values) %>% filter(time %in% time_periods)

  ggplot(melted, aes_string("x", "y", fill = "value")) +
    geom_raster() +
    facet_grid(time~variable) +
    scale_fill_viridis_c() +
    coord_fixed()
}

plot_map_diff(spatial_bias_dat)
######################
######################




######################
# Save results of repeated simulations
######################

library(sdmTMB)
library(ggplot2)
library(dplyr)


# simulate 'true' data
######################


## create fine-scale square 100 x 100 grid to predict on
grid <- expand.grid(X = seq(1:100), Y= seq(1:100))
## or use the boundaries of the Queen Charlotte Sound
# grid <- qcs_grid

sim_runs <- function(x=grid$X, y=grid$Y, time_steps = 9, plot = TRUE,
                    ar1_fields = TRUE,
                    ar1_phi = 0.5,
                    sigma_O = 0.3,
                    sigma_E = 0.3,
                    kappa = 0.05,
                    phi = 0.05,
                    N = 500, n_knots = 200,
                    formula = z ~ 1, family = gaussian(link = "identity")) {

    simdat <- sim(x = x, y = y, time_steps = time_steps, plot = plot, ar1_fields = ar1_fields, ar1_phi = ar1_phi, sigma_O = sigma_O, sigma_E = sigma_E, kappa = kappa, phi = phi, seed = sample.int(1e6, 1))
    dat <- simdat %>% group_by(time) %>% sample_n(N) %>% ungroup() # sub-sample from 'true' data
    spde <- make_spde(x, y, n_knots)
    m <- sdmTMB(silent = FALSE,
                ar1_fields = ar1_fields,
                include_spatial = TRUE,
                data = dat, formula = formula, time = "time",
                family = family, spde = spde
                )
    r <- m$tmb_obj$report()
    estimates <- c( ar1_phi = (2 * plogis(m$model$par[['ar1_phi']]) - 1), # back transformed temporal correlation coefficent
                    sigma_O = r$sigma_O,
                    sigma_E = r$sigma_E, # spatio-temporal random field SD
                    kappa = exp(r$ln_kappa), # Matern for both space and space-time where distance at which ~%10 correlation = sqrt(Type(8.0)) / exp(ln_kappa)
                    phi = exp(r$ln_phi) # SD of observation error (aka sigma)
                    )
    inputs <- c(ar1_phi = ar1_phi,
                sigma_O = sigma_O,
                sigma_E = sigma_E,
                kappa = kappa,
                phi = phi
                )
    diff <- estimates - inputs
    diff
}

sim_results <- replicate(2, sim_runs(x=grid$X, y=grid$Y, time_steps = 9, plot = TRUE,
  ar1_fields = TRUE,
  ar1_phi = 0.5,
  sigma_O = 0.3,
  sigma_E = 0.3,
  kappa = 0.05,
  phi = 0.05,
  N = 500, n_knots = 200,
  formula = z ~ 1, family = gaussian(link = "identity")))


t(sim_results) # transpose data

# PROBLEM: exact same data and draw seem to be being used each time, why?


############ work in progress ##################


   # predict for all points in simulated data
   ######################

   # replicate for each time period
   original_time <- sort(unique(m$data[[m$time]]))
   nd <- do.call("rbind",
     replicate(length(original_time), grid, simplify = FALSE))
   nd[[m$time]] <- rep(original_time, each = nrow(grid))
   # attributes(nd)

   # run TMB with prediction turned on but replace sample 'dat' with new grid 'nd'
   predictions <- predict(m, newdata = nd, xy_cols = c("X", "Y"))$data
   head(predictions)



   # combine true and predicted values for each point in space and time

   spatial_bias_dat <- full_join(simdat, predictions, by=c("x"="X", "y"="Y", "time"="time"), suffix = c("", "_est")) %>% mutate (diff = est - real_z)
   head(spatial_bias_dat)

   hist(spatial_bias_dat$diff, breaks = 40)





     )

   ######################
   ######################
   ######################
   ######################


