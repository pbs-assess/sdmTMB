######################
# Simulation test for bias in AR1 sdmTMB models
######################

library(sdmTMB)
library(ggplot2)
library(dplyr)
library(purrr)
library(future)

# Steps first executed independently and then in functions that repeat simulation and compile results

# 1. Set spatial grid ( also run this section before running functions)
######################

## create fine-scale square 100 x 100 grid to predict on
grid <- expand.grid(X = seq(1:40), Y = seq(1:40))

## or use the boundaries of the Queen Charlotte Sound
# grid <- qcs_grid


# 2. Simulate 'true' data
######################


#' Simulate data in list format with a vector of input values
#'
#' @param x A vector of x coordinates.
#' @param y A vector of y coordinates.
#' @param time_steps The number of time steps.
#' @param ar1_fields Logical for whether or not to include AR1 structure.
#' @param ar1_phi Correlation between years; should be between -1 and 1.
#' @param sigma_O SD of spatial process (Omega).
#' @param sigma_E SD of spatiotemporal process (Epsilon).
#' @param kappa Parameter that controls the decay of spatial correlation (3/kappa is roughly the distance at which points are %10 correlated)
#' @param phi Observation error scale parameter.
#' @param plot Logical for whether or not to produce a plot.
#'
#' @return
#' @export
#'
#' @examples
sim_args_vec <- function(x = stats::runif(400, 0, 10), y = stats::runif(400, 0, 10),
                         time_steps = 1L, ar1_fields = FALSE, ar1_phi = 0.5,
                         sigma_O = 0.4, sigma_E = 0.3, kappa = 1.3, phi = 0.2,
                         plot = FALSE) {
  d <- sim(
    x = x, y = y, time_steps = time_steps, ar1_fields = ar1_fields, plot = plot,
    ar1_phi = ar1_phi, sigma_O = sigma_O, sigma_E = sigma_E, kappa = kappa, phi = phi
  )

  list(d, inputs = c(
    ar1_phi = ar1_phi,
    sigma_O = sigma_O,
    sigma_E = sigma_E,
    kappa = kappa,
    phi = phi
  ))
}

simdat <- sim_args_vec(
  x = grid$X, y = grid$Y,
  time_steps = 9,
  plot = TRUE,
  ar1_fields = TRUE,
  ar1_phi = 0.5,
  sigma_O = 0.3,
  sigma_E = 0.3,
  kappa = 0.05,
  phi = 0.05
)

# 3. Sub-sample from 'true' data
######################

dat <- simdat[[1]] %>% group_by(time) %>% sample_n(1500) %>% ungroup()
# dat <- simdat %>% filter((x <25|x>50) & (y<25|y>50)) # try removing whole chunks of data


# 4. Model from sub-sample
######################

spde <- make_spde(x = dat$x, y = dat$y, n_knots = 150)
plot_spde(spde)

m <- sdmTMB(
  silent = FALSE,
  ar1_fields = TRUE,
  include_spatial = TRUE, # no fixed spatial random field
  data = dat, formula = z ~ 1, time = "time",
  family = gaussian(link = "identity"), spde = spde
)

r <- m$tmb_obj$report()
# r$range # distance at which ~%10 correlation


# 5. Compare parameter estimates with input values
######################

# Input values (from vectorized simdat list element)
inputs <- as.data.frame(reshape::melt(simdat[[2]]))

# back transform parameter estimates from model report() for comparison with input values
estimates <- c(
  ar1_phi = (2 * plogis(m$model$par[["ar1_phi"]]) - 1),
  sigma_O = r$sigma_O,
  sigma_E = r$sigma_E,
  kappa = exp(r$ln_kappa),
  phi = exp(r$ln_phi)
)

estimates <- as.data.frame(reshape::melt(estimates))

inputs
estimates

# 6. Use model to predict for all grid points in simulated data
######################

# replicate for each time period
original_time <- sort(unique(m$data[[m$time]]))
nd <- do.call("rbind", replicate(length(original_time), grid, simplify = FALSE))
nd[[m$time]] <- rep(original_time, each = nrow(grid))

# run TMB with prediction turned on but replace sample 'dat' with new grid 'nd'
predictions <- predict(m, newdata = nd, xy_cols = c("X", "Y"))$data
head(predictions)


# 7. Contrast true and predicted values for each point in space and time
######################

spatial_bias_dat <- full_join(simdat[[1]], predictions, by = c("x" = "X", "y" = "Y", "time" = "time"), suffix = c("", "_est")) %>% mutate(diff = est - real_z)
head(spatial_bias_dat)

hist(spatial_bias_dat$diff, breaks = 40)

#' Function to plot spatial distribution of model estimates beside simulated true values and the difference between them
#'
#' @param dataframe Dataframe containing all simulated and predicted values to be plotted spatially
#' @param id List of columns in dataframe that together identify unique observations (Default = c("x", "y", "time"))
#' @param values List of values to be shown in plotted on separate panels (Default = c("real_z", "est", "diff"))
#' @param time_periods List of time periods to be shown in plot (Default = c("1", "5", "9"))
#'
#' @return
#' @export
#'
#' @examples
#' simdat <- sim_args_vec()
#' dat <- simdat[[1]] %>% group_by(time) %>% sample_n(1500) %>% ungroup()
#' spde <- make_spde(x = dat$x, y = dat$y, n_knots = 150)
#' m <- sdmTMB(data = dat, formula = z ~ 1, time = "time", family = gaussian(link = "identity"), spde = spde)
#' original_time <- sort(unique(m$data[[m$time]]))
#' nd <- do.call("rbind", replicate(length(original_time), grid, simplify = FALSE))
#' nd[[m$time]] <- rep(original_time, each = nrow(grid))
#' predictions <- predict(m, newdata = nd, xy_cols = c("X", "Y"))$data
#' spatial_bias_dat <- full_join(simdat[[1]], predictions, by = c("x" = "X", "y" = "Y", "time" = "time"), suffix = c("", "_est")) %>%
#'                     mutate(diff = est - real_z)
#' plot_map_diff(spatial_bias_dat)

plot_map_diff <- function(dataframe,
                          id = c("x", "y", "time"),
                          values = c("real_z", "est", "diff"),
                          time_periods = c("1", "5", "9")) {

  melted <- reshape2::melt(dataframe, id) %>% # could be replace with tidyr::gather(spatial_bias_dat, "variable", "value",... )?
    filter(variable %in% values) %>%
    filter(time %in% time_periods)

  ggplot(melted, aes_string("x", "y", fill = "value")) +
    geom_raster() +
    facet_grid(time ~ variable) +
    scale_fill_viridis_c() +
    coord_fixed()
}

plot_map_diff(spatial_bias_dat)


######################
######################

######################
# Save parameter estimates of repeated simulations
######################

######################
######################

# simulate 'true' data
######################

#' Function that saves a tibble of parameter inputs and estimates for a single iteration
#'
#' @param iter Iteration id number; default is a random number; used to set.seed
#' @param grid Dataframe of spatial coordinates eg. c(X, Y)
#' @param x A vector of x coordinates.
#' @param y A vector of y coordinates.
#' @param time_steps The number of time steps.
#' @param plot Logical for whether or not to produce a plot.
#' @param ar1_fields Logical for whether or not to include AR1 structure.
#' @param ar1_phi Correlation between years; should be between -1 and 1.
#' @param sigma_O SD of spatial process (Omega).
#' @param sigma_E SD of spatiotemporal process (Epsilon).
#' @param kappa Parameter that controls the decay of spatial correlation (3/kappa is roughly the distance at which points are %10 correlated)
#' @param phi Observation error scale parameter.
#' @param N Sub-sample size = number of observations included in the sdmTMB model
#' @param n_knots Number of knots for spatial process#' @param formula
#' @param formula Define model to be assessed with the simulated data
#' @param family Set family of model to be assessed
#'
#' @return
#' @export
#'
#' @examples
sim_parameters <- function(iter = sample.int(1e3, 1), plot = TRUE,
                           grid = grid, x = grid$X, y = grid$Y,
                           time_steps = 9,
                           ar1_fields = TRUE,
                           ar1_phi = 0.5,
                           sigma_O = 0.3,
                           sigma_E = 0.3,
                           kappa = 0.05,
                           phi = 0.05,
                           N = 1000, n_knots = 150,
                           formula = z ~ 1, family = gaussian(link = "identity")) {
  set.seed(iter * 581267)

  simdat <- sdmTMB::sim(
    x = x, y = y, time_steps = time_steps, plot = plot,
    ar1_fields = ar1_fields, ar1_phi = ar1_phi, sigma_O = sigma_O, sigma_E = sigma_E, kappa = kappa, phi = phi
  )

  dat <- simdat %>% group_by(time) %>% sample_n(N) %>% ungroup() # sub-sample from 'true' data

  spde <- make_spde(dat$x, dat$y, n_knots)
  plot_spde(spde)
  # browser()

  m <- sdmTMB(
    silent = FALSE,
    ar1_fields = ar1_fields,
    include_spatial = TRUE,
    data = dat, formula = formula, time = "time",
    family = family, spde = spde
  )

  r <- m$tmb_obj$report()

  estimates <- c(
    ar1_phi = (2 * plogis(m$model$par[["ar1_phi"]]) - 1),
    sigma_O = r$sigma_O,
    sigma_E = r$sigma_E,
    kappa = exp(r$ln_kappa),
    phi = exp(r$ln_phi)
  )

  inputs <- c(
    ar1_phi = ar1_phi,
    sigma_O = sigma_O,
    sigma_E = sigma_E,
    kappa = kappa,
    phi = phi
  )

  parameter <- names(inputs)
  converg <- m$model$convergence
  # diff <- estimates - inputs
  run <- tibble(parameter = parameter, inputs = inputs, estimates = estimates, iter = iter, converg = converg)
  run
}

single_run <- sim_parameters(
  iter = sample.int(1e3, 1), grid = grid, x = grid$X, y = grid$Y, time_steps = 6, plot = TRUE,
  ar1_fields = TRUE,
  ar1_phi = 0.7,
  sigma_O = 0.2,
  sigma_E = 0.2,
  kappa = 0.05,
  phi = 0.01,
  formula = z ~ 1, family = gaussian(link = "identity")
)

single_run


######################
######################

######################
# Add predictions for all points in simulated data
######################

######################
######################

#' Function that saves a list of tibble of parameter inputs and estimates for a single iteration
#'
#' @param iter Iteration id number; default is a random number; used to set.seed
#' @param grid Dataframe of spatial coordinates eg. c(X, Y)
#' @param x A vector of x coordinates.
#' @param y A vector of y coordinates.
#' @param time_steps The number of time steps.
#' @param plot Logical for whether or not to produce a plot.
#' @param ar1_phi Correlation between years; should be between -1 and 1.
#' @param sigma_O SD of spatial process (Omega).
#' @param sigma_E SD of spatiotemporal process (Epsilon).
#' @param kappa Parameter that controls the decay of spatial correlation (3/kappa is roughly the distance at which points are %10 correlated)
#' @param phi Observation error scale parameter.
#' @param N Sub-sample size = number of observations included in the sdmTMB model
#' @param n_knots Number of knots for spatial process#' @param formula
#' @param formula Define model to be assessed with the simulated data
#' @param family Set family of model to be assessed
#'
#' @return
#' @export
#'
#' @examples
sim_predictions <- function(iter = sample.int(1e3, 1), plot = TRUE,
                            grid = grid, x = grid$X, y = grid$Y,
                            time_steps = 3,
                            ar1_fields = TRUE,
                            ar1_phi = 0.5,
                            sigma_O = 0.3,
                            sigma_E = 0.3,
                            kappa = 0.05,
                            phi = 0.05,
                            N = 1000, n_knots = 150, # iter.max=1e4, eval.max=1e4,
                            formula = z ~ 1, family = gaussian(link = "identity")) {
  set.seed(iter * 581267)

  simdat <- sim(
    plot = plot, x = x, y = y,
    time_steps = time_steps, ar1_fields = ar1_fields,
    ar1_phi = ar1_phi, sigma_O = sigma_O, sigma_E = sigma_E, kappa = kappa, phi = phi
  )

  dat <- simdat %>% group_by(time) %>% sample_n(N) %>% ungroup() # sub-sample from 'true' data
  spde <- make_spde(dat$x, dat$y, n_knots)
  plot_spde(spde)
  # browser()

  m <- sdmTMB(
    silent = FALSE,
    ar1_fields = ar1_fields,
    include_spatial = TRUE,
    data = dat, formula = formula, time = "time",
    family = family, spde = spde
  )

  r <- m$tmb_obj$report()

  estimates <- c(
    ar1_phi = (2 * plogis(m$model$par[["ar1_phi"]]) - 1),
    sigma_O = r$sigma_O,
    sigma_E = r$sigma_E,
    kappa = exp(r$ln_kappa),
    phi = exp(r$ln_phi)
  )

  inputs <- c(
    ar1_phi = ar1_phi,
    sigma_O = sigma_O,
    sigma_E = sigma_E,
    kappa = kappa,
    phi = phi
  )
  parameter <- names(inputs)
  converg <- m$model$convergence

  # replicate grid for each time period
  original_time <- sort(unique(m$data[[m$time]]))
  nd <- do.call("rbind", replicate(length(original_time), grid, simplify = FALSE))
  nd[[m$time]] <- rep(original_time, each = nrow(grid))

  # run TMB with prediction turned on but replace sample 'dat' with new grid 'nd'
  predictions <- predict(m, newdata = nd, xy_cols = c("X", "Y"))$data
  head(predictions)

  # combine true and predicted values for each point in space and time
  spatial_bias <- full_join(simdat, predictions, by = c("x" = "X", "y" = "Y", "time" = "time"), suffix = c("", "_est")) %>% mutate(diff = est - real_z)
  spatial_bias$converg <- converg
  run <- list(par = tibble(parameter = parameter, inputs = inputs, estimates = estimates, iter = iter, converg = converg), predicted = as_tibble(spatial_bias))
  run
}

#' Generate dataframe with parameter arguments for multiple simulation runs
#'
#' Note: This function creates a dataframe of all possible combinations,
#' therefore these parameters can be fixed to a single value, or vary using c(value1, value2, value3, ...)
#' @param ar1_phi Correlation between years; should be between -1 and 1.
#' @param sigma_O SD of spatial process (Omega).
#' @param sigma_E SD of spatiotemporal process (Epsilon).
#' @param kappa Parameter that controls the decay of spatial correlation (3/kappa is roughly the distance at which points are %10 correlated)
#' @param phi Observation error scale parameter.
#' @param time_steps The number of time steps.
#' @param N Sub-sample size = number of observations included in the sdmTMB model
#' @param n_knots Number of knots for spatial process
#' @param repeats Number of runs to be conducted for each unique combination of parameter values
#'
#' @return
#' @export
#'
#' @examples
generate_arg <- function(ar1_phi = 0.4, # c(-0.85, 0.1, 0.85),
                         sigma_O = 0.3,
                         sigma_E = 0.3,
                         kappa = 0.1, # c(0.005, 0.1, 1),
                         phi = 0.05, # c(0.01, 0.1),
                         time_steps = 3,
                         N = 1000,
                         n_knots = 250,
                         repeats = 20L) {
  arguments <- expand.grid(
    ar1_phi = ar1_phi,
    sigma_O = sigma_O,
    sigma_E = sigma_E,
    kappa = kappa,
    phi = phi,
    N = N,
    n_knots = n_knots
  )
  # browser()
  nrow(arguments)
  arguments$count <- repeats
  arguments <- arguments[rep(seq_len(nrow(arguments)), arguments$count), ]
  arguments_apply <- dplyr::select(arguments, -count)
  nrow(arguments_apply)
  arguments_apply$iter <- 1:nrow(arguments_apply)
  arguments_apply
}

args <- generate_arg()
glimpse(args)

all_iter <- purrr::pmap(args, sim_predictions,
  grid = grid, x = grid$X, y = grid$Y,
  formula = z ~ 1, family = gaussian(link = "identity")
)

# Paralell process (but not currently working)
# library(future)
# plan(multisession)
# all_iter <- purrr::pmap(args, ~ future(sim_predictions(.x, grid = grid, x = grid$X, y = grid$Y,
#   plot = TRUE, formula = z ~ 1, family = gaussian(link = "identity")
# )))

######################
######################

######################
# Explore and plot simulation results
######################

######################
######################


# result list of tibbles of predictions or remove # and combine to one tibble
predictions <- all_iter %>% map(~ .x[["predicted"]]) # %>% bind_rows(.id = "iter")

# each run can be plotted
plot_map_diff(predictions[[1]], time_periods = c(1,2,3))


# result tibble of parameter estimates
params <- all_iter %>% map(~ .x[["par"]]) %>% bind_rows()

# Calculate difference between parameter value input into simulation and estimate based on sdmTMB model
par_diff <- params %>% group_by(iter, parameter) %>% mutate(diff = inputs - estimates) %>% ungroup()
# past_run <- par_diff

# par_diff_converged <- par_diff %>% filter(converg == 0) %>% ungroup()
# par_diff_fail <- par_diff %>% filter(converg == 1) %>% ungroup()

ggplot(par_diff, aes(x = diff, fill = as.factor(converg))) + geom_histogram() + xlim(-0.25, 0.25) + facet_wrap(~parameter) # , scales = "free_x")

kappa <- par_diff %>% filter(parameter == "kappa")
ggplot(data = kappa, aes(x = estimates)) + geom_histogram(bins = 50) + xlim(-0.5, 0.5)
ggplot(data = kappa, aes(x = diff)) + geom_histogram(bins = 50) #+ xlim(-0.5,0.5)

phi <- par_diff %>% filter(parameter == "phi")
ggplot(data = phi, aes(x = estimates)) + geom_histogram(bins = 50) #+ xlim(-0.5,0.5)

sigma_E <- par_diff %>% filter(parameter == "sigma_E")
ggplot(data = sigma_E, aes(x = diff)) + geom_histogram(bins = 50) #+ xlim(-0.5,0.5)

# make sure enough time steps are included if this is to be estimated?
ar1_phi <- par_diff %>% filter(parameter == "ar1_phi")
ggplot(data = ar1_phi, aes(x = diff)) + geom_histogram(bins = 50) #+ xlim(-0.5,0.5)
