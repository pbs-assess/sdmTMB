# Simulation test for bias in AR1 sdmTMB models

######################
######################

# Steps first executed independently and then in functions that repeat simulation and compile results

library(sdmTMB)
library(ggplot2)
library(dplyr)
library(purrr)
library(future)


# Set spatial grid (also run this section before running functions)
######################

## create fine-scale square 100 x 100 grid to predict on
grid <- expand.grid(X = seq(1:40), Y = seq(1:40))

## or use the boundaries of the Queen Charlotte Sound
# grid <- qcs_grid


# 1. Simulate 'true' data
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

# 2. Sub-sample from 'true' data
######################

dat <- simdat[[1]] %>% group_by(time) %>% sample_n(1500) %>% ungroup()
# dat <- simdat %>% filter((x <25|x>50) & (y<25|y>50)) # try removing whole chunks of data


# 3. Model from sub-sample
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


# 4. Compare parameter estimates with input values
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

# 5. Use model to predict for all grid points in simulated data
######################

# replicate for each time period
original_time <- sort(unique(m$data[[m$time]]))
nd <- do.call("rbind", replicate(length(original_time), grid, simplify = FALSE))
nd[[m$time]] <- rep(original_time, each = nrow(grid))

# run TMB with prediction turned on but replace sample 'dat' with new grid 'nd'
predictions <- predict(m, newdata = nd, xy_cols = c("X", "Y"))$data
head(predictions)


# 6. Contrast true and predicted values for each point in space and time
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
#' m <- sdmTMB(data = dat,
#'             formula = z ~ 1,
#'             time = "time",
#'             family = gaussian(link = "identity"),
#'             spde = spde)
#' original_time <- sort(unique(m$data[[m$time]]))
#' nd <- do.call("rbind", replicate(length(original_time), grid, simplify = FALSE))
#' nd[[m$time]] <- rep(original_time, each = nrow(grid))
#' predictions <- predict(m, newdata = nd, xy_cols = c("X", "Y"))$data
#' spatial_bias_dat <- full_join(simdat[[1]], predictions,
#'                               by = c("x" = "X", "y" = "Y", "time" = "time"),
#'                               suffix = c("", "_est")) %>%
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

# Save parameter estimates of repeated simulations

######################
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

# Repeat simulations
#           j times for each set of parameter combinations
#           save lists of tibbles

######################
######################

#' Function that saves a list of tibbles of:
#'          [[1]] parameter inputs and estimates &
#'          [[2]] real values and predicted values
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
#'
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
#' @param j Number of runs to be conducted for each unique combination of parameter values
#'
#' @return
#' @export
#'
#' @examples
#' grid <- expand.grid(X = seq(1:40), Y = seq(1:40))
#' args <- generate_arg(kappa = c(0.005, 0.1, 1), time_steps = 3, repeats = 3L)
#' all_iter <- purrr::pmap(args, sim_predictions,
#'                         grid = grid, x = grid$X, y = grid$Y,
#'                         formula = z ~ 1, family = gaussian(link = "identity"))
#' predictions <- all_iter %>% map(~ .x[["predicted"]])
#' spatial_bias_plots <- purrr::map(predictions, plot_map_diff, time_periods = c(1,2,3))
#' params <- all_iter %>% map(~ .x[["par"]]) %>%
#'                        bind_rows()
#' par_diff <- params %>% group_by(iter, parameter) %>%
#'                        mutate(diff = inputs - estimates) %>%
#'                        ungroup()
#'
generate_arg <- function(ar1_phi = 0.4, # c(-0.85, 0.1, 0.85),
                         sigma_O = 0.3,
                         sigma_E = 0.3,
                         kappa = 0.1, # c(0.005, 0.1, 1),
                         phi = 0.05, # c(0.01, 0.1),
                         time_steps = 3,
                         N = 1000,
                         n_knots = 250,
                         j = 20L) {
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
  arguments$count <- j
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

# Explore and plot simulation results

######################
######################


# result list of tibbles of predictions or remove # and combine to one tibble
predictions <- all_iter %>% map(~ .x[["predicted"]]) # %>% bind_rows(.id = "iter")

# each run can be plotted
spatial_bias_plots <- purrr::map(predictions, plot_map_diff, time_periods = c(1,2,3))
pdf("spatial_bias_plots.pdf")
spatial_bias_plots
dev.off()

# result tibble of parameter estimates
params <- all_iter %>% map(~ .x[["par"]]) %>%
                       bind_rows()

# Calculate difference between parameter value input into simulation and estimate based on sdmTMB model
par_diff <- params %>% group_by(parameter) %>%
                       mutate(sd_est = sd(estimates), n = n()) %>%
                       group_by(iter, parameter) %>%
                       mutate(std_diff = (inputs - estimates)/sd_est) %>%
                       ungroup()

#' Plot histograms of parameter estimates from n simulations
#'
#' @param data Dataframe containing all simulated parameter estimates
#' @param x Varible to be plotted (Default = data$std_diff)
#' @param xlabel Description of variable to be plotted for use on x axis label
#' @param fill Varible used to colour bars to indicate if some estimates should be trusted more than others (Default = data$converg)
#' @param notes Description of fill choice or other caveats
#' @param bins Number of bins in histogram (Default = n/4)
#'
#' @return
#' @export
#'
#' @examples
#' params <- all_iter %>% map(~ .x[["par"]]) %>%
#'                        bind_rows()
#' par_diff <- params %>% group_by(parameter) %>%
#'                        mutate(sd_est = sd(estimates), n = n()) %>%
#'                        group_by(iter, parameter) %>%
#'                        mutate(std_diff = (inputs - estimates)/sd_est) %>%
#'                        ungroup()
#' par_error_hist(par_diff)
#'
par_error_hist <- function(data = par_diff,
                           x = data$std_diff,
                           xlabel = "Relative difference from input value",
                           fill = data$converg,
                           notes = "Note: if 2 colours, than some models did not converg",
                           bins = n/4){
  n <- data$n[1]
  fill <- as.factor(fill)
  ggplot(data, aes(x = x, fill = fill)) +
    geom_histogram(bins = n/4)  +
    #scale_fill_viridis_d() +
    geom_vline(xintercept = 0, linetype="dashed") +
    labs(title = "Simulated parameter estimates", x = xlabel, caption = notes) +
    facet_wrap(~parameter, scales = "free_x") +
    theme(legend.position="none", plot.caption=element_text(size=12))
}

par_error_hist(par_diff)
