######################
# Simulation test of "real" vs. predicted values
######################

library(sdmTMB)
library(ggplot2)
library(dplyr)
library(purrr)
library(future)

# set.seed(183228)


# simulate 'true' data
######################

## create fine-scale square 100 x 100 grid to predict on
grid <- expand.grid(X = seq(1:50), Y = seq(1:50))

## or use the boundaries of the Queen Charlotte Sound
# grid <- qcs_grid

sim_args_vec <- function(x = stats::runif(400, 0, 10), y = stats::runif(400, 0, 10),
  time_steps = 1L, ar1_fields = FALSE, ar1_phi = 0.5,
  sigma_O = 0.4, sigma_E = 0.3, kappa = 1.3, phi = 0.2,
   seed = sample.int(1e6, 1), plot = FALSE) {
  d = sim(x, y,
     time_steps, ar1_fields, ar1_phi,
     sigma_O, sigma_E, kappa , phi,

    seed, plot )
  list(d, inputs = c(ar1_phi = ar1_phi,
    sigma_O =sigma_O,
    sigma_E = sigma_E,

    kappa = kappa,
    phi = phi,
    seed = seed))
}

simdat <- sim_args_vec(
  x = grid$X, y = grid$Y, time_steps = 3, plot = TRUE,
  ar1_fields = TRUE,
  ar1_phi = 0.5,
  sigma_O = 0.3,
  sigma_E = 0.3,
  kappa = 0.05,
  phi = 0.05
)

# sub-sample from 'true' data
######################

dat <- simdat[[1]] %>% group_by(time) %>% sample_n(2000) %>% ungroup()
# dat <- simdat %>% filter((x <25|x>50) & (y<25|y>50)) # try removing whole chunks of data


# model from sub-sample
######################

spde <- make_spde(x = dat$x, y = dat$y, n_knots = 200)

plot_spde(spde)

m <- sdmTMB(
  silent = FALSE,
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
inputs <- as.data.frame(reshape::melt(simdat[[2]]))
row.names(inputs)

r <- m$tmb_obj$report()
# r$range # distance at which ~%10 correlation

estimates <- c(
  ar1_phi = (2 * plogis(m$model$par[["ar1_phi"]]) - 1), # back transformed temporal correlation coefficent
  sigma_O = r$sigma_O,
  sigma_E = r$sigma_E, # spatio-temporal random field SD
  kappa = exp(r$ln_kappa), # Matern for both space and space-time where distance at which ~%10 correlation = sqrt(Type(8.0)) / exp(ln_kappa)
  phi = exp(r$ln_phi) # SD of observation error (aka sigma)
)
estimates <- as.data.frame(reshape::melt(estimates))

# predict for all points in simulated data
######################

# replicate for each time period
original_time <- sort(unique(m$data[[m$time]]))
nd <- do.call(
  "rbind",
  replicate(length(original_time), grid, simplify = FALSE)
)
nd[[m$time]] <- rep(original_time, each = nrow(grid))
# attributes(nd)

# run TMB with prediction turned on but replace sample 'dat' with new grid 'nd'
predictions <- predict(m, newdata = nd, xy_cols = c("X", "Y"))$data
head(predictions)

# combine true and predicted values for each point in space and time

spatial_bias_dat <- full_join(simdat[[1]], predictions, by = c("x" = "X", "y" = "Y", "time" = "time"), suffix = c("", "_est")) %>% mutate(diff = est - real_z)
head(spatial_bias_dat)

hist(spatial_bias_dat$diff, breaks = 40)


# plot spatial distribution of model estimates beside simulated true values and the difference between them

plot_map_diff <- function(dat, id = c("x", "y", "time"), values = c("real_z", "est", "diff"), time_periods = c("1", "5", "9")) {
  melted <- reshape::melt(dat, id = id) %>% filter(variable %in% values) %>% filter(time %in% time_periods)

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


library(sdmTMB)
library(ggplot2)
library(dplyr)
library(tibble)

# simulate 'true' data
######################


## create fine-scale square 100 x 100 grid to predict on
grid <- expand.grid(X = seq(1:50), Y = seq(1:50))
## or use the boundaries of the Queen Charlotte Sound
# grid <- qcs_grid

model_sim <- function(iter = sample.int(1e3, 1), grid = grid, x = grid$X, y = grid$Y, time_steps = 3, plot = TRUE,
                      ar1_fields = TRUE,
                      ar1_phi = 0.5,
                      sigma_O = 0.3,
                      sigma_E = 0.3,
                      kappa = 0.05,
                      phi = 0.05,
                      N = 100, n_knots = 100,
                      formula = z ~ 1, family = gaussian(link = "identity")) {
  set.seed(iter * 581267)
  simdat <- sim(
    x = x, y = y, time_steps = time_steps, plot = plot,
    ar1_fields = ar1_fields, ar1_phi = ar1_phi, sigma_O = sigma_O, sigma_E = sigma_E, kappa = kappa, phi = phi
  )
  dat <- simdat %>% group_by(time) %>% sample_n(N) %>% ungroup() # sub-sample from 'true' data
  spde <- make_spde(x, y, n_knots)

  m <- sdmTMB(
    silent = FALSE,
    ar1_fields = ar1_fields,
    include_spatial = TRUE,
    data = dat, formula = formula, time = "time",
    family = family, spde = spde
  )
  r <- m$tmb_obj$report()

  estimates <- c(
    ar1_phi = (2 * plogis(m$model$par[["ar1_phi"]]) - 1), # back transformed temporal correlation coefficent
    sigma_O = r$sigma_O,
    sigma_E = r$sigma_E, # spatio-temporal random field SD
    kappa = exp(r$ln_kappa), # Matern for space and space-time = distance at which ~%10 cor = sqrt(Type(8.0)) / exp(ln_kappa)
    phi = exp(r$ln_phi) # SD of observation error (aka sigma)
  )

  inputs <- c(
    ar1_phi = ar1_phi,
    sigma_O = sigma_O,
    sigma_E = sigma_E,
    kappa = kappa,
    phi = phi
  )
  parameter <- names(inputs)
  # diff <- estimates - inputs
  run <- data_frame(parameter = parameter, inputs = inputs, estimates = estimates, iter = iter)
  run
}

# model_sim()

arguments <- expand.grid(
  ar1_phi = 0.2, # c(-0.85, 0.1, 0.85),
  sigma_O = 0.3,
  sigma_E = 0.3,
  kappa = 0.1, # c(0.005, 0.1, 1),
  phi = 0.05, # c(0.01, 0.1),
  N = 100,
  n_knots = 200
)

nrow(arguments)
arguments$count <- 10L
arguments <- arguments[rep(seq_len(nrow(arguments)), arguments$count), ]
arguments_apply <- dplyr::select(arguments, -count)
nrow(arguments_apply)
arguments_apply$iter <- 1:nrow(arguments_apply)

par_all_iter <- purrr::pmap_dfr(arguments_apply,
  model_sim,
  grid = grid, x = grid$X, y = grid$Y,
  ar1_fields = TRUE, time_steps = 3, plot = TRUE,
  formula = z ~ 1, family = gaussian(link = "identity")
)


######################
######################

######################
# Add predictions for all points in simulated data
######################

######################
######################

sim_predictions <- function(iter = sample.int(1e3, 1), grid = grid, x = grid$X, y = grid$Y, time_steps = 3, plot = TRUE,
                            ar1_fields = TRUE,
                            ar1_phi = 0.5,
                            sigma_O = 0.3,
                            sigma_E = 0.3,
                            kappa = 0.05,
                            phi = 0.05,
                            N = 50, n_knots = 100, #iter.max=1e4, eval.max=1e4,
                            formula = z ~ 1, family = gaussian(link = "identity")) {
  set.seed(iter * 581267)
  simdat <- sim(
    x = x, y = y, time_steps = time_steps, plot = plot,
    ar1_fields = ar1_fields, ar1_phi = ar1_phi, sigma_O = sigma_O, sigma_E = sigma_E, kappa = kappa, phi = phi
  )

  dat <- simdat %>% group_by(time) %>% sample_n(N) %>% ungroup() # sub-sample from 'true' data
  spde <- make_spde(x, y, n_knots)
browser()
  m <- sdmTMB(
    silent = FALSE,
    ar1_fields = ar1_fields,
    include_spatial = TRUE,
    data = dat, formula = formula, time = "time",
    family = family, spde = spde
  )
  r <- m$tmb_obj$report()
  #browser()
  estimates <- c(
    ar1_phi = (2 * plogis(m$model$par[["ar1_phi"]]) - 1), # back transformed temporal correlation coefficent
    sigma_O = r$sigma_O,
    sigma_E = r$sigma_E, # spatio-temporal random field SD
    kappa = exp(r$ln_kappa), # Matern for space and space-time = distance at which ~%10 cor = sqrt(Type(8.0)) / exp(ln_kappa)
    phi = exp(r$ln_phi) # SD of observation error (aka sigma)
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
  nd <- do.call(
    "rbind",
    replicate(length(original_time), grid, simplify = FALSE)
  )
  nd[[m$time]] <- rep(original_time, each = nrow(grid))

  # run TMB with prediction turned on but replace sample 'dat' with new grid 'nd'
  predictions <- predict(m, newdata = nd, xy_cols = c("X", "Y"))$data
  head(predictions)

  # combine true and predicted values for each point in space and time
  spatial_bias <- full_join(simdat, predictions, by = c("x" = "X", "y" = "Y", "time" = "time"), suffix = c("", "_est")) %>% mutate(diff = est - real_z)
  spatial_bias$converg <- converg
  run <- list(par = data_frame(parameter = parameter, inputs = inputs, estimates = estimates, iter = iter, converg=converg), predicted = as_tibble(spatial_bias))
  run
}

#' Generate dataframe with parameter arguments for multiple simulation runs
#' Note: these parameters can be fixed to a single value, or vary using c(value1, value2, value3, ...)
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
generate_arg <- function(ar1_phi = 0.1, # c(-0.85, 0.1, 0.85),
                         sigma_O = 0.3,
                         sigma_E = 0.3,
                         kappa = 0.1, # c(0.005, 0.1, 1),
                         phi = 0.05, # c(0.01, 0.1),
                         time_steps = 3,
                         N = 2000,
                         n_knots = 200,
                         repeats = 2L) {
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

# arguments that can vary between runs are set using generate_arg to create a dataframe of all possible combinations
args <- generate_arg()
glimpse(args)

all_iter1 <- purrr::pmap(args, sim_predictions,
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

# result tibble of parameter estimates
params <- all_iter %>% map(~ .x[["par"]]) %>% bind_rows()

# OR if you want to use output from simulation of parameter values only
# params <- par_all_iter

# result list of tibbles of predictions or remove # and combine to one tibble
predictions <- all_iter %>% map(~ .x[["predicted"]]) # %>% bind_rows(.id = "iter")


# Calculate difference between parameter value input into simulation and estimate based on sdmTMB model
par_diff <- params %>% group_by(iter, parameter) %>% mutate(diff = inputs - estimates)
par_diff_converged <- par_diff %>% filter(converg==0) %>% ungroup()
par_diff_fail <- par_diff %>% filter(converg==1) %>% ungroup()

ggplot(par_diff, aes(x = diff, fill=as.factor(converg))) + geom_histogram() + facet_wrap(~parameter, scales = "free_x")

kappa <- par_diff %>% filter(parameter == "kappa")
ggplot(data = kappa, aes(x = estimates)) + geom_histogram(bins = 50) #+ xlim(-0.5,0.5)
ggplot(data = kappa, aes(x = diff)) + geom_histogram(bins = 50) #+ xlim(-0.5,0.5)

phi <- par_diff %>% filter(parameter == "phi")
ggplot(data = phi, aes(x = estimates)) + geom_histogram(bins = 50) #+ xlim(-0.5,0.5)

sigma_E <- par_diff %>% filter(parameter == "sigma_E")
ggplot(data = sigma_E, aes(x = diff)) + geom_histogram(bins = 50) #+ xlim(-0.5,0.5)

# make sure enough time steps are included if this is to be estimated?
ar1_phi <- par_diff %>% filter(parameter == "ar1_phi")
ggplot(data = ar1_phi, aes(x = diff)) + geom_histogram(bins = 50) #+ xlim(-0.5,0.5)

