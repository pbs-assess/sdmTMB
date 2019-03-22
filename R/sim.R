#' Simulate from a spatial Matern Random Field
#'
#' @param x A vector of x coordinates.
#' @param y A vector of y coordinates.
#' @param time_steps The number of time steps.
#' @param ar1_fields Should random field draws be dependent on the previous year
#'   (`TRUE`) or not (`FALSE`).
#' @param ar1_phi Correlation between years; should be between -1 and 1.
#' @param sigma_O SD of spatial process (Omega).
#' @param sigma_E SD of spatiotemporal process (Epsilon).
#' @param kappa Parameter that controls the decay of spatial correlation.
#' @param phi Observation error scale parameter.
#' @param seed A random seed.
#' @param plot Logical for whether or not to produce a plot.
#'
#' @return A data frame. The column `z` represents the simulated process.
#' @export
#'
#' @examples
#' set.seed(2957278)
#' dat <- sim(time_steps = 9, plot = TRUE)
#' spde <- make_spde(x = dat$x, y = dat$y, n_knots = 200)
#' plot_spde(spde)
#' m <- sdmTMB(
#'   data = dat, formula = z ~ 1, time = "time",
#'   family = gaussian(link = "identity"), spde = spde
#' )
#' r <- m$tmb_obj$report()
#' r$sigma_O
#' r$sigma_E
#' s <- TMB::sdreport(m$tmb_obj)
#' head(summary(s))
#'
#' set.seed(1838)
#' dat <- sim(
#'   time_steps = 9, ar1_fields = TRUE, ar1_phi = 0.5,
#'   plot = TRUE, sigma_O = 0.01, sigma_E = 0.3, phi = 0.01
#' )
#' spde <- make_spde(x = dat$x, y = dat$y, n_knots = 120)
#' m <- sdmTMB(
#'   silent = FALSE, ar1_fields = TRUE, include_spatial = FALSE,
#'   data = dat, formula = z ~ 1, time = "time",
#'   family = gaussian(link = "identity"), spde = spde
#' )
#' r <- m$tmb_obj$report()
#' r$sigma_E
#' exp(r$ln_kappa)
#' 2 * plogis(m$model$par[["ar1_phi"]]) - 1
sim <- function(x = stats::runif(400, 0, 10),
                y = stats::runif(400, 0, 10),

                n_covariates = 0,
                timevarying = TRUE,
                B_initial = 0.2,
                year_sigma = 0.5,

                time_steps = 1L,
                ar1_fields = FALSE,
                ar1_phi = 0.5,
                sigma_O = 0.4,
                sigma_E = 0.3,
                kappa = 1.3,
                phi = 0.2,
                seed = sample.int(1e6, 1),
                plot = FALSE) {
  if (!identical(length(x), length((y)))) {
    stop("`x` and `y` must be of the same length.")
  }

  set.seed(seed)

  # spatial random effects: omega
  rf_omega <- RandomFields::RMmatern(nu = 1, var = sigma_O^2, scale = 1 / kappa)
  omega_s <- rf_sim(model = rf_omega, x, y)

  # spatiotemporal random effects: epsilon
  epsilon_st <- list()

  if (time_steps > 1L) {
    rf_epsilon <- RandomFields::RMmatern(nu = 1, var = sigma_E^2, scale = 1 / kappa)
    for (i in seq_len(time_steps)) {
      if (i == 1 || !ar1_fields) {
        epsilon_st[[i]] <- rf_sim(rf_epsilon, x, y)
      } else { # AR1 and not first time slice:
        epsilon_st[[i]] <- ar1_phi * epsilon_st[[i - 1]] + rf_sim(rf_epsilon, x, y)
      }
    }
  } else {
    epsilon_st <- list(rep(0, length(x)))
  }
  epsilon_st <- do.call("c", epsilon_st)


  # create betas for each covariate j
  if (n_covariates > 0) {
    B <- matrix(ncol = n_covariates, nrow = time_steps)
    B[1,] <- B_initial
    if (timevarying) {
      for (j in seq_len(n_covariates)) {
        for (i in 2:time_steps) {
          B[i,j] <- B[i-1,j] + rnorm(1, 0, year_sigma)
        }
      }
    }
  }

  # # create covariate matrix per time period
  # if (timevarying) {
  #   matrix_per_time <- list()
  #   for (i in 1:time_steps) {
  #     matrix_per_time[i] <- matrix(ncol = n_covariates, nrow = length(x))
  #     for (j in seq_len(n_covariates)) {
  #       matrix_per_time[[i]][[,j]] <- rnorm(length(x), 0, 1)
  #     }
  #   }
  # }
  #



 #


  # creat covariate matrix

  if (timevarying) {
    time_matrix <- model.matrix( ~ t - 1, data.frame( t = gl(time_steps, length(x))))
    matrix_per_covariate <- list()
    for (j in seq_len(n_covariates)) {
      matrix_per_covariate[[j]] <- matrix(1, ncol = time_steps, nrow = length(x)*time_steps)
      for (i in seq_len(time_steps)) {
        matrix_per_covariate[[j]][,i] <- rnorm(length(x)*time_steps, 0, 1)
      }
    }
    matrix_list <- list()
    for (j in seq_len(n_covariates)) {
      #Bj <- t(B[,j])
      matrix_list[[j]] <- t(apply(matrix_per_covariate[[j]],1,function (x) {x*t(B[,j])}))*time_matrix
      old_name <- colnames(matrix_list[[j]])
      colnames(matrix_list[[j]]) <- paste0("B", j, old_name)
      matrix_list[[j]] <- as_tibble(matrix_list[[j]])
    }
    model_matrix <- cbind(as_tibble(time_matrix),bind_cols(matrix_list))

  } else { # for none time-varying...
    # time_matrix <- model.matrix( ~ t - 1, data.frame( t = gl(time_steps, length(x))))
    # covariates <- list()
    # for (j in seq_len(n_covariates)) {
    #   matrix_per_covariate[[j]] <- matrix(1, ncol = time_steps, nrow = length(x)*time_steps)
    #   for (i in seq_len(time_steps)) {
    #     matrix_per_covariate[[j]][,i] <- rnorm(length(x)*time_steps, 0, 1)
    # }
    #   model_matrix <- cbind(as_tibble(time_matrix),bind_cols(matrix_list))

  }

# STILL NEED TO MERGE IN FIXED EFFECTS

  d <- data.frame(x, y, eta = as.vector(model_matrix, mode = "double"),
    time = rep(seq_len(time_steps), each = length(x)),
    omega_s = rep(omega_s, time_steps), epsilon_st = epsilon_st
  )
  d$real_z <- model_matrix + d$omega_s + d$epsilon_st

  # adds in observation error?
  d$z <- stats::rnorm(nrow(d), mean = d$real_z, sd = phi)

  if (plot) {
    g <- ggplot2::ggplot(d, ggplot2::aes_string("x", "y", colour = "z")) +
      ggplot2::geom_point() +
      ggplot2::facet_wrap(~time) +
      ggplot2::scale_color_gradient2()
    print(g)
  }
  d
  # list (d, inputs = c(ar1_phi=ar1_phi,
  # sigma_O=sigma_O, sigma_E=sigma_E, kappa=kappa, phi=phi)
}

## sim_args1 <- function(x = stats::runif(400, 0, 10), y = stats::runif(400, 0, 10),
##   time_steps = 1L, ar1_fields = FALSE, ar1_phi = 0.5,
##   sigma_O = 0.4, sigma_E = 0.3, kappa = 1.3, phi = 0.2,
##   seed = sample.int(1e6, 1), plot = FALSE) {
##
##   d = sim(x, y, time_steps, ar1_fields, ar1_phi, sigma_O, sigma_E, kappa, phi, seed, plot )
##
##   j = as.list(args("sim"))
##   j = j[!unlist(lapply(j,is.null))]
##   list(d,j)
## }
##
##
## # sim_args_list <- function(x = stats::runif(400, 0, 10), y = stats::runif(400, 0, 10),
## #   time_steps = 1L, ar1_fields = FALSE, ar1_phi = 0.5,
## #   sigma_O = 0.4, sigma_E = 0.3, kappa = 1.3, phi = 0.2,
## #   seed = sample.int(1e6, 1), plot = FALSE) {
## #
## #   d = sim(x  , y ,
## #     time_steps , ar1_fields , ar1_phi ,
## #     sigma_O , sigma_E , kappa , phi,
## #     seed , plot )
## #
## #   list(d, inputs = list(ar1_phi = ar1_phi,
## #                         sigma_O =sigma_O,
## #                         sigma_E = sigma_E,
## #                         kappa = kappa,
## #                         phi = phi,
## #                         seed = seed))
## # }
##
## sim_args_vec <- function(x = stats::runif(400, 0, 10), y = stats::runif(400, 0, 10),
##   time_steps = 1L, ar1_fields = FALSE, ar1_phi = 0.5,
##   sigma_O = 0.4, sigma_E = 0.3, kappa = 1.3, phi = 0.2,
##   seed = sample.int(1e6, 1), plot = FALSE) {
##
##   d = sim(x, y,
##     time_steps, ar1_fields, ar1_phi,
##     sigma_O, sigma_E, kappa , phi,
##     seed, plot )
##
##   list(d, inputs = c(ar1_phi = ar1_phi,
##     sigma_O =sigma_O,
##     sigma_E = sigma_E,
##     kappa = kappa,
##     phi = phi,
##     seed = seed))
## }

rf_sim <- function(model, x, y) {
  set.seed(sample.int(1e5L, 1L))
  suppressMessages(
    RandomFields::RFsimulate(model = model, x = x, y = y)$variable1
  )
}
