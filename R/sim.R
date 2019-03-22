#' Simulate from a spatial Matern Random Field
#'
#' @param x A vector of x coordinates.
#' @param y A vector of y coordinates.
#' @param X An optional covariate design matrix formatted as a list with each
#'   element of the list representing a slice in time.
#' @param time_steps The number of time steps.
#' @param ar1_fields Should random field draws be dependent on the previous year
#'   (`TRUE`) or not (`FALSE`).
#' @param ar1_phi Correlation between years; should be between -1 and 1.
#' @param sigma_O SD of spatial process (Omega).
#' @param sigma_E SD of spatiotemporal process (Epsilon).
#' @param kappa Parameter that controls the decay of spatial correlation.
#' @param phi Observation error scale parameter.
#' @param initial_betas Provide initial beta values, if model will include covariates.
#' @param year_sigma Standard deviation of time-varying random walk on
#'   parameters. Set to 0 for parameters that should not vary through time.
#' @param seed A random seed.
#' @param plot Logical for whether or not to produce a plot.
#' @param list Logical for whether output is in list format:
#'    data in list element [1] and input values in [2].
#'
#' @return A data fram36e. The column `z` represents the simulated process.
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
                X = NULL,
                initial_betas = NULL,
                year_sigma = 0,
                time_steps = 1L,
                ar1_fields = FALSE,
                ar1_phi = 0.5,
                sigma_O = 0.4,
                sigma_E = 0.3,
                kappa = 1.3,
                phi = 0.2,
                seed = sample.int(1e6, 1),
                plot = FALSE,
                list = FALSE) {
  if (!identical(length(x), length((y)))) {
    stop("`x` and `y` must be of the same length.")
  }

  set.seed(seed)
  n_covariates <- length(initial_betas)
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

  # create betas for each covariate k
  if (n_covariates > 0) {
    B <- matrix(ncol = n_covariates, nrow = time_steps)
    B[1, ] <- initial_betas
    if (time_steps > 1) {
      for (k in seq_len(n_covariates)) {
        for (i in seq(2, time_steps)) {
          B[i, k] <- B[i - 1, k] + stats::rnorm(1, 0, year_sigma[k])
        }
      }
    }

    # creat covariate matrix
    eta <- list()
    cov_mat <- list()
    for (i in seq_len(time_steps)) {
      if (is.null(X)) {
        cov_mat[[i]] <- matrix(stats::rnorm(length(x)),
          ncol = n_covariates, nrow = length(x)
        )
      } else {
        cov_mat[[i]] <- X[[i]]
      }
      eta[[i]] <- cov_mat[[i]] %*% B[i, ]
    }
    cov_mat <- do.call("rbind", cov_mat)
    eta <- do.call("c", eta)
  } else {
    eta <- vector("numeric", length = length(x) * time_steps)
  }

  d <- data.frame(
    time = rep(seq_len(time_steps), each = length(x)),
    x, y,
    omega_s = rep(omega_s, time_steps), epsilon_st = epsilon_st,
    eta = eta
  )
  d$real <- d$eta + d$omega_s + d$epsilon_st
  d$observed <- stats::rnorm(nrow(d), mean = d$real, sd = phi)

  B <- as.data.frame(B[rep(seq_len(time_steps), each = length(x)),])
  names(B) <- gsub("V", "b", names(B))
  d <- cbind(d, B)

  cov_mat <- as.data.frame(cov_mat)
  names(cov_mat) <- gsub("V", "cov_", names(cov_mat))
  d <- cbind(d, cov_mat)

  if (plot) {
    g <- ggplot2::ggplot(d, ggplot2::aes_string("x", "y", colour = "observed")) +
      ggplot2::geom_point() +
      ggplot2::facet_wrap(~time) +
      ggplot2::scale_color_gradient2()
    print(g)
  }

  if (list) {
    sim_out <- list(d, inputs = list(
      ar1_phi = ar1_phi,
      sigma_O = sigma_O,
      sigma_E = sigma_E,
      kappa = kappa,
      phi = phi,
      initial_betas = initial_betas,
      year_sigma = year_sigma
    ))
    sim_out
  } else {
    d
  }
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
