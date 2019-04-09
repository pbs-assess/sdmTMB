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
#' @param sigma_V Standard deviation of time-varying random walk on
#'   parameters. Set to 0 for parameters that should not vary through time.
#' @param seed A random seed.
#' @param plot Logical for whether or not to produce a plot.
#' @param list Logical for whether output is in list format:
#'    data in list element 1 and input values in element 2.
#'
#' @return A data frame where:
#'    `omega_s` represents the simulated spatial random effects.
#'    `epsilon_st` represents the simulated spatiotemporal random effects.
#'    `eta` is the estimate based on fixed effects for each point in space and time.
#'    `real` represents the simulated process without observation error.
#'    `observed` represents the simulated process with random observation error.
#'    `b*` contain the beta values for each covariate used to simulate each time slice.
#'    `cov*` covariate residuals for each observation.

#' @export
#'
#' @examples
#' set.seed(2957278)
#' dat <- sim(time_steps = 9, plot = TRUE)
#' spde <- make_spde(x = dat$x, y = dat$y, n_knots = 50)
#' plot_spde(spde)
#' m <- sdmTMB(
#'   data = dat, formula = observed ~ 1, time = "time",
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
#' spde <- make_spde(x = dat$x, y = dat$y, n_knots = 80)
#' m <- sdmTMB(
#'   silent = FALSE, ar1_fields = TRUE, include_spatial = FALSE,
#'   data = dat, formula = observed ~ 1, time = "time",
#'   family = gaussian(link = "identity"), spde = spde
#' )
#' r <- m$tmb_obj$report()
#' r$sigma_E
#' 2 * plogis(m$model$par[["ar1_phi"]]) - 1
#'
#' d <- sim(x = runif(90), y = runif(90), initial_betas = c(-0.2, 0.2),
#'   sigma_V = c(0.1, 0.1), time_steps = 10, phi = 0.1, ar1_fields = TRUE,
#'   ar1_phi = 0.5, plot = TRUE, sigma_O = 0.001, sigma_E = 0.3)
#' spde <- make_spde(d$x, d$y, n_knots = 50)
#' m <- sdmTMB(data = d, formula = observed ~ 1, time = "time",
#'   time_varying = ~ 0 + cov1 + cov2,
#'   silent = FALSE, ar1_fields = TRUE,
#'   include_spatial = FALSE, spde = spde)
#' m$model$par
#' r <- m$tmb_obj$report()
#' r$b_rw_t
#' unique(d[,c("b1", "b2")])
#' exp(r$ln_tau_V)

sim <- function(x = stats::runif(100, 0, 10),
                y = stats::runif(100, 0, 10),
                X = NULL,
                initial_betas = NULL,
                sigma_V = 0,
                time_steps = 1L,
                ar1_fields = FALSE,
                ar1_phi = 0.5,
                sigma_O = 0.01,
                sigma_E = 0.01,
                kappa = 0.01,
                phi = 0.01,
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
    b <- matrix(ncol = n_covariates, nrow = time_steps)
    b[1, ] <- initial_betas
    if (time_steps > 1) {
      for (k in seq_len(n_covariates)) {
        for (i in seq(2, time_steps)) {
          b[i, k] <- b[i - 1, k] + stats::rnorm(1, 0, sigma_V[k])
        }
      }
    }

    if (is.null(X)) {
      cov_mat <-  matrix(stats::rnorm(length(x) * n_covariates * time_steps, sd = 1),
        ncol = n_covariates, nrow = length(x) * time_steps)
      cov_mat <- as.data.frame(cov_mat)
      names(cov_mat) <- gsub("V", "cov", names(cov_mat))
    } else {
      cov_mat <- X
    }

    V <- b[rep(seq_len(time_steps), each = length(x)),]
    B <- as.data.frame(V)
    names(B) <- gsub("V", "b", names(B))
    cov_values <- B * cov_mat
    names(cov_values) <- gsub("b", "cov", names(cov_values))
    eta <- rowSums(cov_values)

  } else {
    # empty vector for eta in simulations with no fixed effects
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

  if (n_covariates > 0) {
    d <- cbind(d, B)
    d <- cbind(d, cov_mat)
  }

  if (plot) {
    g <- ggplot2::ggplot(d, ggplot2::aes_string("x", "y", colour = "observed")) +
      ggplot2::geom_point() +
      ggplot2::facet_wrap(~time) +
      ggplot2::scale_color_gradient2()
    print(g)
  }

  if (list) {
    sim_out <- list(data = d, inputs = list(
      ar1_phi = ar1_phi,
      sigma_O = sigma_O,
      sigma_E = sigma_E,
      kappa = kappa,
      phi = phi,
      initial_betas = initial_betas,
      sigma_V = sigma_V
    ))
    sim_out
  } else {
    d
  }
}

rf_sim <- function(model, x, y) {
  set.seed(sample.int(1e5L, 1L))
  suppressMessages(
    RandomFields::RFsimulate(model = model, x = x, y = y)$variable1
  )
}
