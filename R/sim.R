#' Simulate from a spatial Matern Random Field
#'
#' @param x A vector of x coordinates.
#' @param y A vector of y coordinates.
#' @param ar1_fields should random field draws be dependent on the previous year (TRUE) or not (FALSE).
#' @param ar1_phi correlation between years, must be between -1 and 1, default set to 0.5.
#' @param sigma_O SD of spatial process (Omega).
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
sim <- function(x = stats::runif(400, 0, 10), y = stats::runif(400, 0, 10),
                time_steps = 1L, ar1_fields = FALSE, ar1_phi = 0.5,
                sigma_O = 0.4, sigma_E = 0.3, kappa = 1.3, phi = 0.2,
                seed = sample.int(1e6, 1), plot = FALSE) {
  if (!identical(length(x), length((y)))) {
    stop("`x` and `y` must be of the same length.")
  }

  set.seed(seed)

  # spatial random effects: = omega
  rf_omega <- RandomFields::RMmatern(nu = 1, var = sigma_O^2, scale = 1 / kappa) # defines correlation matrix
  omega_s <- suppressMessages(
    RandomFields::RFsimulate(model = rf_omega, x = x, y = y)$variable1
  )

  # spatiotemporal random effects: = epsilon
  epsilon_st <- list()
  if (time_steps > 1L) {
    rf_epsilon <- RandomFields::RMmatern(nu = 1, var = sigma_E^2, scale = 1 / kappa)
    for (i in seq_len(time_steps)) {
      if (i == 1 || !ar1_fields) {
        epsilon_st[[i]] <- suppressMessages(
          RandomFields::RFsimulate(model = rf_epsilon, x = x, y = y)$variable1
        )
      } else { # AR1 and not first time slice:
        epsilon_st[[i]] <- ar1_phi * epsilon_st[[i - 1]] + suppressMessages(
          RandomFields::RFsimulate(model = rf_epsilon, x = x, y = y)$variable1
        )
      }
    }
  } else {
    epsilon_st <- list(rep(0, length(x)))
  }

  epsilon_st <- do.call("c", epsilon_st)


  d <- data.frame(x, y,
    time = rep(seq_len(time_steps), each = length(x)),
    omega_s = rep(omega_s, time_steps), epsilon_st = epsilon_st
  )
  d$real_z <- d$omega_s + d$epsilon_st
  d$z <- stats::rnorm(nrow(d), mean = d$real_z, sd = phi)

  if (plot) {
    g <- ggplot2::ggplot(d, ggplot2::aes_string("x", "y", colour = "z")) +
      ggplot2::geom_point() +
      ggplot2::facet_wrap(~time) +
      ggplot2::scale_color_gradient2()
    print(g)
  }
  d
}
