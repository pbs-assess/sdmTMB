#' Simulate from a spatial Matern Random Field
#'
#' @param x A vector of x coordinates.
#' @param y A vector of y coordinates.
#' @param sigma_O SD of spatial process (Omega).
#' @param kappa Parameter that controls the decay of spatial correlation.
#' @param phi Observation error scale parameter.
#' @param seed A random seed.
#' @param plot Logical for whether or not to produce a plot.
#' @param B A vector of parameters. The first element is the intercept
#' @param phi The auto regressive parameter on the mean of the random field knots
#' @param X The model matrix
#
#' @return A data frame. The column `z` represents the simulated process.
#' @export
#'
#' @examples
#' set.seed(2957278)
#' dat <- simAR1(time_steps = 9, plot = TRUE)
#' spde <- make_spde(x = dat$x, y = dat$y, n_knots = 200)
#' plot_spde(spde)
#' m <- sdmTMB(data = dat, formula = z ~ 1, time = "time",
#'   family = gaussian(link = "identity"), spde = spde)
#' r <- m$tmb_obj$report()
#' r$sigma_O
#' r$sigma_E
#' s <- TMB::sdreport(m$tmb_obj)
#' head(summary(s))

simAR1 <- function(n_data_points = 400,
                time_steps = 1L,
                sigma_O = 0.4,
                sigma_E = 0.3,
                kappa = 1.3,
                phi = 0.2,
                B = c(0),                         # input fixed effects?
                X = rep(1, time_steps * n_data_points),     # input predictor matrix?
                x = stats::runif(n_data_points, 0, 10),
                y = stats::runif(n_data_points, 0, 10),
                seed = sample.int(1e6, 1),
                plot = FALSE
                ) {
  if (!identical(length(x), length((y))))
    stop("`x` and `y` must be of the same length.")

  set.seed(seed)

  # multiply coefficients by the design matrix
  eta <- as.vector(B %*% t(X), mode = "double")
          #eta_mat <- matrix(eta, nrow = nrow(proj), byrow = TRUE) # code from glmmfields that not sure i need

  # spatial random effects:
  rf_omega <- RandomFields::RMmatern(nu = 1, var = sigma_O^2, scale = 1 / kappa)
  omega_s <- suppressMessages(
    RandomFields::RFsimulate(model = rf_omega, x = x, y = y)$variable1)
          #browser()

  # spatiotemporal random effects:
  epsilon_st <- list()
  if (time_steps > 1L) {
    rf_epsilon <- RandomFields::RMmatern(nu = 1, var = sigma_E^2, scale = 1 / kappa)
    for (i in seq_len(time_steps)) {
      epsilon_st[[i]] <- suppressMessages(
        RandomFields::RFsimulate(model = rf_epsilon, x = x, y = y)$variable1)
    }
  } else {
    epsilon_st <- list(rep(0, length(x)))
  }

  epsilon_st <- do.call("c", epsilon_st)

  d <- data.frame(x, y, time = rep(seq_len(time_steps), each = length(n_data_points)),
    omega_s = rep(omega_s, time_steps), epsilon_st = epsilon_st)

  d$real_z <- d$omega_s + d$epsilon_st + eta

  d$z <- stats::rnorm(nrow(d), mean = d$real_z, sd = phi)

  # # potentially with AR process:
  # if (time_steps > 1) {
  #   for (i in seq(2, time_steps)) {
  #     if (mvt) {
  #       re_knots[i, ] <- mvtnorm::rmvt(1,
  #                                      # delta = ar * (re_knots[i - 1, ] - mean(re_knots[i - 1, ])),
  #                                      delta = phi * (re_knots[i - 1, ]),
  #                                      sigma = sigma_knots, df = df
  #       )
  #     }
  #     if (!mvt) {
  #       re_knots[i, ] <- mvtnorm::rmvnorm(1,
  #                                         # mean = ar * (re_knots[i - 1, ] - mean(re_knots[i - 1, ])),
  #                                         mean = phi * (re_knots[i - 1, ]),
  #                                         sigma = sigma_knots
  #       )
  #     }
  #   }
  # }

  if (plot) {
    gg <- ggplot2::ggplot(d, ggplot2::aes_string("x", "y", colour = "z")) +
      ggplot2::geom_point() +
      ggplot2::facet_wrap(~time) +
      ggplot2::scale_color_gradient2()
    print(gg)
  }
  d
}

