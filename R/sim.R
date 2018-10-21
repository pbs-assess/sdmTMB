#' Simulate from a spatial Matern Random Field
#'
#' @param x A vector of x coordinates.
#' @param y A vector of y coordinates.
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
#' dat <- sim()
#' spde <- make_spde(x = dat$x, y = dat$y, n_knots = 150)
#' plot_spde(spde)
#' m <- sdmTMB(data = dat, formula = z ~ 1, time = "time",
#'   family = gaussian(link = "identity"), spde = spde)
#' r <- m$tmb_obj$report()
#' r$sigma_O
#' exp(r$ln_kappa)
#' exp(r$ln_phi)
#' s <- TMB::sdreport(m$tmb_obj)
#' head(summary(s))
sim <- function(x = runif(400, 0, 10), y = runif(400, 0, 10),
                sigma_O = 0.4, kappa = 1.3, phi = 0.2,
                seed = sample.int(1e6, 1), plot = FALSE) {
  set.seed(seed)
  rf_omega <- RandomFields::RMmatern(nu = 1, var = sigma_O^2, scale = 1 / kappa)
  omega_s <- suppressMessages(
    RandomFields::RFsimulate(model = rf_omega, x = x, y = y)$variable1)
  d <- data.frame(x, y)
  d$z <- stats::rnorm(nrow(d), mean = omega_s, sd = phi)
  if (plot) {
    g <- ggplot2::ggplot(d, ggplot2::aes_string("x", "y", colour = "z")) +
      ggplot2::geom_point() +
      ggplot2::scale_color_gradient2()
    print(g)
  }
  d$time <- 1
  d
}
