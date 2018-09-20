#' Plot anisotropy
#'
#' @param object An object from [sdmTMB()].
#' @param arrow_length Arrow head length.
#'
#' @importFrom ggplot2 ggplot aes_string geom_segment coord_equal
#' @export
#' @rdname plot_anisotropy
#' @examples
#' set.seed(1) # affects INLA mesh
#' m <- sdmTMB(subset(pcod, year >= 2013),
#'   density ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
#'   time = "year", spde = make_spde(pcod$X, pcod$Y, n_knots = 100),
#'   family = tweedie(link = "log"), anisotropy = TRUE)
#' plot_anisotropy(m)
plot_anisotropy <- function(object, arrow_length = 10) {
  assertthat::are_equal(class(object), "sdmTMB")
  report <- object$tmb_obj$report()
  eig = eigen(report$H)
  dat <- data.frame(
    x0 = c(0, 0),
    y0 = c(0, 0),
    x1 = eig$vectors[1, , drop = TRUE] * eig$values,
    y1 = eig$vectors[2, , drop = TRUE] * eig$values
  )
  ggplot(dat, aes_string(x = "x0", y = "y0", xend = "x1", yend = "y1")) +
    geom_segment(arrow = grid::arrow(length = grid::unit(arrow_length, "points"))) +
    coord_equal()
}

