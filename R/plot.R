#' Plot anisotropy
#'
#' @param object An object from [sdmTMB()].
#' @param arrow_length Arrow head length.
#'
#' @importFrom ggplot2 ggplot aes_string geom_segment coord_equal
#' @export
#' @rdname plot_anisotropy
#' @examples
#' \donttest{
#' d <- pcod
#' m <- sdmTMB(data = d,
#'   formula = density ~ 0 + as.factor(year),
#'   time = "year", spde = make_spde(d$X, d$Y, n_knots = 80),
#'   family = tweedie(link = "log"), anisotropy = TRUE,
#'   include_spatial = FALSE)
#' plot_anisotropy(m)
#' }
plot_anisotropy <- function(object, arrow_length = 10) {
  stopifnot(identical(class(object), "sdmTMB"))
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

