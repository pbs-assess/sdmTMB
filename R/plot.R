#' Plot anisotropy
#'
#' @param object An object from [sdmTMB()].
#'
#' @export
#' @rdname plot_anisotropy
#' @examples
#' \donttest{
#' if (inla_installed()) {
#' d <- pcod
#' m <- sdmTMB(data = d,
#'   formula = density ~ 0 + as.factor(year),
#'   time = "year", mesh = make_mesh(d, c("X", "Y"), n_knots = 80, type = "kmeans"),
#'   family = tweedie(link = "log"), anisotropy = TRUE,
#'   spatial = "off")
#' plot_anisotropy(m)
#' }
#' }
plot_anisotropy <- function(object) {
  stopifnot(identical(class(object), "sdmTMB"))
  report <- object$tmb_obj$report()
  eig <- eigen(report$H)
  dat <- data.frame(
    x0 = c(0, 0),
    y0 = c(0, 0),
    x1 = eig$vectors[1, , drop = TRUE] * eig$values,
    y1 = eig$vectors[2, , drop = TRUE] * eig$values
  )
  plot(0, xlim = range(c(dat$x0, dat$x1)),
    ylim = range(c(dat$y0, dat$y1)), type = "n", asp = 1, xlab = "", ylab = "")
  graphics::arrows(dat$x0, dat$y0, dat$x1, dat$y1)
  invisible(list(eig = eig, dat = dat, H = report$H))
}

