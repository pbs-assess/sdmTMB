#' Plot anisotropy
#'
#' @param object An object from [sdmTMB()].
#' @param model Which model if a delta model.
#'
#' @export
#' @rdname plot_anisotropy
#'
#' @return A plot of eigenvectors illustrating the estimated anisotropy. A list
#'   of the plotted data is invisibly returned.
#' @references Code adapted from VAST
#' @examples
#' \donttest{
#' if (inla_installed()) {
#'   d <- pcod
#'   m <- sdmTMB(
#'     data = d,
#'     formula = density ~ 0 + as.factor(year),
#'     time = "year", mesh = make_mesh(d, c("X", "Y"), n_knots = 80, type = "kmeans"),
#'     family = tweedie(link = "log"), anisotropy = TRUE,
#'     spatial = "off"
#'   )
#'   plot_anisotropy(m)
#' }
#' }
plot_anisotropy <- function(object, model = 1) {
  stopifnot(inherits(object, "sdmTMB"))
  report <- object$tmb_obj$report(object$tmb_obj$env$last.par.best)
  if (model == 1) eig <- eigen(report$H)
  if (model == 2) eig <- eigen(report$H2)
  dat <- data.frame(
    x0 = c(0, 0),
    y0 = c(0, 0),
    x1 = eig$vectors[1, , drop = TRUE] * eig$values,
    y1 = eig$vectors[2, , drop = TRUE] * eig$values
  )
  plot(0,
    xlim = range(c(dat$x0, dat$x1)),
    ylim = range(c(dat$y0, dat$y1)),
    type = "n", asp = 1, xlab = "", ylab = ""
  )
  graphics::arrows(dat$x0, dat$y0, dat$x1, dat$y1)
  invisible(list(eig = eig, dat = dat, H = report$H))
}

#' Plot a smooth term from an sdmTMB model
#'
#' @param object An [sdmTMB()] model.
#' @param select The smoother term to plot.
#' @param n The number of equally spaced points to evaluate the smoother along.
#' @param level The confidence level.
#' @param ggplot Logical: use the \pkg{ggplot2} package?
#' @param rug Logical: add rug lines along the lower axis?
#' @param return_data Logical: return the predicted data instead of making a plot?
#' @export
#'
#' @details
#' Note:
#' * Any numeric predictor is set to its mean
#' * Any factor predictor is set to its first-level value
#' * The time element (if present) is set to its minimum value
#' * The x and y coordinates are set to their mean values
#'
#' @examples
#' if (inla_installed()) {
#'   d <- subset(pcod, year >= 2000 & density > 0)
#'   pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 30)
#'   m <- sdmTMB(
#'     data = d,
#'     formula = log(density) ~ s(depth_scaled) + s(year, k = 5),
#'     mesh = pcod_spde
#'   )
#'   plot_smooth(m)
#' }
plot_smooth <- function(object, select = 1, n = 100, level = 0.95,
                        ggplot = FALSE, rug = TRUE, return_data = FALSE) {
  se <- TRUE
  if (isTRUE(object$delta))
    nice_stop("This function doesn't work with delta models yet")

  assert_that(inherits(object, "sdmTMB"))
  assert_that(is.logical(ggplot))
  assert_that(is.logical(return_data))
  assert_that(is.logical(se))
  assert_that(is.numeric(n))
  assert_that(is.numeric(level))
  assert_that(length(level) == 1L)
  assert_that(length(select) == 1L)
  assert_that(length(n) == 1L)
  assert_that(is.numeric(select))
  assert_that(level > 0 & level < 1)
  assert_that(n < 500)

  if (ggplot) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      nice_stop("ggplot2 not installed")
    }
  }

  sm <- parse_smoothers(object$formula[[1]], object$data)
  sm_names <- unlist(lapply(sm$Zs, function(x) attr(x, "s.label")))
  sm_names <- gsub("\\)$", "", gsub("s\\(", "", sm_names))

  fe_names <- colnames(object$tmb_data$X_ij)
  fe_names <- fe_names[!fe_names == "offset"] # FIXME
  fe_names <- fe_names[!fe_names == "(Intercept)"]

  all_names <- c(sm_names, fe_names)
  if (select > length(sm_names)) {
    nice_stop("`select` is greater than the number of smooths")
  }
  sel_name <- sm_names[select]
  non_select_names <- all_names[!all_names %in% sel_name]

  x <- object$data[[sel_name]]
  nd <- data.frame(x = seq(min(x), max(x), length.out = n))
  names(nd)[1] <- sel_name

  dat <- object$data
  .t <- terms(object$formula[[1]])
  .t <- labels(.t)
  checks <- c("^as\\.factor\\(", "^factor\\(")
  for (ch in checks) {
    if (any(grepl(ch, .t))) { # any factors from formula? if so, explicitly switch class
      ft <- grep(ch, .t)
      for (i in ft) {
        x <- gsub(ch, "", .t[i])
        x <- gsub("\\)$", "", x)
        dat[[x]] <- as.factor(dat[[x]])
      }
    }
  }
  dat[, object$spde$xy_cols] <- NULL
  dat[[object$time]] <- NULL
  for (i in seq_len(ncol(dat))) {
    if (names(dat)[i] != sel_name) {
      if (is.factor(dat[, i, drop = TRUE])) {
        nd[[names(dat)[[i]]]] <- sort(dat[, i, drop = TRUE])[[1]] # TODO note!
      } else {
        nd[[names(dat)[[i]]]] <- mean(dat[, i, drop = TRUE], na.rm = TRUE) # TODO note!
      }
    }
  }
  nd[object$time] <- min(object$data[[object$time]], na.rm = TRUE) # TODO note!
  nd[[object$spde$xy_cols[1]]] <- mean(object$data[[object$spde$xy_cols[1]]], na.rm = TRUE) # TODO note!
  nd[[object$spde$xy_cols[2]]] <- mean(object$data[[object$spde$xy_cols[2]]], na.rm = TRUE) # TODO note!

  p <- predict(object, newdata = nd, se_fit = se, re_form = NA)
  if (return_data) {
    return(p)
  }
  inv <- object$family$linkinv
  qv <- stats::qnorm(1 - (1 - level) / 2)

  if (!ggplot) {
    if (se) {
      lwr <- inv(p$est - qv * p$est_se)
      upr <- inv(p$est + qv * p$est_se)
      ylim <- range(c(lwr, upr))
    } else {
      ylim <- range(p$est)
    }
    plot(nd[[sel_name]], inv(p$est),
      type = "l", ylim = ylim,
      xlab = sel_name, ylab = paste0("s(", sel_name, ")")
    )
    if (se) {
      graphics::lines(nd[[sel_name]], lwr, lty = 2)
      graphics::lines(nd[[sel_name]], upr, lty = 2)
    }
    if (rug) rug(object$data[[sel_name]])
  } else {
    g <- ggplot2::ggplot(p, ggplot2::aes_string(sel_name, "inv(est)",
      ymin = "inv(est - qv * est_se)", ymax = "inv(est + qv * est_se)"
    )) +
      ggplot2::geom_line() +
      ggplot2::geom_ribbon(alpha = 0.4) +
      ggplot2::labs(x = sel_name, y = paste0("s(", sel_name, ")"))
    if (rug) {
      g <- g +
        ggplot2::geom_rug(
          data = object$data, mapping = ggplot2::aes_string(x = sel_name),
          sides = "b", inherit.aes = FALSE, alpha = 0.3
        )
    }
    return(g)
  }
}
