# dat <- sim(time_steps = 10, plot = TRUE, initial_betas = c(0.1, 0.2, -0.1), sigma_V = c(0,  0.2, 0.2))
# spde <- make_spde(x = dat$x, y = dat$y, n_knots = 50)
# plot_spde(spde)
# m <- sdmTMB(
#   data = dat, formula = observed ~ cov1, time = "time", include_spatial = T,
#   time_varying = ~ 0 + cov2 + cov2,
#   family = gaussian(link = "identity"), spde = spde
# )

#' @export
#' @import methods
print.sdmTMB <- function(x, ...) {
  r <- suppressWarnings(tryCatch(x$tmb_obj$report(), error = function(e) NA))
  if (all(is.na(r))) {
    x <- update_model(x)
    r <- x$tmb_obj$report()
  }

  spatial_only <- !is.null(r$r$sigma_E) && !is.null(r$r$sigma_O_trend)

  fit_by <- "ML"
  if ("reml" %in% names(x)) { # for backwards compatibility
    if (isTRUE(x$reml)) fit_by <- "REML" else "ML"
  }

  if (isTRUE(spatial_only)) {
    title <- paste0("Spatial model fit by ", fit_by, " ['sdmTMB']\n")
  } else {
    title <- paste0("Spatiotemporal model fit by ", fit_by, " ['sdmTMB']\n")
  }
  formula <- paste0("Formula: ", deparse(x$call$formula), "\n")
  # time <- paste0("Time column: ", deparse(x$call$time), "\n")
  spde <- paste0("SPDE: ", deparse(x$call$spde), "\n")
  data <- paste0("Data: ", deparse(x$call$data), "\n")
  family <- paste0("Family: ", paste0(x$family$family, "(link = '", x$family$link, "')"), "\n")
  criterion <- paste0(fit_by, " criterion at convergence: ", mround(x$model$objective, 3), "\n")

  if (!"mgcv" %in% names(x)) x[["mgcv"]] <- FALSE
  if (isFALSE(x$mgcv)) {
    fe_names <- colnames(model.matrix(x$formula, x$data))
  } else {
    fe_names <- colnames(model.matrix(mgcv::gam(x$formula, data = x$data)))
  }
  fe_names <- fe_names[!fe_names == "offset"]

  pars <- x$model$par
  b_j <- round(unname(pars[grep("b_j", names(pars))]), 2L)

  if ("ln_phi" %in% names(as.list(pars))) {
    phi <- mround(exp(as.list(pars)$ln_phi), 2L)
    phi <- paste0("Dispersion parameter: ", phi, "\n")
  } else {
    phi <- ""
  }
  range <- mround(r$range, 2L)

  pre <- "Spatial SD (sigma_O): "
  if (!is.null(r$sigma_O)) {
    sigma_O <- paste0(pre, mround(r$sigma_O, 2L), "\n")
  } else {
    sigma_O <- ""
  }

  pre <- "Spatiotemporal SD (sigma_E): "
  if (x$tmb_data$spatial_only == 0) {
    sigma_E <- paste0(pre, mround(r$sigma_E, 2L), "\n")
  } else {
    sigma_E <- paste0(pre, "not estimated\n")
  }

  pre <- "Spatiotemporal AR1 correlation (rho): "
  if (!is.null(r$rho) && r$rho != 0L) {
    rho <- paste0(pre, mround(r$rho, 2L), "\n")
  } else {
    rho <- ""
  }

  sr <- x$sd_report
  sr_se <- summary(sr)[, "Std. Error"]
  sr_est <- summary(sr)[, "Estimate"]
  b_j_se <- unname(round(sr_se[grep("b_j", names(sr_se))], 2L))
  b_j <- unname(round(sr_est[grep("b_j", names(sr_est))], 2L))

  mm <- cbind(b_j, b_j_se)
  colnames(mm) <- c("coef.est", "coef.se")
  row.names(mm) <- fe_names

  # Add pretty-printing of threshold parameters FIXME

  cat(title,
    formula,
    # time,
    spde,
    # data,
    family,
    sep = ""
  )

  print(mm)

  cat("\n",
    paste0("Matern range parameter: ", range, "\n"),
    phi,
    sigma_O,
    sigma_E,
    rho,
    criterion,
    sep = ""
  )
}

#' @export
summary.sdmTMB <- function(object, ..., digits) {
  print(object, ...)
}

mround <- function(x, digits) {
  sprintf(paste0("%.", digits, "f"), round(x, digits))
}

#' Extract the number of observations of a sdmTMB model
#'
#' @param object The fitted sdmTMB model object
#' @importFrom stats nobs
#' @export
#' @noRd
nobs.sdmTMB <- function(object, ...)
  sum(!is.na(object$data[all.vars(object$formula)[1]]))

#' Extract the log likelihood of a sdmTMB model
#'
#' @param object The fitted sdmTMB model object
#' @importFrom stats logLik
#' @export
#' @noRd
logLik.sdmTMB <- function(object, ...) {
  val <- -object$model$objective

  nobs <- nobs.sdmTMB(object)
  df <- length(object$model$par) # fixed effects only
  structure(val,
    nobs = nobs, nall = nobs, df = df,
    class = "logLik"
  )
}

#' Extract the AIC of a sdmTMB model
#'
#' @param fit The fitted sdmTMB model
#' @param scale The scale (note used)
#' @param k Penalization parameter, defaults to 2
#' @param ... Anything else
#' @noRd
#'
#' @export
extractAIC.sdmTMB <- function(fit, scale, k = 2, ...) {
  L <- logLik(fit)
  edf <- attr(L, "df")
  return(c(edf, c(-2 * L + k * edf)))
}
