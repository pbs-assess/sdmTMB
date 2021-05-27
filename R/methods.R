#' @export
#' @import methods
print.sdmTMB <- function(x, ...) {
  r <- suppressWarnings(tryCatch(x$tmb_obj$report(), error = function(e) NA))
  if (all(is.na(r))) {
    stop("It looks like the model was built with a different version of sdmTMB. ",
      "Please fit your model again with the current sdmTMB.", call. = FALSE)
  }
  # need to initialize the new TMB object once:
  sink(tempfile())
  x$tmb_obj$fn(x$model$par)
  lp <- x$tmb_obj$env$last.par.best
  r <- x$tmb_obj$report(lp)
  sink()

  spatial_only <- !is.null(r$sigma_E) && !is.null(r$sigma_O_trend)

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
  time <- paste0("Time column: ", deparse(x$call$time), "\n")
  spde <- paste0("SPDE: ", deparse(x$call$spde), "\n")
  data <- paste0("Data: ", deparse(x$call$data), "\n")
  family <- paste0("Family: ", paste0(x$family$family, "(link = '", x$family$link, "')"), "\n")
  criterion <- paste0(fit_by, " criterion at convergence: ", mround(x$model$objective, 3), "\n")

  # .formula <- check_and_parse_thresh_params(x$formula, x$data)$formula
  .formula <- x$split_formula$fixedFormula
  if (!"mgcv" %in% names(x)) x[["mgcv"]] <- FALSE
  if (isFALSE(x$mgcv)) {
    fe_names <- colnames(model.matrix(.formula, x$data))
  } else {
    fe_names <- colnames(model.matrix(mgcv::gam(.formula, data = x$data)))
  }
  fe_names <- fe_names[!fe_names == "offset"]

  pars <- x$model$par
  b_j_exact <- unname(pars[grep("b_j", names(pars))])
  b_j <- round(b_j_exact, 2L)

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
  if (x$tmb_data$spatial_only == 0L) {
    if (!isTRUE(is.na(x$tmb_map$b_epsilon_logit))) {
      sigma_E <- paste0(pre, mround(r$sigma_E, 2L), "\n")
    } else {
      sigma_E <- paste0(pre, mround(r$sigma_E[1], 2L), "\n")
    }
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

  sr_se <- as.list(sr, "Std. Error")
  sr_est <- as.list(sr, "Estimate")

  if (x$tmb_data$threshold_func > 0) {
    mm_thresh <- cbind(sr_est$b_threshold, sr_se$b_threshold)
    if (x$threshold_function == 1L) {
      row.names(mm_thresh) <- paste0(x$threshold_parameter, c("-slope", "-breakpt"))
    } else {
      row.names(mm_thresh) <- paste0(x$threshold_parameter, c("-s50", "-s95", "-smax"))
    }
    colnames(mm_thresh) <- c("coef.est", "coef.se")
    mm_thresh[,1] <- round(mm_thresh[,1], 2)
    mm_thresh[,2] <- round(mm_thresh[,2], 2)

    mm <- rbind(mm, mm_thresh)
  } else {
    mm_thresh <- NULL
  }

  .tidy <- tidy(x, "ran_pars")
  if ("tau_G" %in% .tidy$term) {
    re_int_names <- barnames(x$split_formula$reTrmFormulas)
    re_int_mat <- matrix(NA_real_, nrow = length(re_int_names), ncol = 1)
    re_int_mat[,1] <- round(.tidy$estimate[.tidy$term == "tau_G"], 2)
    rownames(re_int_mat) <- paste(re_int_names, "(Intercept)")
    colnames(re_int_mat) <- "Std. Dev."
  } else {
    re_int_mat <- NULL
  }

  if (!is.null(x$time_varying)) {
    tv_names <- colnames(model.matrix(x$time_varying, x$data))
    mm_tv <- cbind(round(as.numeric(sr_est$b_rw_t), 2), round(as.numeric(sr_se$b_rw_t), 2))
    colnames(mm_tv) <- c("coef.est", "coef.se")
    time_slices <- sort(unique(x$data[[x$time]]))
    row.names(mm_tv) <- paste(rep(tv_names, each = length(time_slices)), time_slices, sep = "-")
  } else {
    mm_tv <- NULL
  }

  cat(title,
    formula,
    time,
    spde,
    data,
    family,
    sep = ""
  )

  print(mm)

  if (!is.null(re_int_mat)) {
    cat("\n")
    print(re_int_mat)
  }

  if (!is.null(x$time_varying)) {
    cat("\nTime-varying parameters:\n" )
    print(mm_tv)
  }

  cat("\n",
    paste0("Matern range parameter: ", range, "\n"),
    phi,
    sigma_O,
    sigma_E,
    rho,
    criterion,
    "\nSee ?tidy.sdmTMB to extract these values as a data frame.\n",
    sep = ""
  )

  invisible(list(fe = mm, tv = mm_tv, thresh = mm_thresh))
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
