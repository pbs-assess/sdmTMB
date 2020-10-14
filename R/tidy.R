#' Turn sdmTMB model output into a tidy data frame
#'
#' @param x Output from [sdmTMB()].
#' @param effects A character vector including one or more of "fixed"
#'   (fixed-effect parameters), "ran_pars" (standard deviations, spatial range,
#'   and other random effect terms).
#' @param conf.int Include a confidence interval?
#' @param conf.level Confidence level for CI.
#' @param exponentiate Whether to exponentiate the fixed-effect coefficient
#'   estimates and confidence intervals.
#' @param ... Extra arguments (not used).
#'
#' @return A data frame
#' @details
#' Follows the conventions of the \pkg{broom} and \pkg{broom.mixed} packages.
#'
#' Note that the standard errors for variance terms are not in natural
#' and not log space and so are a very rough approximation.
#' @export
#'
#' @importFrom assertthat assert_that
#' @examples
#' # See ?sdmTMB
tidy.sdmTMB <- function(x, effects = c("fixed", "ran_pars"),
                 conf.int = FALSE, conf.level = 0.95, exponentiate = FALSE, ...) {
  effects <- match.arg(effects)
  assert_that(is.logical(exponentiate))
  assert_that(is.logical(conf.int))
  if (conf.int) {
    assert_that(is.numeric(conf.level),
      conf.level > 0, conf.level < 1,
      length(conf.level) == 1,
      msg = "`conf.level` must be length 1 and between 0 and 1")
  }

  .formula <- check_and_parse_thresh_params(x$formula, x$data)$formula
  if (!"mgcv" %in% names(x)) x[["mgcv"]] <- FALSE
  if (isFALSE(x$mgcv)) {
    fe_names <- colnames(model.matrix(.formula, x$data))
  } else {
    fe_names <- colnames(model.matrix(mgcv::gam(.formula, data = x$data)))
  }
  fe_names <- fe_names[!fe_names == "offset"]
  se_rep <- as.list(x$sd_report, "Std. Error", report = TRUE)
  est_rep <- as.list(x$sd_report, "Estimate", report = TRUE)
  se <- as.list(x$sd_report, "Std. Error", report = FALSE)
  est <- as.list(x$sd_report, "Estimate", report = FALSE)
  b_j <- est$b_j
  b_j_se <- se$b_j
  out <- data.frame(term = fe_names, estimate = b_j, std.error = b_j_se, stringsAsFactors = FALSE)
  crit <- stats::qnorm(1 - (1 - conf.level) / 2)
  if (exponentiate) trans <- exp else trans <- I

  if (x$tmb_data$threshold_func > 0) {
    if (x$threshold_function == 1L) {
      par_name <- paste0(x$threshold_parameter, c("-slope", "-breakpt"))
    } else {
      par_name <- paste0(x$threshold_parameter, c("-s50", "-s95", "-smax"))
    }
    out <- rbind(
      out,
      data.frame(
        term = par_name, estimate = est$b_threshold,
        std.error = se$b_threshold, stringsAsFactors = FALSE
      )
    )
  }

  if (conf.int) {
    out$conf.low <- as.numeric(trans(out$estimate - crit * out$std.error))
    out$conf.high <- as.numeric(trans(out$estimate + crit * out$std.error))
  }


  se <- c(se, se_rep)
  est <- c(est, est_rep)
  ii <- 1

  out_re <- list()
  for (i in c("sigma_O", "sigma_E", "sigma_O_trend", "ln_tau_V", "range", "ln_phi")) {
    if (i %in% names(est) && est[[i]] != 0) {
      out_re[[i]] <- data.frame(
        term = i, estimate = est[[i]], std.error = se[[i]],
        conf.low = NA_real_, conf.high = NA_real_, stringsAsFactors = FALSE
      )
      if (i == "sigma_O_trend") out_re[[i]]$term <- "sigma_Z"
      ii <- ii + 1
    }
  }

  r <- x$tmb_obj$report()
  if (!is.null(r$rho) && r$rho != 0L) {
    ar_phi <- est$ar1_phi
    ar_phi_se <- se$ar1_phi
    ar_phi_est <- 2 * stats::plogis(ar_phi) - 1
    ar_phi_lwr <- 2 * stats::plogis(ar_phi - crit * ar_phi_se) - 1
    ar_phi_upr <- 2 * stats::plogis(ar_phi + crit * ar_phi_se) - 1
    out_re[[ii]] <- data.frame(
      term = "ar1_phi", estimate = ar_phi_est, std.error = NA_real_,
      conf.low = ar_phi_lwr, conf.high = ar_phi_upr, stringsAsFactors = FALSE
    )
    ii <- ii + 1
  }

  out_re <- do.call("rbind", out_re)
  row.names(out_re) <- NULL

  if (!conf.int) {
    out_re[["conf.low"]] <- NULL
    out_re[["conf.high"]] <- NULL
  }

  if (effects == "fixed") {
    return(out)
  } else {
    return(out_re)
  }
}

#' @importFrom generics tidy
#' @export
generics::tidy
