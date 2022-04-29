#' Sanity check of sdmTMB model
#'
#' @param fit Fitted model from [sdmTMB()]
#' @param se_ratio SE ratio to abs(parameter values) to issue warning
#' @param gradient_thresh Gradient threshold to issue warning
#'
#' @return An invisible named list of checks
#' @export
#'
#' @examples
#' if (inla_installed()) {
#'   fit <- sdmTMB(
#'     present ~ s(depth),
#'     data = pcod_2011, mesh = pcod_mesh_2011,
#'     family = binomial()
#'   )
#'   sanity(fit)
#' }

sanity <- function(fit, se_ratio = 10, gradient_thresh = 0.01) {

  hessian_ok <- eigen_values_ok <- gradients_ok <- se_magnitude_ok <- FALSE
  nlminb_ok <- FALSE

  if (identical(fit$model$convergence, 0L)) {
    msg <- "Non-linear minimizer suggests successful convergence"
    cli::cli_alert_success(msg)
    nlminb_ok <- TRUE
  }

  simplify_msg <- "Try simplifying the model, adjusting the mesh, or adding priors"
  if (isFALSE(fit$pos_def_hessian)) {
    msg <- "Non-positive-definite Hessian matrix: model may not have converged"
    cli::cli_alert_danger(msg)
    cli::cli_alert_info(simplify_msg)
  } else {
    msg <- "Hessian matrix is positive definite"
    cli::cli_alert_success(msg)
    hessian_ok <- TRUE
  }

  if (isTRUE(fit$bad_eig)) {
    msg <- "Extreme or very small eigen values detected: model may not have converged"
    cli::cli_alert_danger(msg)
    cli::cli_alert_info(simplify_msg)
  } else {
    msg <- "No extreme or very small eigen values detected"
    cli::cli_alert_success(msg)
    eigen_values_ok <- TRUE
  }

  g <- fit$gradients
  np <- names(fit$model$par)
  for (i in seq_along(g)) {
    if (g[i] > gradient_thresh) {
      cli::cli_alert_danger(c(
        "`", np[i],
        paste0("` gradient > ", gradient_thresh)
      ))
      msg <- "See `?run_extra_optimization` or `?sdmTMBcontrol`"
      cli::cli_alert_info(msg)
    }
  }

  if (all(g <= gradient_thresh)) {
    msg <- "No gradients with respect to fixed effects are >= "
    cli::cli_alert_success(paste0(msg, gradient_thresh))
    gradients_ok <- TRUE
  }

  obj <- fit$tmb_obj
  random <- unique(names(obj$env$par[obj$env$random]))
  s <- summary(fit$sd_report)
  se <- s[,"Std. Error"]
  fixed_se <- !names(se) %in% random
  se <- se[fixed_se]
  np <- names(se)
  se_ok <- TRUE
  for (i in seq_along(se)) {
    if (is.na(se[i])) {
      cli::cli_alert_danger(c("`", np[i], paste0("` standard error is NA")))
      cli::cli_alert_info(simplify_msg)
      se_ok <- FALSE
    }
  }
  if (se_ok) {
    msg <- "No fixed-effect standard errors are NA"
    cli::cli_alert_success(msg)
  }

  est <- as.list(fit$sd_report, "Estimate")
  se <- as.list(fit$sd_report, "Std. Error")
  fixed <- !(names(est) %in% random)
  est <- est[fixed]
  se <- se[fixed]
  too_big <- function(est, se) {
    if (any(!is.na(se))) {
      ratio <- se[!is.na(se)] / abs(est[!is.na(se)])
      if (any(ratio > se_ratio)) TRUE
    }
  }
  se_big <- mapply(too_big, est, se)
  for (i in seq_along(se_big)) {
    if (isTRUE(se_big[[i]])) {
      msg <- "` standard error may be large (> 10x parameter estimate)"
      cli::cli_alert_danger(c("`", names(se_big)[i], msg))
      cli::cli_alert_info(simplify_msg)
    }
  }
  if (all(unlist(lapply(se_big, is.null)))) {
    msg <- "No fixed-effect standard errors look unreasonably large"
    cli::cli_alert_success(msg)
    se_magnitude_ok <- TRUE
  }

  b <- tidy(fit, conf.int = TRUE)
  b2 <- tidy(fit, "ran_pars", conf.int = TRUE)
  b <- rbind(b, b2)
  s <- grep("sigma", b$term)
  sigmas_ok <- TRUE
  if (length(s)) {
    for (i in s) {
      if (b$estimate[i] < 1e-3) {
        msg <- "` is smaller than 0.001"
        cli::cli_alert_danger(c("`", b$term[i], msg))
        msg <- "Consider omitting this part of the model"
        cli::cli_alert_info(msg)
        sigmas_ok <- FALSE
      }
    }
  }
  if (sigmas_ok) {
    msg <- "No sigma parameters are < 0.001"
    cli::cli_alert_success(msg)
  }

  ret <- named_list(
    hessian_ok, eigen_values_ok, nlminb_ok,
    gradients_ok, se_magnitude_ok, se_ok, sigmas_ok
  )
  all_ok <- all(unlist(ret))
  ret <- c(ret, all_ok = all_ok)
  invisible(ret)
}
