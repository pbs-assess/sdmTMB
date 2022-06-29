#' Sanity check of sdmTMB model
#'
#' @param fit Fitted model from [sdmTMB()]
#' @param se_ratio SE ratio to abs(parameter values) to issue warning
#' @param gradient_thresh Gradient threshold to issue warning
#'
#' @return An invisible named list of checks
#' @export
#'
#' @examplesIf inla_installed()
#' fit <- sdmTMB(
#'   present ~ s(depth),
#'   data = pcod_2011, mesh = pcod_mesh_2011,
#'   family = binomial()
#' )
#' sanity(fit)
#'
#' s <- sanity(fit)
#' s

sanity <- function(fit, se_ratio = 10, gradient_thresh = 0.001) {

  hessian_ok <- eigen_values_ok <- gradients_ok <- se_magnitude_ok <- FALSE
  nlminb_ok <- FALSE

  simplify_msg <- "Try simplifying the model, adjusting the mesh, or adding priors"

  if (identical(fit$model$convergence, 0L)) {
    msg <- "Non-linear minimizer suggests successful convergence"
    cli::cli_alert_success(msg)
    nlminb_ok <- TRUE
  } else {
    msg <- "Non-linear minimizer did not converge: do not trust this model!"
    cli::cli_alert_danger(msg)
    cli::cli_alert_info(simplify_msg)
    cat("\n")
  }

  if (isFALSE(fit$pos_def_hessian)) {
    msg <- "Non-positive-definite Hessian matrix: model may not have converged"
    cli::cli_alert_danger(msg)
    cli::cli_alert_info(simplify_msg)
    cat("\n")
  } else {
    msg <- "Hessian matrix is positive definite"
    cli::cli_alert_success(msg)
    hessian_ok <- TRUE
  }

  if (isTRUE(fit$bad_eig)) {
    msg <- "Extreme or very small eigen values detected: model may not have converged"
    cli::cli_alert_danger(msg)
    cli::cli_alert_info(simplify_msg)
    cat("\n")
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
      msg <- "See `?run_extra_optimization()`"
      cli::cli_alert_info(msg)
      msg <- "Or refit with `control = sdmTMBcontrol(newton_loops = 1)`"
      cli::cli_alert_info(msg)
      cat("\n")
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
  se_na_ok <- TRUE
  for (i in seq_along(se)) {
    if (is.na(se[i])) {
      cli::cli_alert_danger(c("`", np[i], paste0("` standard error is NA")))
      par_message(np[i])
      cli::cli_alert_info(simplify_msg)
      cat("\n")
      se_na_ok <- FALSE
    }
  }
  if (se_na_ok) {
    msg <- "No fixed-effect standard errors are NA"
    cli::cli_alert_success(msg)
  }

  est <- as.list(fit$sd_report, "Estimate")
  se <- as.list(fit$sd_report, "Std. Error")
  fixed <- !(names(est) %in% random)
  est <- est[fixed]
  se <- se[fixed]
  too_big <- function(est, se, divide = TRUE) {
    if (any(!is.na(se))) {
      if (divide) {
        ratio <- se[!is.na(se)] / abs(est[!is.na(se)])
        if (any(ratio > se_ratio)) return(TRUE)
      } else {
        se_max <- max(se, na.rm = TRUE)
        if (any(se_max > 3)) return(TRUE)
      }
    }
  }
  # log vars don't make a lot of sense to check like this:

  est <- as.list(fit$sd_report, "Estimate")
  se <- as.list(fit$sd_report, "Std. Error")

  bji <- grepl("^b_j", names(est))
  est <- est[bji]
  se <- se[bji]
  se_big <- mapply(too_big, est, se)

  # range and sigma pars:
  estr <- as.list(fit$sd_report, "Estimate", report = TRUE)
  ser <- as.list(fit$sd_report, "Std. Error", report = TRUE)
  se_big2 <- mapply(too_big, estr, ser, divide = FALSE)

  se_big <- c(se_big, se_big2)

  for (i in seq_along(se_big)) {
    if (isTRUE(se_big[[i]])) {
      msg <- paste0(
        "` standard error may be large")
      # (> ",
      #   se_ratio,
      #   "x parameter estimate)"
      # )
      cli::cli_alert_danger(c("`", names(se_big)[i], msg))
      par_message(names(se_big)[i])
      cli::cli_alert_info(simplify_msg)
      cat("\n")
    }
  }
  if (all(unlist(lapply(se_big, is.null)))) {
    msg <- "No fixed-effect standard errors look unreasonably large"
    cli::cli_alert_success(msg)
    se_magnitude_ok <- TRUE
  }

  b <- tidy(fit, conf.int = TRUE)
  br <- tidy(fit, "ran_pars", conf.int = TRUE)
  b <- rbind(b, br)

  if (isTRUE(fit$family$delta)) {
    b2 <- tidy(fit, conf.int = TRUE, model = 2)
    br2 <- tidy(fit, "ran_pars", conf.int = TRUE, model = 2)
    b2 <- rbind(b2, br2)
    b <- rbind(b, b2)
  }
  s <- grep("sigma", b$term)
  sigmas_ok <- TRUE
  if (length(s)) {
    for (i in s) {
      if (b$estimate[i] < 0.01) {
        msg <- "` is smaller than 0.01"
        cli::cli_alert_danger(c("`", b$term[i], msg))
        par_message(b$term[i])
        msg <- "Consider omitting this part of the model"
        cli::cli_alert_info(msg)
        cat("\n")
        sigmas_ok <- FALSE
      }
    }
  }
  if (sigmas_ok) {
    msg <- "No sigma parameters are < 0.01"
    cli::cli_alert_success(msg)
  }

  if (length(s)) {
    for (i in s) {
      if (b$estimate[i] > 100) {
        msg <- "` is larger than 100"
        cli::cli_alert_danger(c("`", b$term[i], msg))
        par_message(b$term[i])
        msg <- "Consider simplifying the model or adding priors"
        cli::cli_alert_info(msg)
        cat("\n")
        sigmas_ok <- FALSE
      }
    }
  }
  if (sigmas_ok) {
    msg <- "No sigma parameters are > 100"
    cli::cli_alert_success(msg)
  }

  r1 <- diff(range(fit$data[[fit$mesh$xy_cols[1]]]))
  r2 <- diff(range(fit$data[[fit$mesh$xy_cols[2]]]))
  r <- max(r1, r2)
  range_ok <- TRUE
  if ("range" %in% b$term) {
    if (max(b$estimate[b$term == "range"]) > r) {
      msg <- "A `range` parameter looks fairly large (> greatest distance in data)"
      cli::cli_alert_danger(msg)
      cli::cli_alert_info(simplify_msg)
      cli::cli_alert_info("Also make sure your spatial coordinates are not too big or small (e.g., work in UTM km instead of UTM m)", wrap = TRUE)
      cat("\n")
      range_ok <- FALSE
    } else {
      nr <- length(grep("range", b$term))
      if (nr == 1L) msg <- "Range parameter doesn't look unreasonably large"
      if (nr > 1L) msg <- "Range parameters don't look unreasonably large"
      cli::cli_alert_success(msg)
    }
  }

  ret <- named_list(
    hessian_ok, eigen_values_ok, nlminb_ok, range_ok,
    gradients_ok, se_magnitude_ok, se_na_ok, sigmas_ok
  )
  all_ok <- all(unlist(ret))
  ret <- c(ret, all_ok = all_ok)
  invisible(ret)
}

par_df <- function() {
  data.frame(
    internal = c("ln_tau_O", "ln_tau_E", "ln_tau_G", "ln_tau_V", "ln_tau_Z", "ln_kappa"),
    external = c("sigma_O", "sigma_E", "sigma_G", "sigma_V", "sigma_Z", "range"),
    meaning = c(
      "spatial standard deviation",
      "spatiotemporal standard deviation",
      "random intercept standard deviation",
      "time-varying coefficient standard deviation",
      "spatially varying coefficient standard deviation",
      "distance at which data are effectively independent"
      ),
    stringsAsFactors = FALSE
  )
}

par_message <- function(par) {
  df <- par_df()
  if (par %in% df$internal) {
    par_clean <- df$external[df$internal == par]
    meaning <- df$meaning[df$internal == par]
    cli::cli_alert_info(paste0("`", par, "` is an internal parameter affecting `", par_clean, "`"))
    cli::cli_alert_info(paste0("`", par_clean, "` is the ", meaning))
  }

  if (par %in% "bs") {
    msg <- "`bs` is an internal parameter presenting the linear compoment of a smoother"
    cli::cli_alert_info(msg)
    msg <- "It's not uncommon for this parameter to have large standard errors,\nwhich can possibly be ignored if the rest of the model looks OK"
    cli::cli_alert_info(msg)
  }
}

