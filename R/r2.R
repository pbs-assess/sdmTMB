#' @importFrom performance r2
#' @export
performance::r2

#' Nakagawa R2
#'
#' @param x An [sdmTMB()] model object.
#'
#' @return A data frame of proportion variance explained.
#'
#' @export
r2.sdmTMB <- function(x) {
  if (!inherits(x, "sdmTMB")) {
    cli_abort("'x' must be a model of class sdmTMB.")
  }
  if (isTRUE(x$reml)) {
    cli_abort("r2.sdmTMB() does not yet work with REML")
  }
  if (length(x$family$family) > 1) {
    cli_abort("r2.sdmTMB() does not work for delta (hurdle) models yet.")
  }
  if (x$family$family == "student") {
    cli_inform("Family is student, but the variance does not (yet) account for the degrees of freedom.")
  }
  if (!x$family$family %in% c("student", "gaussian", "binomial", "tweedie", "Gamma", "poisson")) {
    cli_abort("r2.sdmTMB() currently only works for Gaussian, binomial, Gamma, Poisson, and Tweedie models.")
  }
  if (!is.null(x$spatial_varying)) {
    cli_abort("r2.sdmTMB() currently does not work with spatially varying coefficient models.")
  }

  varF <- .get_var_fe(x)
  varSmooths <- .get_var_sm(x)
  varO <- .get_var_sp(x)
  varE <- .get_var_spt(x)
  varV <- .get_var_tv(x)
  varG <- .get_var_re_iid(x)
  varR <- .get_var_resid(x)

  denominator <- varF + varSmooths + varO + varE + varR + varV + varG
  varF_all <- varF + varSmooths
  marginal <- varF_all / denominator
  cond_rf_sp <- cond_rf_spt <- cond_tv <- cond_re <- cond_all <- cond_smooth <- cond_fixed <- NULL
  if (varO != 0) {
    cond_rf_sp <- varO / denominator
  }
  if (varE != 0) {
    cond_rf_spt <- varE / denominator
  }
  if (varV != 0) {
    cond_tv <- varV / denominator
  }
  if (varG != 0) {
    cond_re <- varG / denominator
  }
  if (varSmooths != 0) {
    cond_smooth <- varSmooths / denominator
  }
  if (varF != 0) {
    cond_fixed <- varF / denominator
  }
  cond_all <- (denominator - varR) / denominator

  out <- list(
    conditional = cond_all,
    marginal = marginal,
    partial_smoothers = cond_smooth,
    partial_fixed = if (!is.null(cond_smooth)) cond_fixed else NULL,
    partial_spatial = cond_rf_sp,
    partial_spatiotemporal = cond_rf_spt,
    partial_time_varying = cond_tv,
    partial_random_intercepts = cond_re
  )
  out[vapply(out, is.null, logical(1L))] <- NULL
  ret <- t(as.data.frame(lapply(out, `[`, 1L)))
  out <- data.frame(component = row.names(ret), R2 = ret[, 1L, drop = TRUE], stringsAsFactors = FALSE)
  row.names(out) <- NULL
  if (requireNamespace("tibble", quietly = TRUE)) {
    out <- tibble::as_tibble(out)
  }
  out
}

.get_lp_vals <- function(x) {
  lp <- x$tmb_obj$env$last.par.best
  x$tmb_obj$report(lp)
}

.get_var_fe <- function(x) {
  r <- .get_lp_vals(x)
  var(r$eta_fixed_i[, 1L]) # FIXME delta
}

.get_var_sm <- function(x) {
  if (isTRUE(x$smoothers$has_smooths)) {
    r <- .get_lp_vals(x)
    varSmooths <- stats::var(r$eta_smooth_i)
  } else {
    varSmooths <- 0
  }
  varSmooths
}

.get_var_sp <- function(x) {
  if (x$tmb_data$include_spatial == 1L && !x$tmb_data$no_spatial) {
    b <- tidy(x, "ran_par")
    varO <- b$estimate[b$term == "sigma_O"]^2 # spatial variance
  } else {
    varO <- 0
  }
  varO
}

.get_var_spt <- function(x) {
  if (x$tmb_data$spatial_only == 0L && !x$tmb_data$no_spatial) {
    b <- tidy(x, "ran_par")
    varE <- b$estimate[b$term == "sigma_E"]^2 # spatiotemporal variance
  } else {
    varE <- 0
  }
  varE
}

.get_var_tv <- function(x) {
  if (x$tmb_data$random_walk %in% c(1L, 2L) || x$tmb_data$ar1_time) {
    if (!identical(x$time_varying, ~1)) {
      cli_abort("r2.sdmTMB() currently only works with time-varying intercepts.")
    }
    b <- tidy(x, "ran_par")
    varV <- b$estimate[b$term == "sigma_V"]^2 # time-varying variance
  } else {
    varV <- 0
  }
  varV
}

.get_var_re_iid <- function(x) {
  if (x$tmb_data$nobs_RE > 0) {
    b <- tidy(x, "ran_par")
    varG <- b$estimate[b$term == "sigma_G"]^2 # random effect variance
  } else {
    varG <- 0
  }
  varG
}

.get_variance_tweedie <- function(x, mu, phi) {
  p <- unname(stats::plogis(get_pars(x)$thetaf) + 1)
  phi * mu^p
}

.get_distribution_variance <- function(x) {
  phi <- exp(get_pars(x)$ln_phi)
  if (x$family$family %in% "gaussian") {
    return(phi^2)
  }
  if (x$family$family %in% "Gamma") {
    return(stats::family(x)$variance(exp(-0.5 * log(phi))))
  }
  if (x$family$family %in% c("tweedie", "Gamma", "poisson")) {
    re <- x$split_formula[[1]][[2]]
    if (!is.null(re)) {
      rterms <- paste0("(1 | ", re, ")") # FIXME: works for multiple!?
      nullform <- reformulate(rterms, response = ".")
    } else {
      nullform <- ". ~ 1"
    }
    null_model <- update(x, nullform)
    mu <- null_model$family$linkinv(unname(fixef(null_model)))
  }
  cvsquared <- tryCatch(
    {
      vv <- switch(x$family$family,
        poisson = family(x)$variance(mu),
        # Gamma = stats::family(x)$variance(phi),
        # nbinom1 = ,
        # nbinom2 = .variance_family_nbinom(x, mu, sig, faminfo),
        # truncated_nbinom2 = stats::family(x)$variance(mu, sig),
        tweedie = .get_variance_tweedie(x, mu, phi)
        # beta = .variance_family_beta(x, mu, sig),
      )
      if (vv < 0) {
        cli_warn("Model's distribution-specific variance is negative. Results are not reliable.")
      }
      vv / mu^2
    },
    error = function(x) {
      cli_warn("Can't calculate model's distribution-specific variance. Results are not reliable.")
      0
    }
  )
  log1p(cvsquared)
}

.get_sigma <- function(x) {
  exp(get_pars(x)$ln_phi)
}

.get_var_resid <- function(x) {
  if (x$family$family %in% c("student", "gaussian")) {
    varR <- .get_sigma(x)^2
  } else if (x$family$family == "binomial" && x$family$link == "logit") {
    # 'theoretical' method of Nakagawa supp. row 115
    varR <- pi^2 / 3
  } else if (x$family$family %in% c("tweedie", "Gamma", "poisson")) {
    varR <- .get_distribution_variance(x)
  } else {
    cli_abort("Family and/or link not implemented", call. = FALSE)
  }
  varR
}
