#' Extract parameter simulations from the joint precision matrix
#'
#' `spread_sims()` returns a wide-format data frame. `gather_sims()` returns a
#' long-format data frame. The format matches the format in the \pkg{tidybayes}
#' `spread_draws()` and `gather_draws()` functions.
#'
#' @param object Output from [sdmTMB()].
#' @param n_sims The number of simulation draws.
#'
#' @references
#' Code for simulating from the joint precision matrix adapted from:
#'
#' Osgood-Zimmerman, A. and Wakefield, J. 2021. A Statistical introduction
#' to Template Model Builder: a flexible tool for spatial modeling.
#' arXiv 2103.09929. <https://arxiv.org/abs/2103.09929>.
#'
#' @export
#' @rdname gather_sims
#' @examples
#' if (inla_installed()) {
#' m <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2,
#'   data = pcod_2011, spde = pcod_mesh_2011, family = tweedie(link = "log"),
#'   spatiotemporal = "AR1", time = "year")
#' head(spread_sims(m, n_sims = 10))
#' head(gather_sims(m, n_sims = 10))
#' samps <- gather_sims(m, n_sims = 1000)
#'
#' if (require("ggplot2", quietly = TRUE)) {
#'   ggplot(samps, aes(.value)) + geom_histogram() +
#'     facet_wrap(~.variable, scales = "free_x")
#' }
#' }

spread_sims <- function(object, n_sims = 200) {
  if (!"jointPrecision" %in% names(object$sd_report)) {
    stop("TMB::sdreport() must be run with the joint precision returned.", call. = FALSE)
  }
  tmb_sd <- object$sd_report
  samps <- rmvnorm_prec(object$tmb_obj$env$last.par.best, tmb_sd, n_sims)
  pars <- c(tmb_sd$par.fixed, tmb_sd$par.random)
  pn <- names(pars)
  pars <- pars[pn %in% names(tmb_sd$par.fixed)]
  samps <- samps[pn %in% names(tmb_sd$par.fixed), , drop = FALSE]
  pn <- pn[pn %in% names(tmb_sd$par.fixed)]
  .formula <- object$split_formula$fixedFormula
  if (isFALSE(object$mgcv)) {
    fe_names <- colnames(model.matrix(.formula, object$data))
  } else {
    fe_names <- colnames(model.matrix(mgcv::gam(.formula, data = object$data)))
  }
  fe_names <- fe_names[!fe_names == "offset"]
  row.names(samps) <- pn
  row.names(samps)[row.names(samps) == "b_j"] <- fe_names
  out <- as.data.frame(t(samps))
  if ("ln_kappa" %in% names(out)) {
    out$range <- sqrt(8) / exp(out$ln_kappa)
  }
  if ("ln_phi" %in% names(out)) {
    out$phi <- exp(out$ln_phi)
  }
  if ("thetaf" %in% names(out)) {
    out$tweedie_p <- stats::plogis(out$thetaf) + 1
  }
  if ("ar1_phi" %in% names(out)) {
    out$ar1_rho <- 2 * stats::plogis(out$ar1_phi) - 1
  }
  if ("ln_tau_O" %in% names(out)) {
    out$sigma_O <- 1 / sqrt(4 * pi * exp(2 * out$ln_tau_O + 2 * out$ln_kappa))
  }
  if ("ln_tau_E" %in% names(out)) {
    out$sigma_E <- 1 / sqrt(4 * pi * exp(2 * out$ln_tau_E + 2 * out$ln_kappa))
  }
  if ("ln_tau_O_trend" %in% names(out)) {
    out$sigma_O_trend <- 1 / sqrt(4 * pi * exp(2 * out$ln_tau_O_trend + 2 * out$ln_kappa))
  }
  out$ln_kappa <- out$ln_tau_O <- out$ln_tau_E <- out$ln_tau_O_trend <-
    out$ar1_phi <- out$thetaf <- out$ln_phi <- NULL
  data.frame(.iteration = seq_len(n_sims), out)
}

#' @export
#' @rdname gather_sims
gather_sims <- function(object, n_sims = 200) {
  out_wide <- spread_sims(object, n_sims)
  out_wide$.iteration <- NULL
  out <- stats::reshape(out_wide, direction = "long", varying = list(names(out_wide)),
    idvar = ".iteration", timevar = "variable_num")
  names(out)[2] <- ".value"
  row.names(out) <- NULL
  par_names <- data.frame(variable_num = unique(out$variable), .variable = names(out_wide))
  out <- base::merge(out, par_names)
  out$variable_num <- NULL
  out[ , c(".iteration", ".variable", ".value"), drop = FALSE]
}

rmvnorm_prec <- function(mu, tmb_sd, n_sims) {
  # mu <- c(tmb_sd$par.fixed, tmb_sd$par.random)
  z <- matrix(stats::rnorm(length(mu) * n_sims), ncol = n_sims)
  L <- Matrix::Cholesky(tmb_sd[["jointPrecision"]], super = TRUE)
  z <- Matrix::solve(L, z, system = "Lt")
  z <- Matrix::solve(L, z, system = "Pt")
  z <- as.matrix(z)
  mu + z
}
