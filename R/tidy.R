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
#' @param model Which model to tidy if a delta model (1 or 2).
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
#' @importFrom stats plogis
#' @examples
#' # See ?sdmTMB
tidy.sdmTMB <- function(x, effects = c("fixed", "ran_pars"), model = 1,
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

  .formula <- x$split_formula$fixedFormula
  .formula <- remove_s_and_t2(.formula)
  if (!"mgcv" %in% names(x)) x[["mgcv"]] <- FALSE
  fe_names <- colnames(model.matrix(.formula, x$data))

  se_rep <- as.list(x$sd_report, "Std. Error", report = TRUE)
  est_rep <- as.list(x$sd_report, "Estimate", report = TRUE)
  se <- as.list(x$sd_report, "Std. Error", report = FALSE)
  est <- as.list(x$sd_report, "Estimate", report = FALSE)
  b_j <- est$b_j[!fe_names == "offset"]
  b_j_se <- se$b_j[!fe_names == "offset"]
  fe_names <- fe_names[!fe_names == "offset"]
  out <- data.frame(term = fe_names, estimate = b_j[model], std.error = b_j_se[model], stringsAsFactors = FALSE)
  crit <- stats::qnorm(1 - (1 - conf.level) / 2)
  if (exponentiate) trans <- exp else trans <- I
  if (exponentiate) out$estimate <- trans(out$estimate)

  delta <- isTRUE(x$family$delta)
  assert_that(is.numeric(model))
  assert_that(length(model) == 1L)
  if (delta) assert_that(model %in% c(1, 2), msg = "`model` must be 1 or 2.")
  if (!delta) assert_that(model == 1, msg = "Only one model: `model` must be 1.")

  if (x$tmb_data$threshold_func > 0) {
    if (delta) stop("not implemented for threshold delta models yet.")
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
  # cleanup:
  est$epsilon_st <- NULL
  est$zeta_s <- NULL
  est$omega_s <- NULL
  est$ln_H_input <- NULL

  se$epsilon_st <- NULL
  se$zeta_s <- NULL
  se$omega_s <- NULL
  se$ln_H_input <- NULL

  subset_pars <- function(p, model) {
    p$b_j <- p$b_j[,model]
    p$ln_tau_O <- p$ln_tau_O[model]
    p$ln_tau_Z <- p$ln_tau_Z[model]
    p$ln_tau_E <- p$ln_tau_E[model]
    p$ln_kappa <- p$ln_kappa[,model]
    p$ln_phi <- p$ln_phi[model]
    p$ln_tau_V <- p$ln_tau_V[,model]
    p$ar1_phi <- p$ar1_phi[model]
    p$ln_tau_G <- p$ln_tau_G[,model]
    p$log_sigma_O <- p$log_sigma_O[model]
    p$log_sigma_E <- p$log_sigma_E[model]
    p$log_sigma_Z <- p$log_sigma_Z[model]
    p$log_range <- p$log_range[,model]
    p
  }
  est <- subset_pars(est, model)
  se <- subset_pars(se, model)

  if (x$family$family[[model]] %in% c("binomial", "poisson")) {
    se$ln_phi <- NULL
    est$ln_phi <- NULL
  }

  # not ADREPORTed for speed:
  optional_assign <- function(est, from, to) {
    if (from %in% names(est)) est[[to]] <- exp(est[[from]])
    est
  }
  est <- optional_assign(est, "log_sigma_E", "sigma_E")
  est <- optional_assign(est, "log_sigma_O", "sigma_O")
  est <- optional_assign(est, "log_sigma_Z", "sigma_Z")
  est <- optional_assign(est, "log_range", "range")
  est <- optional_assign(est, "ln_phi", "phi")
  est <- optional_assign(est, "ln_tau_G", "sigma_G")
  est <- optional_assign(est, "ln_tau_V", "sigma_V")

  ii <- 1
  if (length(unique(est$sigma_E)) == 1L) {
    se$sigma_E <- se$sigma_E[1]
    est$sigma_E <- est$sigma_E[1]
    se$log_sigma_E <- se$log_sigma_E[1]
    est$log_sigma_E <- est$log_sigma_E[1]
  }

  out_re <- list()
  log_name <- c("log_range")
  name <- c("range")
  if (!isTRUE(is.na(x$tmb_map$ln_phi))) {
    log_name <- c(log_name, "ln_phi")
    name <- c(name, "phi")
  }
  if (x$tmb_data$include_spatial) {
    log_name <- c(log_name, "log_sigma_O")
    name <- c(name, "sigma_O")
  }
  if (!x$tmb_data$spatial_only) {
    log_name <- c(log_name, "log_sigma_E")
    name <- c(name, "sigma_E")
  }
  if (x$tmb_data$spatial_covariate) {
    log_name <- c(log_name, "log_sigma_Z")
    name <- c(name, "sigma_Z")
  }
  if (x$tmb_data$random_walk) {
    log_name <- c(log_name, "ln_tau_V")
    name <- c(name, "tau_V")
  }
  if (length(est$ln_tau_G) > 0L) {
    log_name <- c(log_name, "ln_tau_G")
    name <- c(name, "sigma_G")
  }

  j <- 0
  if (!"log_range" %in% names(est)) {
    warning("This model was fit with an old version of sdmTMB. Some parameters may not be available to the tidy() method. Re-fit the model with the current version of sdmTMB if you need access to any missing parameters.", call. = FALSE)
  }
  for (i in name) {
    j <- j + 1
    if (i %in% names(est)) {
      .e <- est[[log_name[j]]]
      .se <- se[[log_name[j]]]
      .e <- if (is.null(.e)) NA else .e
      .se <- if (is.null(.se)) NA else .se
      out_re[[i]] <- data.frame(
        term = i, estimate = est[[i]], std.error = NA,
        conf.low = exp(.e - crit * .se),
        conf.high = exp(.e + crit * .se),
        stringsAsFactors = FALSE
      )

      ii <- ii + 1
    }
    out_re[[i]]$std.error <- NA
  }
  discard <- unlist(lapply(out_re, function(x) length(x) == 1L)) # e.g. old models and phi
  out_re[discard] <- NULL

  if ("tweedie" %in% x$family$family) {
    out_re$tweedie_p <- data.frame(
      term = "tweedie_p", estimate = plogis(est$thetaf) + 1,
      std.error = NA, stringsAsFactors = FALSE)
    out_re$tweedie_p$conf.low <- plogis(est$thetaf - crit * se$thetaf) + 1
    out_re$tweedie_p$conf.high <- plogis(est$thetaf + crit * se$thetaf) + 1
  }

  r <- x$tmb_obj$report()
  if (!is.null(r$rho) && r$rho != 0L) {
    ar_phi <- est$ar1_phi
    ar_phi_se <- se$ar1_phi
    rho_est <- 2 * stats::plogis(ar_phi) - 1
    rho_lwr <- 2 * stats::plogis(ar_phi - crit * ar_phi_se) - 1
    rho_upr <- 2 * stats::plogis(ar_phi + crit * ar_phi_se) - 1
    out_re[[ii]] <- data.frame(
      term = "rho", estimate = rho_est, std.error = NA,
      conf.low = rho_lwr, conf.high = rho_upr, stringsAsFactors = FALSE
    )
    ii <- ii + 1
  }

  if (x$control$map_rf) out_re$range <- NULL
  out_re <- do.call("rbind", out_re)
  row.names(out_re) <- NULL

  if (identical(est$ln_tau_E, 0)) out_re <- out_re[out_re$term != "sigma_E", ]
  if (identical(est$ln_tau_V, 0)) out_re <- out_re[out_re$term != "sigma_V", ]
  if (identical(est$ln_tau_G, 0)) out_re <- out_re[out_re$term != "sigma_G", ]
  if (identical(est$ln_tau_O, 0)) out_re <- out_re[out_re$term != "sigma_O", ]
  if (identical(est$ln_tau_Z, 0)) out_re <- out_re[out_re$term != "sigma_Z", ]

  if (!conf.int) {
    out_re[["conf.low"]] <- NULL
    out_re[["conf.high"]] <- NULL
  }

  out <- unique(out) # range can be duplicated
  out_re <- unique(out_re)

  if (effects == "fixed") {
    return(out)
  } else {
    return(out_re)
  }
}

#' @importFrom generics tidy
#' @export
generics::tidy
